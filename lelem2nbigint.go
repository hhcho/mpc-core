package mpc_core

import (
	"math/big"

	"github.com/hhcho/frand"
)

type LElem2NBigInt struct {
	Val *big.Int
}

const LElem2NBigIntUniqueID uint8 = 7

var LElem2NBigIntZero LElem2NBigInt = LElem2NBigInt{big.NewInt(0)}
var LElem2NBigIntMod LElem2NBigInt     // Variable set on the fly
var LElem2NBigIntModHalf LElem2NBigInt // Variable set on the fly
var LElem2NBigIntModBitLen int         // Variable set on the fly

/* LElem2NBigInt */

func (a LElem2NBigInt) SetModulusPowerOf2(k uint) {
	LElem2NBigIntMod = LElem2NBigInt{big.NewInt(0).Lsh(big.NewInt(1), k)}
	LElem2NBigIntModHalf = LElem2NBigInt{big.NewInt(0).Lsh(big.NewInt(1), k-1)}
	LElem2NBigIntModBitLen = int(k + 1)
}

func (a LElem2NBigInt) Mul(b interface{}) RElem {
	m := new(big.Int).Mul(a.Val, b.(LElem2NBigInt).Val)
	return LElem2NBigInt{new(big.Int).Mod(m, LElem2NBigIntMod.Val)}
}

func (a LElem2NBigInt) Add(b interface{}) RElem {
	m := new(big.Int).Add(a.Val, b.(LElem2NBigInt).Val)
	return LElem2NBigInt{new(big.Int).Mod(m, LElem2NBigIntMod.Val)}
}

func (a LElem2NBigInt) Sub(b interface{}) RElem {
	return a.Add(b.(LElem2NBigInt).Neg())
}

func (a LElem2NBigInt) Neg() RElem {
	m := new(big.Int).Neg(a.Val)
	return LElem2NBigInt{new(big.Int).Mod(m, LElem2NBigIntMod.Val)}
}

func (a LElem2NBigInt) Inv() RElem {
	bInv := big.NewInt(0).ModInverse(a.Val, LElem2NBigIntMod.Val)
	if bInv == nil {
		panic("ModInverse does not exist")
	}
	return LElem2NBigInt{bInv}
}

func (a LElem2NBigInt) ToBigInt() *big.Int {
	return new(big.Int).Set(a.Val)
}

func (a LElem2NBigInt) FromBigInt(n *big.Int) RElem {
	return LElem2NBigInt{new(big.Int).Set(n)}
}

func (a LElem2NBigInt) ToBigFloat(fracBits int) *big.Float {
	out := new(big.Float).SetInt(a.Val)
	out.Mul(out, new(big.Float).SetMantExp(big.NewFloat(1), -fracBits))
	return out
}

func (a LElem2NBigInt) ToSignedBigInt() *big.Int {
	if a.Val.Cmp(LElem2NBigIntModHalf.Val) >= 0 {
		return new(big.Int).Sub(a.Val, LElem2NBigIntMod.Val)
	}
	return new(big.Int).Set(a.Val)
}

func (a LElem2NBigInt) ToSignedBigFloat(fracBits int) *big.Float {
	out := new(big.Float).SetInt(a.ToSignedBigInt())
	out.Mul(out, new(big.Float).SetMantExp(big.NewFloat(1), -fracBits))
	return out
}

func (a LElem2NBigInt) FromBigFloat(n *big.Float, fracBits int) RElem {
	f := new(big.Float).Mul(n, new(big.Float).SetMantExp(big.NewFloat(1), fracBits))
	return a.FromBigInt(roundFloat(f))
}

func (a LElem2NBigInt) AssertTypeFor(n RElem) RElem {
	return n.(LElem2NBigInt)
}
func (a LElem2NBigInt) FromInt(n int) RElem {
	if n >= 0 {
		return LElem2NBigInt{new(big.Int).SetInt64(int64(n))}
	}
	return LElem2NBigInt{new(big.Int).Sub(LElem2NBigIntMod.Val, new(big.Int).SetInt64(-int64(n)))}
}

func (a LElem2NBigInt) FromUint64(n uint64) RElem {
	return LElem2NBigInt{new(big.Int).SetUint64(n)}
}

func (a LElem2NBigInt) FromFloat64(n float64, fracBits int) RElem {
	if n < 0 {
		return LElem2NBigInt{LElem2NBigIntZero.FromFloat64(-n, fracBits).(LElem2NBigInt).ToBigInt()}.Neg()
	} else {
		return LElem2NBigInt{LElem2NBigIntZero.FromFloat64(n, fracBits).(LElem2NBigInt).ToBigInt()}
	}
}

func (a LElem2NBigInt) Float64(fracBits int) float64 {
	var sgn int
	var b LElem2NBigInt
	if a.Val.Cmp(LElem2NBigIntModHalf.Val) < 0 {
		sgn, b = 1, a
	} else {
		sgn, b = -1, a.Neg().(LElem2NBigInt)
	}

	out, _ := new(big.Rat).SetFrac(b.Val, new(big.Int).Lsh(big.NewInt(1), uint(fracBits))).Float64()
	return float64(sgn) * out
}

func (a LElem2NBigInt) FromBytes(buf []byte) RElem {
	return LElem2NBigInt{new(big.Int).SetBytes(buf[:a.NumBytes()])}
}

func (a LElem2NBigInt) ToBytes(buf []byte) {
	if a.Val.Sign() < 0 {
		a.Val.Mod(a.Val, LElem2NBigIntMod.Val)
	}
	a.Val.FillBytes(buf[:a.NumBytes()])
}

func (a LElem2NBigInt) Uint64() uint64 {
	return a.Val.Uint64()
}
func (a LElem2NBigInt) Zero() RElem {
	return LElem2NBigInt{big.NewInt(0)}
}

func (a LElem2NBigInt) One() RElem {
	return LElem2NBigInt{big.NewInt(1)}
}

func (a LElem2NBigInt) Modulus() *big.Int {
	return new(big.Int).Set(LElem2NBigIntMod.Val)
}

func (a LElem2NBigInt) NumBytes() uint32 {
	return 1 + ((uint32(LElem2NBigIntModBitLen) - 2) >> 3)
}

func (a LElem2NBigInt) ModBitLength() int {
	return LElem2NBigIntModBitLen
}

func (a LElem2NBigInt) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB >= a.ModBitLength() {
		panic("Invalid bit position")
	}
	return a.Val.Bit(posFromLSB)
}
func (a LElem2NBigInt) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return LElem2NBigInt{new(big.Int).Set(a.Val)}
	}

	return LElem256Big{new(big.Int).Rem(a.Val, new(big.Int).Lsh(big.NewInt(1), uint(nBits)))}
}

func (a LElem2NBigInt) TypeID() uint8 {
	return LElem2NBigIntUniqueID
}

func (a LElem2NBigInt) Rand(prg *frand.RNG) RElem {
	bitLen := a.ModBitLength() - 1
	numBytes := 1 + ((bitLen - 1) >> 3)
	buf := make([]byte, numBytes)
	headLen := bitLen - ((numBytes - 1) << 3)

	prg.Read(buf)
	if headLen < 8 { // SetBytes is big-endian
		buf[0] = uint8(buf[0]) % (1 << headLen)
	}

	return LElem2NBigInt{new(big.Int).SetBytes(buf)}
}

func (a LElem2NBigInt) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits >= a.ModBitLength() {
		panic("Requested bit length is larger than modulus")
	}

	buf := make([]byte, 1+((nbits-1)/8))
	prg.Read(buf)

	return LElem2NBigInt{new(big.Int).Rem(new(big.Int).SetBytes(buf), new(big.Int).Lsh(big.NewInt(1), uint(nbits)))}
}

func (a LElem2NBigInt) Copy() RElem {
	return LElem2NBigInt{new(big.Int).Set(a.Val)}
}
