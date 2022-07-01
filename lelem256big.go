package mpc_core

import (
	"fmt"
	"math/big"

	"github.com/hhcho/frand"
)

type LElem256Big struct {
	Val *big.Int
}

const LElem256BigBytes uint32 = 32 // Extra byte for sign bit for ease of marshalling
var LElem256BigZero LElem256Big = LElem256Big{big.NewInt(0)}
var LElem256BigMod = LElem256Big{new(big.Int).Sub(big.NewInt(0).Lsh(big.NewInt(1), 256), big.NewInt(189))}    // 2^256 - 189
var LElem256BigModHalf = LElem256Big{new(big.Int).Sub(big.NewInt(0).Lsh(big.NewInt(1), 255), big.NewInt(94))} // 2^255 - 94
var LElem256BigModBitLen int = (LElem256BigMod).Val.BitLen()

const LElem256BigUniqueID uint8 = 8

/* LElem256Big */

func (a LElem256Big) Mul(b interface{}) RElem {
	m := new(big.Int).Mul(a.Val, b.(LElem256Big).Val)
	fmt.Println("m", m)
	return LElem256Big{new(big.Int).Mod(m, LElem256BigMod.Val)}
}

func (a LElem256Big) Add(b interface{}) RElem {
	m := new(big.Int).Add(a.Val, b.(LElem256Big).Val)
	return LElem256Big{new(big.Int).Mod(m, LElem256BigMod.Val)}
}

func (a LElem256Big) Sub(b interface{}) RElem {
	return a.Add(b.(LElem256Big).Neg())
}

func (a LElem256Big) Neg() RElem {
	m := new(big.Int).Neg(a.Val)
	return LElem256Big{new(big.Int).Mod(m, LElem256BigMod.Val)}
}

func (a LElem256Big) Inv() RElem {
	bInv := big.NewInt(0).ModInverse(a.Val, LElem256BigMod.Val)
	if bInv == nil {
		panic("ModInverse does not exist")
	}
	return LElem256Big{bInv}
}

func (a LElem256Big) ToBigInt() *big.Int {
	return new(big.Int).Set(a.Val)
}

func (a LElem256Big) FromBigInt(n *big.Int) RElem {
	//return LElem256Big{new(big.Int).Mod(new(big.Int).Set(n), LElem256BigMod.Val)}
	return LElem256Big{new(big.Int).Set(n)}
}

func (a LElem256Big) ToBigFloat(fracBits int) *big.Float {
	out := new(big.Float).SetInt(a.Val)
	out.Mul(out, new(big.Float).SetMantExp(big.NewFloat(1), -fracBits))
	return out
}

func (a LElem256Big) ToSignedBigInt() *big.Int {
	if a.Val.Cmp(LElem256BigModHalf.Val) >= 0 {
		return new(big.Int).Sub(a.Val, LElem256BigMod.Val)
	}
	return new(big.Int).Set(a.Val)
}

func (a LElem256Big) ToSignedBigFloat(fracBits int) *big.Float {
	out := new(big.Float).SetInt(a.ToSignedBigInt())
	out.Mul(out, new(big.Float).SetMantExp(big.NewFloat(1), -fracBits))
	return out
}

func (a LElem256Big) FromBigFloat(n *big.Float, fracBits int) RElem {
	f := new(big.Float).Mul(n, new(big.Float).SetMantExp(big.NewFloat(1), fracBits))
	return a.FromBigInt(roundFloat(f))
}
func (a LElem256Big) AssertTypeFor(n RElem) RElem {
	return n.(LElem256Big)
}

func (a LElem256Big) FromInt(n int) RElem {
	if n >= 0 {
		return LElem256Big{new(big.Int).SetInt64(int64(n))}
	}
	return LElem256Big{new(big.Int).Sub(LElem256BigMod.Val, new(big.Int).SetInt64(-int64(n)))}
}

func (a LElem256Big) FromUint64(n uint64) RElem {
	return LElem256Big{new(big.Int).SetUint64(n)}
}

func (a LElem256Big) FromFloat64(n float64, fracBits int) RElem {
	if n < 0 {
		return LElem256Big{LElem128Zero.FromFloat64(-n, fracBits).(LElem128).ToBigInt()}.Neg()
	} else {
		return LElem256Big{LElem128Zero.FromFloat64(n, fracBits).(LElem128).ToBigInt()}
	}
}

func (a LElem256Big) Float64(fracBits int) float64 {
	var sgn int
	var b LElem256Big
	if a.Val.Cmp(LElem256BigModHalf.Val) < 0 {
		sgn, b = 1, a
	} else {
		sgn, b = -1, a.Neg().(LElem256Big)
	}

	out, _ := new(big.Rat).SetFrac(b.Val, new(big.Int).Lsh(big.NewInt(1), uint(fracBits))).Float64()
	return float64(sgn) * out
}

func (a LElem256Big) FromBytes(buf []byte) RElem {
	return LElem256Big{new(big.Int).SetBytes(buf[:a.NumBytes()])}
}

func (a LElem256Big) ToBytes(buf []byte) {
	if a.Val.Sign() < 0 {
		a.Val.Mod(a.Val, LElem256BigMod.Val)
	}
	a.Val.FillBytes(buf[:a.NumBytes()])
}

func (a LElem256Big) Uint64() uint64 {
	return a.Val.Uint64()
}

func (a LElem256Big) Zero() RElem {
	return LElem256Big{big.NewInt(0)}
}

func (a LElem256Big) One() RElem {
	return LElem256Big{big.NewInt(1)}
}

func (a LElem256Big) Modulus() *big.Int {
	return new(big.Int).Set(LElem256BigMod.Val)
}

func (a LElem256Big) NumBytes() uint32 {
	return LElem256BigBytes
}

func (a LElem256Big) ModBitLength() int {
	return LElem256BigModBitLen
}

func (a LElem256Big) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB >= a.ModBitLength() {
		panic("Invalid bit position")
	}
	return a.Val.Bit(posFromLSB)
}

func (a LElem256Big) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return LElem256Big{new(big.Int).Set(a.Val)}
	}

	return LElem256Big{new(big.Int).Rem(a.Val, new(big.Int).Lsh(big.NewInt(1), uint(nBits)))}
}

func (a LElem256Big) TypeID() uint8 {
	return LElem256BigUniqueID
}

func (a LElem256Big) Rand(prg *frand.RNG) RElem {
	bitLen := a.ModBitLength()
	numBytes := 1 + ((bitLen - 1) >> 3)
	buf := make([]byte, numBytes)
	headLen := bitLen - ((numBytes - 1) << 3)

again:
	prg.Read(buf)
	if headLen < 8 { // SetBytes is big-endian
		buf[0] = uint8(buf[0]) % (1 << headLen)
	}
	r := new(big.Int).SetBytes(buf)
	if r.Cmp(LElem256BigMod.Val) >= 0 {
		goto again
	}
	return LElem256Big{new(big.Int).SetBytes(buf)}
}

func (a LElem256Big) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits >= a.ModBitLength() {
		panic("Requested bit length is larger than modulus")
	}

	buf := make([]byte, 1+((nbits-1)/8))
	prg.Read(buf)

	return LElem256Big{new(big.Int).Rem(new(big.Int).SetBytes(buf), new(big.Int).Lsh(big.NewInt(1), uint(nbits)))}
}

func (a LElem256Big) Copy() RElem {
	return LElem256Big{new(big.Int).Set(a.Val)}
}
