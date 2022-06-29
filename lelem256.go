package mpc_core

import (
	"math/big"

	"github.com/hhcho/frand"
)

type LElem256 struct {
	Val *big.Int
}

const LElem256Bytes uint32 = 32 // Extra byte for sign bit for ease of marshalling
var LElem256Zero LElem256 = LElem256{big.NewInt(0)}
var LElem256Mod = LElem256{new(big.Int).Sub(big.NewInt(0).Lsh(big.NewInt(1), 256), big.NewInt(189))}    // 2^256 - 189
var LElem256ModHalf = LElem256{new(big.Int).Sub(big.NewInt(0).Lsh(big.NewInt(1), 255), big.NewInt(94))} // 2^255 - 94
var LElem256ModBitLen int = (LElem256Mod).Val.BitLen()

const LElem256UniqueID uint8 = 5

/* LElem256 */

func (a LElem256) Mul(b interface{}) RElem {
	m := new(big.Int).Mul(a.Val, b.(LElem256).Val)
	return LElem256{new(big.Int).Mod(m, LElem256Mod.Val)}
}

func (a LElem256) Add(b interface{}) RElem {
	m := new(big.Int).Add(a.Val, b.(LElem256).Val)
	return LElem256{new(big.Int).Mod(m, LElem256Mod.Val)}
}

func (a LElem256) Sub(b interface{}) RElem {
	return a.Add(b.(LElem256).Neg())
}

func (a LElem256) Neg() RElem {
	m := new(big.Int).Neg(a.Val)
	return LElem256{new(big.Int).Mod(m, LElem256Mod.Val)}
}

func (a LElem256) Inv() RElem {
	bInv := big.NewInt(0).ModInverse(a.Val, LElem256Mod.Val)
	if bInv == nil {
		panic("ModInverse does not exist")
	}
	return LElem256{bInv}
}

func (a LElem256) ToBigInt() *big.Int {
	return new(big.Int).Set(a.Val)
}

func (a LElem256) FromBigInt(n *big.Int) RElem {
	//return LElem256{new(big.Int).Mod(new(big.Int).Set(n), LElem256Mod.Val)}
	return LElem256{new(big.Int).Set(n)}
}

func (a LElem256) ToBigFloat(fracBits int) *big.Float {
	out := new(big.Float).SetInt(a.Val)
	out.Mul(out, new(big.Float).SetMantExp(big.NewFloat(1), -fracBits))
	return out
}

func (a LElem256) ToSignedBigInt() *big.Int {
	if a.Val.Cmp(LElem256ModHalf.Val) >= 0 {
		return new(big.Int).Sub(a.Val, LElem256Mod.Val)
	}
	return new(big.Int).Set(a.Val)
}

func (a LElem256) ToSignedBigFloat(fracBits int) *big.Float {
	out := new(big.Float).SetInt(a.ToSignedBigInt())
	out.Mul(out, new(big.Float).SetMantExp(big.NewFloat(1), -fracBits))
	return out
}

func (a LElem256) FromBigFloat(n *big.Float, fracBits int) RElem {
	f := new(big.Float).Mul(n, new(big.Float).SetMantExp(big.NewFloat(1), fracBits))
	return a.FromBigInt(roundFloat(f))
}
func (a LElem256) AssertTypeFor(n RElem) RElem {
	return n.(LElem256)
}

func (a LElem256) FromInt(n int) RElem {
	if n >= 0 {
		return LElem256{new(big.Int).SetInt64(int64(n))}
	}
	return LElem256{new(big.Int).Sub(LElem256Mod.Val, new(big.Int).SetInt64(-int64(n)))}
}

func (a LElem256) FromUint64(n uint64) RElem {
	return LElem256{new(big.Int).SetUint64(n)}
}

func (a LElem256) FromFloat64(n float64, fracBits int) RElem {
	if n < 0 {
		return LElem256{LElem128Zero.FromFloat64(-n, fracBits).(LElem128).ToBigInt()}.Neg()
	} else {
		return LElem256{LElem128Zero.FromFloat64(n, fracBits).(LElem128).ToBigInt()}
	}
}

func (a LElem256) Float64(fracBits int) float64 {
	var sgn int
	var b LElem256
	if a.Val.Cmp(LElem256ModHalf.Val) < 0 {
		sgn, b = 1, a
	} else {
		sgn, b = -1, a.Neg().(LElem256)
	}

	out, _ := new(big.Rat).SetFrac(b.Val, new(big.Int).Lsh(big.NewInt(1), uint(fracBits))).Float64()
	return float64(sgn) * out
}

func (a LElem256) FromBytes(buf []byte) RElem {
	return LElem256{new(big.Int).SetBytes(buf[:a.NumBytes()])}
}

func (a LElem256) ToBytes(buf []byte) {
	if a.Val.Sign() < 0 {
		a.Val.Mod(a.Val, LElem256Mod.Val)
	}
	a.Val.FillBytes(buf[:a.NumBytes()])
}

func (a LElem256) Uint64() uint64 {
	return a.Val.Uint64()
}

func (a LElem256) Zero() RElem {
	return LElem256{big.NewInt(0)}
}

func (a LElem256) One() RElem {
	return LElem256{big.NewInt(1)}
}

func (a LElem256) Modulus() *big.Int {
	return new(big.Int).Set(LElem256Mod.Val)
}

func (a LElem256) NumBytes() uint32 {
	return LElem256Bytes
}

func (a LElem256) ModBitLength() int {
	return LElem256ModBitLen
}

func (a LElem256) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB >= a.ModBitLength() {
		panic("Invalid bit position")
	}
	return a.Val.Bit(posFromLSB)
}

func (a LElem256) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return LElem256{new(big.Int).Set(a.Val)}
	}

	return LElem256{new(big.Int).Rem(a.Val, new(big.Int).Lsh(big.NewInt(1), uint(nBits)))}
}

func (a LElem256) TypeID() uint8 {
	return LElem256UniqueID
}

func (a LElem256) Rand(prg *frand.RNG) RElem {
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
	if r.Cmp(LElem256Mod.Val) >= 0 {
		goto again
	}
	return LElem256{new(big.Int).SetBytes(buf)}
}

func (a LElem256) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits >= a.ModBitLength() {
		panic("Requested bit length is larger than modulus")
	}

	buf := make([]byte, 1+((nbits-1)/8))
	prg.Read(buf)

	return LElem256{new(big.Int).Rem(new(big.Int).SetBytes(buf), new(big.Int).Lsh(big.NewInt(1), uint(nbits)))}
}

func (a LElem256) Copy() RElem {
	return LElem256{new(big.Int).Set(a.Val)}
}
