package mpc_core

import (
	"encoding/binary"
	"math"
	"math/big"
	"math/bits"
	"unsafe"

	"github.com/hhcho/frand"
)

// Large ring with prime modulus
type LElemP uint64

const LElemPBytes uint32 = uint32(unsafe.Sizeof(LElemP(0)))
const LElemPRandBnd LElemP = ^LElemP(0) - ((^LElemP(0)) % LElemPMod)
const LElemPMod LElemP = 2305843009213693951     // 2^61 - 1
const LElemPModHalf LElemP = 1152921504606846976 // 2^60
const LElemPMod64 uint64 = uint64(LElemPMod)
const LElemPUniqueID uint8 = 1

var LElemPModBitLen int = bits.Len64(uint64(LElemPMod))
var LElemPModBig *big.Int = big.NewInt(0).SetUint64(uint64(LElemPMod))

func (a LElemP) Mul(b interface{}) RElem {
	hi, lo := bits.Mul64(uint64(a), uint64(b.(LElemP)))

	if hi < LElemPMod64 {
		_, r := bits.Div64(hi, lo, LElemPMod64)
		return LElemP(r)
	}

	_, r := bits.Div64(0, hi, LElemPMod64)
	_, r = bits.Div64(r, lo, LElemPMod64)
	return LElemP(r)
}

func (a LElemP) Add(b interface{}) RElem {
	c, carry := bits.Add64(uint64(a), uint64(b.(LElemP)), 0)
	_, r := bits.Div64(carry, c, LElemPMod64)
	return LElemP(r)
}

func (a LElemP) Sub(b interface{}) RElem {
	return a.Add(b.(LElemP).Neg())
}

func (a LElemP) Neg() RElem {
	if a == 0 {
		return a
	}
	return LElemP(LElemPMod) - a
}

func (a LElemP) Inv() RElem {
	b := big.NewInt(0).SetUint64(uint64(a))
	bInv := b.ModInverse(b, LElemPModBig)
	if bInv == nil {
		panic("ModInverse does not exist")
	}
	return LElemP(bInv.Uint64())
}

func (a LElemP) AssertTypeFor(n RElem) RElem {
	return n.(LElemP)
}

func (a LElemP) FromInt(n int) RElem {
	if n >= 0 {
		return LElemP(0).Add(LElemP(n))
	}
	return LElemP(0).Sub(LElemP(-n))
}

func (a LElemP) FromUint64(n uint64) RElem {
	return LElemP(n % LElemPMod64)
}

func (a LElemP) FromFloat64(n float64, fracBits int) RElem {
	sgn := 1
	if n < 0 {
		sgn = -1
		n = -n
	}

	x := math.Round(n * float64(uint64(1)<<fracBits))
	out := LElemP(x)
	if sgn < 0 {
		return out.Neg()
	}
	return out
}

func (a LElemP) Float64(fracBits int) float64 {
	var sgn int
	var b LElemP
	if a < LElemPModHalf {
		sgn, b = 1, a
	} else {
		sgn, b = -1, a.Neg().(LElemP)
	}
	return float64(sgn) * (float64(b) / float64(uint64(1)<<fracBits))
}

func (a LElemP) FromBytes(buf []byte) RElem {
	return LElemP(binary.LittleEndian.Uint64(buf))
}

func (a LElemP) ToBytes(buf []byte) {
	binary.LittleEndian.PutUint64(buf, uint64(a))
}

func (a LElemP) Uint64() uint64 {
	return uint64(a)
}

func (a LElemP) Zero() RElem {
	return LElemP(0)
}

func (a LElemP) One() RElem {
	return LElemP(1)
}

func (a LElemP) Modulus() *big.Int {
	return new(big.Int).Set(LElemPModBig)
}

func (a LElemP) NumBytes() uint32 {
	return LElemPBytes
}

func (a LElemP) ModBitLength() int {
	return LElemPModBitLen
}

func (a LElemP) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB > a.ModBitLength() {
		panic("Invalid bit position")
	}
	return boolToUint((a & (LElemP(1) << posFromLSB)) > 0)
}

func (a LElemP) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return a
	}
	return a % (LElemP(1) << nBits)
}

func (a LElemP) TypeID() uint8 {
	return LElemPUniqueID
}

func (a LElemP) Rand(prg *frand.RNG) RElem {
	buf := make([]byte, a.NumBytes())
again:
	r := randBytes(a, buf, prg).(LElemP)
	if r >= LElemPRandBnd {
		goto again
	}
	return r % LElemPMod
}

func (a LElemP) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits >= a.ModBitLength() {
		panic("Requested bit length is larger than modulus")
	}
	buf := make([]byte, a.NumBytes())
	return randBytes(a, buf, prg).(LElemP) % (LElemP(1) << nbits)
}

func (a LElemP) FromBigInt(n *big.Int) RElem {
	panic("Not implemented")
	return nil
}
