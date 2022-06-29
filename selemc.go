package mpc_core

import (
	"math/big"
	"math/bits"
	"unsafe"

	"github.com/hhcho/frand"
)

// Small ring used for comparison
type SElemC uint8

const SElemCBytes uint32 = uint32(unsafe.Sizeof(SElemC(0)))
const SElemCRandBnd SElemC = ^SElemC(0) - ((^SElemC(0)) % SElemCMod)
const SElemCMod SElemC = 199
const SElemCMod16 uint16 = uint16(SElemCMod)
const SElemCUniqueID uint8 = 3

var SElemCModBitLen int = bits.Len64(uint64(SElemCMod))
var SElemCModBig *big.Int = big.NewInt(0).SetUint64(uint64(SElemCMod))

func (a SElemC) Mul(b interface{}) RElem {
	return SElemC((uint16(a) * uint16(b.(SElemC))) % SElemCMod16)
}

func (a SElemC) Add(b interface{}) RElem {
	return SElemC((uint16(a) + uint16(b.(SElemC))) % SElemCMod16)
}

func (a SElemC) Sub(b interface{}) RElem {
	return a.Add(b.(SElemC).Neg())
}

func (a SElemC) Neg() RElem {
	if a == 0 {
		return a
	}
	return SElemC(SElemCMod) - a
}

func (a SElemC) Inv() RElem {
	b := big.NewInt(0).SetUint64(uint64(a))
	bInv := b.ModInverse(b, SElemCModBig)
	if bInv == nil {
		panic("ModInverse does not exist")
	}
	return SElemC(bInv.Uint64())
}

func (a SElemC) AssertTypeFor(n RElem) RElem {
	return n.(SElemC)
}

func (a SElemC) FromInt(n int) RElem {
	if n >= 0 {
		return SElemC(0).Add(SElemC(n))
	}
	return SElemC(0).Sub(SElemC(-n))
}

func (a SElemC) FromUint64(n uint64) RElem {
	return SElemC(n % uint64(SElemCMod))
}

func (a SElemC) FromFloat64(n float64, fracBits int) RElem {
	panic("SElemC is not meant for storing fractional values")
	return a
}

func (a SElemC) Float64(fracBits int) float64 {
	panic("SElemC is not meant for storing fractional values")
	return 0
}
func (a SElemC) FromBytes(buf []byte) RElem {
	return SElemC(buf[0])
}

func (a SElemC) ToBytes(buf []byte) {
	buf[0] = byte(a)
}

func (a SElemC) Uint64() uint64 {
	return uint64(a)
}

func (a SElemC) Zero() RElem {
	return SElemC(0)
}

func (a SElemC) One() RElem {
	return SElemC(1)
}

func (a SElemC) Modulus() *big.Int {
	return new(big.Int).Set(SElemCModBig)
}

func (a SElemC) NumBytes() uint32 {
	return SElemCBytes
}

func (a SElemC) ModBitLength() int {
	return SElemCModBitLen
}

func (a SElemC) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB > a.ModBitLength() {
		panic("Invalid bit position")
	}
	return boolToUint((a & (SElemC(1) << posFromLSB)) > 0)
}

func (a SElemC) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return a
	}
	return a % (SElemC(1) << nBits)
}

func (a SElemC) TypeID() uint8 {
	return SElemCUniqueID
}

func (a SElemC) Rand(prg *frand.RNG) RElem {
	buf := make([]byte, a.NumBytes())
again:
	r := randBytes(a, buf, prg).(SElemC)
	if r >= SElemCRandBnd {
		goto again
	}
	return r % SElemCMod
}

func (a SElemC) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits >= a.ModBitLength() {
		panic("Requested bit length is larger than modulus")
	}
	buf := make([]byte, a.NumBytes())
	return randBytes(a, buf, prg).(SElemC) % (SElemC(1) << nbits)
}

func (a SElemC) FromBigInt(n *big.Int) RElem {
	panic("Not implemented")
	return nil
}
