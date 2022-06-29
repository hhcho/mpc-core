package mpc_core

import (
	"math/big"
	"math/bits"
	"unsafe"

	"github.com/hhcho/frand"
)

// Small ring used for division and square root
type SElemDS uint8

const SElemDSBytes uint32 = uint32(unsafe.Sizeof(SElemDS(0)))
const SElemDSRandBnd SElemDS = ^SElemDS(0) - ((^SElemDS(0)) % SElemDSMod)
const SElemDSMod SElemDS = 199
const SElemDSMod16 uint16 = uint16(SElemDSMod)
const SElemDSUniqueID uint8 = 2

var SElemDSModBitLen int = bits.Len64(uint64(SElemDSMod))
var SElemDSModBig *big.Int = big.NewInt(0).SetUint64(uint64(SElemDSMod))

// Mod64 returns r = u%v.
//func (u Uint128) QuoRem64(v uint64) (q Uint128, r uint64) {
//	if u.Hi < v {
//		q.Lo, r = bits.Div64(u.Hi, u.Lo, v)
//	} else {
//		q.Hi, r = bits.Div64(0, u.Hi, v)
//		q.Lo, r = bits.Div64(r, u.Lo, v)
//	}
//	return
//}

func (a SElemDS) Mul(b interface{}) RElem {
	return SElemDS((uint16(a) * uint16(b.(SElemDS))) % SElemDSMod16)
}

func (a SElemDS) Add(b interface{}) RElem {
	return SElemDS((uint16(a) + uint16(b.(SElemDS))) % SElemDSMod16)
}

func (a SElemDS) Sub(b interface{}) RElem {
	return a.Add(b.(SElemDS).Neg())
}

func (a SElemDS) Neg() RElem {
	if a == 0 {
		return a
	}
	return SElemDS(SElemDSMod) - a
}

func (a SElemDS) Inv() RElem {
	b := big.NewInt(0).SetUint64(uint64(a))
	bInv := b.ModInverse(b, SElemDSModBig)
	if bInv == nil {
		panic("ModInverse does not exist")
	}
	return SElemDS(bInv.Uint64())
}

func (a SElemDS) AssertTypeFor(n RElem) RElem {
	return n.(SElemDS)
}
func (a SElemDS) FromInt(n int) RElem {
	if n >= 0 {
		return SElemDS(0).Add(SElemDS(n))
	}
	return SElemDS(0).Sub(SElemDS(-n))
}

func (a SElemDS) FromUint64(n uint64) RElem {
	return SElemDS(n % uint64(SElemDSMod))
}

func (a SElemDS) FromFloat64(n float64, fracBits int) RElem {
	panic("SElemDS is not meant for storing fractional values")
	return a
}
func (a SElemDS) Float64(fracBits int) float64 {
	panic("SElemDS is not meant for storing fractional values")
	return 0
}

func (a SElemDS) FromBytes(buf []byte) RElem {
	return SElemDS(buf[0])
}
func (a SElemDS) ToBytes(buf []byte) {
	buf[0] = byte(a)
}

func (a SElemDS) Uint64() uint64 {
	return uint64(a)
}

func (a SElemDS) Zero() RElem {
	return SElemDS(0)
}

func (a SElemDS) One() RElem {
	return SElemDS(1)
}

func (a SElemDS) Modulus() *big.Int {
	return new(big.Int).Set(SElemDSModBig)
}

func (a SElemDS) NumBytes() uint32 {
	return SElemDSBytes
}

func (a SElemDS) ModBitLength() int {
	return SElemDSModBitLen
}

func (a SElemDS) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB > a.ModBitLength() {
		panic("Invalid bit position")
	}
	return boolToUint((a & (SElemDS(1) << posFromLSB)) > 0)
}

func (a SElemDS) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return a
	}
	return a % (SElemDS(1) << nBits)
}

func (a SElemDS) TypeID() uint8 {
	return SElemDSUniqueID
}

func (a SElemDS) Rand(prg *frand.RNG) RElem {
	buf := make([]byte, a.NumBytes())
again:
	r := randBytes(a, buf, prg).(SElemDS)
	if r >= SElemDSRandBnd {
		goto again
	}
	return r % SElemDSMod
}

func (a SElemDS) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits >= a.ModBitLength() {
		panic("Requested bit length is larger than modulus")
	}
	buf := make([]byte, a.NumBytes())
	return randBytes(a, buf, prg).(SElemDS) % (SElemDS(1) << nbits)
}
func (a SElemDS) FromBigInt(n *big.Int) RElem {
	panic("Not implemented")
	return nil
}

func (a SElemDS) Copy() RElem {
	return a
}
