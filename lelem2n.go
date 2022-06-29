package mpc_core

import (
	"encoding/binary"
	"math"
	"math/big"
	"unsafe"

	"github.com/hhcho/frand"
)

// Large ring with power of two modulus
type LElem2N uint64

const LElem2NBytes uint32 = uint32(unsafe.Sizeof(LElem2N(0)))
const LElem2NUniqueID uint8 = 0

var LElem2NModHalf LElem2N = LElem2N(uint64(1) << 63)
var LElem2NModBitLen int = 65
var LElem2NModBig *big.Int = big.NewInt(0).Add(big.NewInt(0).SetUint64(^uint64(0)), big.NewInt(1))

func (a LElem2N) Mul(b interface{}) RElem {
	return a * b.(LElem2N)
}

func (a LElem2N) Add(b interface{}) RElem {
	return a + b.(LElem2N)
}

func (a LElem2N) Sub(b interface{}) RElem {
	return a - b.(LElem2N)
}

func (a LElem2N) Neg() RElem {
	return -a
}

func (a LElem2N) Inv() RElem {
	panic("Modular inverse attempted for LElem2N")
	return LElem2N(0)
}

func (a LElem2N) AssertTypeFor(n RElem) RElem {
	return n.(LElem2N)
}

func (a LElem2N) FromInt(n int) RElem {
	return LElem2N(n)
}

func (a LElem2N) FromUint64(n uint64) RElem {
	return LElem2N(n)
}

func (a LElem2N) FromFloat64(n float64, fracBits int) RElem {
	sgn := 1
	if n < 0 {
		sgn = -1
		n = -n
	}

	x := math.Round(n * float64(uint64(1)<<fracBits))
	out := LElem2N(x)
	if sgn < 0 {
		return out.Neg()
	}
	return out
}

func (a LElem2N) Float64(fracBits int) float64 {
	var sgn int
	var b LElem2N
	if a < LElem2NModHalf {
		sgn, b = 1, a
	} else {
		sgn, b = -1, a.Neg().(LElem2N)
	}
	return float64(sgn) * (float64(b) / float64(uint64(1)<<fracBits))
}

func (a LElem2N) FromBytes(buf []byte) RElem {
	return LElem2N(binary.LittleEndian.Uint64(buf))
}

func (a LElem2N) ToBytes(buf []byte) {
	binary.LittleEndian.PutUint64(buf, uint64(a))
}

func (a LElem2N) Uint64() uint64 {
	return uint64(a)
}

func (a LElem2N) Zero() RElem {
	return LElem2N(0)
}

func (a LElem2N) One() RElem {
	return LElem2N(1)
}

func (a LElem2N) Modulus() *big.Int {
	return new(big.Int).Set(LElem2NModBig)
}

func (a LElem2N) NumBytes() uint32 {
	return LElem2NBytes
}

func (a LElem2N) ModBitLength() int {
	return LElem2NModBitLen
}

func (a LElem2N) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB >= a.ModBitLength() {
		panic("Invalid bit position")
	}
	mask := LElem2N(1) << posFromLSB
	return boolToUint((a & mask) > 0)
}

func (a LElem2N) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return a
	}
	return a % (LElem2N(1) << nBits)
}

func (a LElem2N) TypeID() uint8 {
	return LElem2NUniqueID
}
func (a LElem2N) Rand(prg *frand.RNG) RElem {
	buf := make([]byte, a.NumBytes())
	return randBytes(a, buf, prg)
}

func (a LElem2N) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits >= a.ModBitLength() {
		panic("Requested bit length is larger than modulus")
	}
	if nbits == a.ModBitLength()-1 { // Full range
		return a.Rand(prg)
	}

	buf := make([]byte, a.NumBytes())
	return randBytes(a, buf, prg).(LElem2N) % (LElem2N(1) << nbits)
}

func (a LElem2N) FromBigInt(n *big.Int) RElem {
	panic("Not implemented")
	return nil
}

func (a LElem2N) Copy() RElem {
	return a
}
