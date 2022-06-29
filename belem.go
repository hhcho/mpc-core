package mpc_core

import (
	"encoding/binary"
	"math/big"
	"unsafe"

	"github.com/hhcho/frand"
)

// Bitwise sharing
type BElem uint64

const BElemUniqueID uint8 = 6

var BElemBytes uint32 = uint32(unsafe.Sizeof(BElem(0)))
var BElemBits int = int(8 * BElemBytes)

func (a BElem) Mul(b interface{}) RElem {
	return BElem(uint64(a) & uint64(b.(BElem)))
}

func (a BElem) Add(b interface{}) RElem {
	return BElem(uint64(a) ^ uint64(b.(BElem)))
}

func (a BElem) Sub(b interface{}) RElem {
	return BElem(uint64(a) ^ uint64(b.(BElem)))
}

func (a BElem) Neg() RElem {
	return a
}

func (a BElem) Inv() RElem {
	panic("Modular inverse attempted for BElem")
	return BElem(0)
}

func (a BElem) AssertTypeFor(n RElem) RElem {
	return n.(BElem)
}

func (a BElem) FromInt(n int) RElem {
	return BElem(n)
}

func (a BElem) FromUint64(n uint64) RElem {
	return BElem(n)
}

func (a BElem) FromFloat64(n float64, fracBits int) RElem {
	panic("FromFloat64 for BElem undefined")
	return a
}

func (a BElem) Float64(fracBits int) float64 {
	panic("Float64 for BElem undefined")
	return 0
}

func (a BElem) FromBytes(buf []byte) RElem {
	return BElem(binary.LittleEndian.Uint64(buf))
}

func (a BElem) ToBytes(buf []byte) {
	binary.LittleEndian.PutUint64(buf, uint64(a))
}

func (a BElem) Uint64() uint64 {
	return uint64(a)
}

func (a BElem) Zero() RElem {
	return BElem(0)
}

func (a BElem) One() RElem {
	return BElem(^uint64(0))
}

func (a BElem) Modulus() *big.Int {
	panic("Modulus for BElem undefined")
	return new(big.Int)
}

func (a BElem) NumBytes() uint32 {
	return BElemBytes
}

func (a BElem) ModBitLength() int {
	panic("ModBitLength for BElem undefined")
	return 0
}

func (a BElem) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB >= 8*int(a.NumBytes()) {
		panic("Invalid bit position")
	}
	mask := uint64(1) << posFromLSB
	return boolToUint((uint64(a) & mask) > 0)
}

func (a BElem) Trunc(nBits int) RElem {
	if nBits < 0 || nBits >= 8*int(a.NumBytes()) {
		panic("Invalid number of bits")
	} else if nBits == 8*int(a.NumBytes()) {
		return a
	}
	return BElem(uint64(a) % (uint64(1) << nBits))
}

func (a BElem) TypeID() uint8 {
	return BElemUniqueID
}
func (a BElem) Rand(prg *frand.RNG) RElem {
	buf := make([]byte, a.NumBytes())
	return randBytes(a, buf, prg)
}

func (a BElem) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits > 8*int(a.NumBytes()) {
		panic("Requested bit length is larger than BElem size")
	}
	return a.Rand(prg).Trunc(nbits)
}

func (a BElem) FromBigInt(n *big.Int) RElem {
	panic("Not implemented")
	return nil
}

func (a BElem) Copy() RElem {
	return a
}
