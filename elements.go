package mpc_core

import (
	"encoding/binary"
	"math"
	"math/big"
	"math/bits"
	"unsafe"

	"github.com/hhcho/frand"
)

type RElem interface {
	Mul(interface{}) RElem
	Add(interface{}) RElem
	Sub(interface{}) RElem
	Neg() RElem
	Inv() RElem
	Zero() RElem
	One() RElem
	Modulus() *big.Int
	Uint64() uint64
	Float64(int) float64
	FromInt(int) RElem
	FromBigInt(*big.Int) RElem
	FromUint64(uint64) RElem
	FromFloat64(float64, int) RElem
	FromBytes([]byte) RElem
	ToBytes([]byte)
	AssertTypeFor(RElem) RElem
	NumBytes() uint32
	Rand(*frand.RNG) RElem
	RandBits(*frand.RNG, int) RElem
	TypeID() uint8
	ModBitLength() int
	GetBit(int) uint
	Trunc(int) RElem
	Copy() RElem
}

type RVec []RElem
type RMat []RVec

// Bitwise sharing
type BElem uint64

const BElemUniqueID uint8 = 6

var BElemBytes uint32 = uint32(unsafe.Sizeof(BElem(0)))
var BElemBits int = int(8 * BElemBytes)

type LElem2NBigInt struct {
	Val *big.Int
}

const LElem2NBigIntUniqueID uint8 = 7

var LElem2NBigIntZero LElem2NBigInt = LElem2NBigInt{big.NewInt(0)}
var LElem2NBigIntMod LElem2NBigInt     // Variable set on the fly
var LElem2NBigIntModHalf LElem2NBigInt // Variable set on the fly
var LElem2NBigIntModBitLen int         // Variable set on the fly

type LElem256 struct {
	Val *big.Int
}

const LElem256Bytes uint32 = 32 // Extra byte for sign bit for ease of marshalling
var LElem256Zero LElem256 = LElem256{big.NewInt(0)}
var LElem256Mod = LElem256{new(big.Int).Sub(big.NewInt(0).Lsh(big.NewInt(1), 256), big.NewInt(189))}    // 2^256 - 189
var LElem256ModHalf = LElem256{new(big.Int).Sub(big.NewInt(0).Lsh(big.NewInt(1), 255), big.NewInt(94))} // 2^255 - 94
var LElem256ModBitLen int = (LElem256Mod).Val.BitLen()

const LElem256UniqueID uint8 = 5

type Uint128 struct {
	Hi, Lo uint64
}

type LElem128 Uint128

const LElem128Bytes uint32 = 16

func initLElem128RandBnd() LElem128 {
	mOne := Uint128{^uint64(0), ^uint64(0)}
	_, rem := div128(Uint128{0, 0}, mOne, Uint128(LElem128Mod))
	out, _ := sub128(mOne, rem)
	return LElem128(out)
}

var LElem128RandBnd LElem128 = initLElem128RandBnd()
var LElem128Zero LElem128 = LElem128{0, 0}
var LElem128Mod LElem128 = LElem128{9223372036854775807, 18446744073709551615} // 2^127 - 1
var LElem128ModHalf LElem128 = LElem128{4611686018427387904, 0}                // 2^126
const LElem128UniqueID uint8 = 4

var LElem128ModBitLen int = bits.Len64(LElem128Mod.Hi) + 64
var LElem128ModBig *big.Int = big.NewInt(0).Sub(big.NewInt(0).Exp(big.NewInt(2), big.NewInt(127), nil), big.NewInt(1))
var LElem128ModHalfBig *big.Int = big.NewInt(0).Exp(big.NewInt(2), big.NewInt(126), nil)

// Large ring with power of two modulus
type LElem2N uint64

const LElem2NBytes uint32 = uint32(unsafe.Sizeof(LElem2N(0)))
const LElem2NUniqueID uint8 = 0

var LElem2NModHalf LElem2N = LElem2N(uint64(1) << 63)
var LElem2NModBitLen int = 65
var LElem2NModBig *big.Int = big.NewInt(0).Add(big.NewInt(0).SetUint64(^uint64(0)), big.NewInt(1))

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

// Small ring used for division and square root
type SElemDS uint8

const SElemDSBytes uint32 = uint32(unsafe.Sizeof(SElemDS(0)))
const SElemDSRandBnd SElemDS = ^SElemDS(0) - ((^SElemDS(0)) % SElemDSMod)
const SElemDSMod SElemDS = 199
const SElemDSMod16 uint16 = uint16(SElemDSMod)
const SElemDSUniqueID uint8 = 2

var SElemDSModBitLen int = bits.Len64(uint64(SElemDSMod))
var SElemDSModBig *big.Int = big.NewInt(0).SetUint64(uint64(SElemDSMod))

// Small ring used for comparison
type SElemC uint8

const SElemCBytes uint32 = uint32(unsafe.Sizeof(SElemC(0)))
const SElemCRandBnd SElemC = ^SElemC(0) - ((^SElemC(0)) % SElemCMod)
const SElemCMod SElemC = 199
const SElemCMod16 uint16 = uint16(SElemCMod)
const SElemCUniqueID uint8 = 3

var SElemCModBitLen int = bits.Len64(uint64(SElemCMod))
var SElemCModBig *big.Int = big.NewInt(0).SetUint64(uint64(SElemCMod))

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

/* LElem128 */

func (a LElem128) ToBigInt() *big.Int {
	r := new(big.Int)
	new(big.Int).QuoRem(new(big.Int).Add(new(big.Int).Lsh(new(big.Int).SetUint64(a.Hi), 64), new(big.Int).SetUint64(a.Lo)), LElem128ModBig, r)
	return r
}

func (a LElem128) FromBigInt(n *big.Int) RElem {
	lo := big.NewInt(0).Mod(n, big.NewInt(0).Lsh(big.NewInt(1), 64))
	hi := big.NewInt(0).Rsh(n, 64)
	return LElem128(Uint128{hi.Uint64(), lo.Uint64()})
}

func (a LElem128) ToBigFloat(fracBits int) *big.Float {
	out := new(big.Float).SetInt(a.ToBigInt())
	out.Mul(out, new(big.Float).SetMantExp(big.NewFloat(1), -fracBits))
	return out
}

func (a LElem128) ToSignedBigInt() *big.Int {
	tmp := a.ToBigInt()
	half := new(big.Int).Rsh(LElem128ModBig, 1)
	if tmp.Cmp(half) >= 0 {
		tmp.Sub(tmp, LElem128ModBig)
	}
	return tmp
}
func (a LElem128) ToSignedBigFloat(fracBits int) *big.Float {
	out := new(big.Float).SetInt(a.ToSignedBigInt())
	out.Mul(out, new(big.Float).SetMantExp(big.NewFloat(1), -fracBits))
	return out
}

var halfFloat = big.NewFloat(0.5)

func roundFloat(a *big.Float) *big.Int {
	var i *big.Int
	if a.Signbit() {
		i, _ = new(big.Float).Sub(a, halfFloat).Int(nil)
	} else {
		i, _ = new(big.Float).Add(a, halfFloat).Int(nil)
	}
	return i
}
func (a LElem128) FromBigFloat(n *big.Float, fracBits int) RElem {
	f := new(big.Float).Mul(n, new(big.Float).SetMantExp(big.NewFloat(1), fracBits))
	return a.FromBigInt(roundFloat(f))
}

// From: https://github.com/lukechampine/uint128/blob/master/uint128.go
// QuoRem returns q = u/v and r = u%v.
//func (u Uint128) quoRem128(v Uint128) (q, r Uint128) {
//	if v.Hi == 0 {
//		var r64 uint64
//		q, r64 = u.QuoRem64(v.Lo)
//		r = From64(r64)
//	} else {
//		// generate a "trial quotient," guaranteed to be within 1 of the actual
//		// quotient, then adjust.
//		n := uint(bits.LeadingZeros64(v.Hi))
//		v1 := v.Lsh(n)
//		u1 := u.Rsh(1)
//		tq, _ := bits.Div64(u1.Hi, u1.Lo, v1.Hi)
//		tq >>= 63 - n
//		if tq != 0 {
//			tq--
//		}
//		q = From64(tq)
//		// calculate remainder using trial quotient, then adjust if remainder is
//		// greater than divisor
//		r = u.Sub(v.Mul64(tq))
//		if r.Cmp(v) >= 0 {
//			q = q.Add64(1)
//			r = r.Sub(v)
//		}
//	}
//	return
//}

func mul128(a, b Uint128) (Hi, Lo Uint128) {
	// (h11), (l11, h10, h01), (l10, l01, h00), (r0)
	h00, r0 := bits.Mul64(a.Lo, b.Lo)
	h01, l01 := bits.Mul64(a.Lo, b.Hi)
	h10, l10 := bits.Mul64(a.Hi, b.Lo)
	h11, l11 := bits.Mul64(a.Hi, b.Hi)

	t, carry1 := bits.Add64(h00, l01, 0)
	r1, carry2 := bits.Add64(t, l10, 0)
	t, carry1 = bits.Add64(l11, h01, carry1)
	r2, carry2 := bits.Add64(t, h10, carry2)
	r3 := h11 + carry1 + carry2

	Hi = Uint128{r3, r2}
	Lo = Uint128{r1, r0}
	return
}

func leadingZeros128(x Uint128) int {
	return bits.LeadingZeros64(x.Hi) + bits.LeadingZeros64(x.Lo)
}

// Lsh returns u<<n.
func lsh128(u Uint128, n uint) (s Uint128) {
	if n > 64 {
		s.Lo = 0
		s.Hi = u.Lo << (n - 64)
	} else {
		s.Lo = u.Lo << n
		s.Hi = u.Hi<<n | u.Lo>>(64-n)
	}
	return
}

// Rsh returns u>>n.
func rsh128(u Uint128, n uint) (s Uint128) {
	if n > 64 {
		s.Lo = u.Hi >> (n - 64)
		s.Hi = 0
	} else {
		s.Lo = u.Lo>>n | u.Hi<<(64-n)
		s.Hi = u.Hi >> n
	}
	return
}

func or128(a, b Uint128) (c Uint128) {
	c.Hi = a.Hi | b.Hi
	c.Lo = a.Lo | b.Lo
	return
}

func lessThan128(a, b Uint128) bool {
	if a.Hi == b.Hi {
		return a.Lo < b.Lo
	}
	return a.Hi < b.Hi
}

func add128(a, b Uint128) (c Uint128, carry uint64) {
	c.Lo, carry = bits.Add64(a.Lo, b.Lo, 0)
	c.Hi, carry = bits.Add64(a.Hi, b.Hi, carry)
	return
}

func sub128(a, b Uint128) (c Uint128, borrow uint64) {
	c.Lo, borrow = bits.Sub64(a.Lo, b.Lo, 0)
	c.Hi, borrow = bits.Sub64(a.Hi, b.Hi, borrow)
	return
}

// Requires Hi < y
func div128(Hi, Lo, y Uint128) (quo, rem Uint128) {
	//const (
	//	two32  = 1 << 32
	//	mask32 = two32 - 1
	//)
	if y.Hi == 0 && y.Lo == 0 {
		panic("Divide by zero error")
	}
	if !lessThan128(Hi, y) {
		panic("Overflow error")
	}

	s := uint(leadingZeros128(y))
	y = lsh128(y, s)

	yn1 := Uint128{0, y.Hi}                         //yn1 := y >> 32
	yn0 := Uint128{0, y.Lo}                         //yn0 := y & mask32
	un32 := or128(lsh128(Hi, s), rsh128(Lo, 128-s)) //un32 := Hi<<s | Lo>>(64-s)
	un10 := lsh128(Lo, s)                           //un10 := Lo << s
	un1 := Uint128{0, un10.Hi}                      //un1 := un10 >> 32
	un0 := Uint128{0, un10.Lo}                      //un0 := un10 & mask32
	q1hi, r := bits.Div64(0, un32.Hi, yn1.Lo)
	q1lo, _ := bits.Div64(r, un32.Lo, yn1.Lo)
	q1 := Uint128{q1hi, q1lo} //q1 := un32 / yn1
	_, t := mul128(q1, yn1)
	rhat, _ := sub128(un32, t) //rhat := un32 - q1*yn1

	_, t = mul128(q1, yn0)
	for q1.Hi != 0 || lessThan128(Uint128{rhat.Lo, un1.Lo}, t) {
		//for q1 >= two32 || q1*yn0 > two32*rhat+un1 {
		q1, _ = sub128(q1, Uint128{0, 1}) //	q1--
		rhat, _ = add128(rhat, yn1)       //rhat += yn1
		if rhat.Hi != 0 {                 //if rhat >= two32 {
			break
		}
	}

	_, q1y := mul128(q1, y)
	un21, _ := sub128(Uint128{un32.Lo, un1.Lo}, q1y) //un21 := un32*two32 + un1 - q1*y
	q0hi, r := bits.Div64(0, un21.Hi, yn1.Lo)
	q0lo, _ := bits.Div64(r, un21.Lo, yn1.Lo) //q0 := un21 / yn1
	q0 := Uint128{q0hi, q0lo}
	_, q0yn1 := mul128(q0, yn1)
	rhat, _ = sub128(un21, q0yn1) //rhat = un21 - q0*yn1

	_, t = mul128(q0, yn0)
	for q0.Hi != 0 || lessThan128(Uint128{rhat.Lo, un0.Lo}, t) {
		//for q0 >= two32 || q0*yn0 > two32*rhat+un0 {
		q0, _ = sub128(q0, Uint128{0, 1}) //	q0--
		rhat, _ = add128(rhat, yn1)       //rhat += yn1
		if rhat.Hi != 0 {                 //if rhat >= two32 {
			break
		}
	}

	quo = Uint128{q1.Lo, q0.Lo}

	_, q0y := mul128(q0, y)
	rem, _ = sub128(Uint128{un21.Lo, un0.Lo}, q0y)
	rem = rsh128(rem, s)

	return quo, rem
	//return q1*two32 + q0, (un21*two32 + un0 - q0*y) >> s
}

var LElem128MFormNInv = LElem128Zero.FromBigInt(new(big.Int).Sub(new(big.Int).Lsh(big.NewInt(1), 128), new(big.Int).ModInverse(LElem128ModBig, new(big.Int).Lsh(big.NewInt(1), 128))))
var LElem128MFormR2 = LElem128Zero.FromBigInt(new(big.Int).Mod(new(big.Int).Lsh(big.NewInt(1), 2*128), LElem128ModBig))

// REDC algorithm for modular arithmetic in Montgomery form
func (a LElem128) REDC(hi, lo Uint128) RElem {
	_, m := mul128(lo, Uint128(LElem128MFormNInv.(LElem128)))
	mNHi, mNLo := mul128(m, Uint128(LElem128Mod))
	_, carry := add128(mNLo, lo)
	t, carry1 := add128(mNHi, hi)
	t, carry2 := add128(t, Uint128{0, carry})

	if carry1+carry2 > 0 || !lessThan128(t, Uint128(LElem128Mod)) {
		res, _ := sub128(t, Uint128(LElem128Mod))
		return LElem128(res)
	} else {
		return LElem128(t)
	}
}

func (a LElem128) MulMForm(b interface{}) RElem {
	THi, TLo := mul128(Uint128(a), Uint128(b.(LElem128)))
	return a.REDC(THi, TLo)
}

func (a LElem128) ExitMForm() RElem {
	return a.REDC(Uint128{0, 0}, Uint128(a))
}

func (a LElem128) MForm() RElem {
	// Multiply by R^2 mod N
	Hi, Lo := mul128(Uint128(a), Uint128(LElem128MFormR2.(LElem128)))
	return a.REDC(Hi, Lo)
}

func (a LElem128) Mul(b interface{}) RElem {
	Hi, Lo := mul128(Uint128(a), Uint128(b.(LElem128)))

	if lessThan128(Hi, Uint128(LElem128Mod)) {
		_, r := div128(Hi, Lo, Uint128(LElem128Mod))
		return LElem128(r)
	}

	_, r := div128(Uint128{0, 0}, Hi, Uint128(LElem128Mod))
	_, r = div128(r, Lo, Uint128(LElem128Mod))
	return LElem128(r)
}

func (a LElem128) Add(b interface{}) RElem {
	out, carry := add128(Uint128(a), Uint128(b.(LElem128)))
	_, rem := div128(Uint128{0, carry}, out, Uint128(LElem128Mod))
	return LElem128(rem)
}

func (a LElem128) Sub(b interface{}) RElem {
	return a.Add(b.(LElem128).Neg())
}

func (a LElem128) Neg() RElem {
	if a.Hi == 0 && a.Lo == 0 {
		return a
	}
	out, _ := sub128(Uint128(LElem128Mod), Uint128(a))
	return LElem128(out)
}

func (a LElem128) Inv() RElem {
	bInv := big.NewInt(0).ModInverse(a.ToBigInt(), LElem128ModBig)
	if bInv == nil {
		panic("ModInverse does not exist")
	}
	return a.FromBigInt(bInv)
}

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
func (a BElem) AssertTypeFor(n RElem) RElem {
	return n.(BElem)
}
func (a LElem2NBigInt) AssertTypeFor(n RElem) RElem {
	return n.(LElem2NBigInt)
}
func (a LElem256) AssertTypeFor(n RElem) RElem {
	return n.(LElem256)
}
func (a LElem128) AssertTypeFor(n RElem) RElem {
	return n.(LElem128)
}
func (a LElem2N) AssertTypeFor(n RElem) RElem {
	return n.(LElem2N)
}
func (a LElemP) AssertTypeFor(n RElem) RElem {
	return n.(LElemP)
}
func (a SElemC) AssertTypeFor(n RElem) RElem {
	return n.(SElemC)
}
func (a SElemDS) AssertTypeFor(n RElem) RElem {
	return n.(SElemDS)
}
func (a BElem) FromInt(n int) RElem {
	return BElem(n)
}
func (a LElem2NBigInt) FromInt(n int) RElem {
	if n >= 0 {
		return LElem2NBigInt{new(big.Int).SetInt64(int64(n))}
	}
	return LElem2NBigInt{new(big.Int).Sub(LElem2NBigIntMod.Val, new(big.Int).SetInt64(-int64(n)))}
}
func (a LElem256) FromInt(n int) RElem {
	if n >= 0 {
		return LElem256{new(big.Int).SetInt64(int64(n))}
	}
	return LElem256{new(big.Int).Sub(LElem256Mod.Val, new(big.Int).SetInt64(-int64(n)))}
}
func (a LElem128) FromInt(n int) RElem {
	if n >= 0 {
		return LElem128(Uint128{0, 0}).Add(LElem128(Uint128{0, uint64(n)}))
	}
	return LElem128(Uint128{0, 0}).Sub(LElem128(Uint128{0, uint64(-n)}))
}
func (a LElem2N) FromInt(n int) RElem {
	return LElem2N(n)
}
func (a LElemP) FromInt(n int) RElem {
	if n >= 0 {
		return LElemP(0).Add(LElemP(n))
	}
	return LElemP(0).Sub(LElemP(-n))
}
func (a SElemC) FromInt(n int) RElem {
	if n >= 0 {
		return SElemC(0).Add(SElemC(n))
	}
	return SElemC(0).Sub(SElemC(-n))
}
func (a SElemDS) FromInt(n int) RElem {
	if n >= 0 {
		return SElemDS(0).Add(SElemDS(n))
	}
	return SElemDS(0).Sub(SElemDS(-n))
}
func (a BElem) FromUint64(n uint64) RElem {
	return BElem(n)
}
func (a LElem2NBigInt) FromUint64(n uint64) RElem {
	return LElem2NBigInt{new(big.Int).SetUint64(n)}
}
func (a LElem256) FromUint64(n uint64) RElem {
	return LElem256{new(big.Int).SetUint64(n)}
}
func (a LElem128) FromUint64(n uint64) RElem {
	if LElem128Mod.Hi == 0 {
		return LElem128(Uint128{0, n % LElem128Mod.Lo})
	}
	return LElem128(Uint128{0, n})
}
func (a LElem2N) FromUint64(n uint64) RElem {
	return LElem2N(n)
}
func (a LElemP) FromUint64(n uint64) RElem {
	return LElemP(n % LElemPMod64)
}
func (a SElemC) FromUint64(n uint64) RElem {
	return SElemC(n % uint64(SElemCMod))
}
func (a SElemDS) FromUint64(n uint64) RElem {
	return SElemDS(n % uint64(SElemDSMod))
}
func (a BElem) FromFloat64(n float64, fracBits int) RElem {
	panic("FromFloat64 for BElem undefined")
	return a
}
func (a LElem2NBigInt) FromFloat64(n float64, fracBits int) RElem {
	if n < 0 {
		return LElem2NBigInt{LElem2NBigIntZero.FromFloat64(-n, fracBits).(LElem2NBigInt).ToBigInt()}.Neg()
	} else {
		return LElem2NBigInt{LElem2NBigIntZero.FromFloat64(n, fracBits).(LElem2NBigInt).ToBigInt()}
	}
}
func (a LElem256) FromFloat64(n float64, fracBits int) RElem {
	if n < 0 {
		return LElem256{LElem128Zero.FromFloat64(-n, fracBits).(LElem128).ToBigInt()}.Neg()
	} else {
		return LElem256{LElem128Zero.FromFloat64(n, fracBits).(LElem128).ToBigInt()}
	}
}
func (a LElem128) FromFloat64(n float64, fracBits int) RElem {
	sgn := 1
	if n < 0 {
		sgn = -1
		n = -n
	}

	const shift64 = 1 << 64
	x := math.Round(n * float64(uint64(1)<<fracBits))
	hi := uint64(x / shift64)
	lo := uint64(x - float64(hi)*shift64)
	out := LElem128{hi, lo}
	if sgn < 0 {
		return out.Neg()
	}
	return out
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
func (a SElemC) FromFloat64(n float64, fracBits int) RElem {
	panic("SElemC is not meant for storing fractional values")
	return a
}
func (a SElemDS) FromFloat64(n float64, fracBits int) RElem {
	panic("SElemDS is not meant for storing fractional values")
	return a
}
func (a BElem) Float64(fracBits int) float64 {
	panic("Float64 for BElem undefined")
	return 0
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
func (a LElem128) Float64(fracBits int) float64 {
	var sgn int
	var b LElem128
	if lessThan128(Uint128(a), Uint128(LElem128ModHalf)) {
		sgn, b = 1, a
	} else {
		sgn, b = -1, a.Neg().(LElem128)
	}

	shift := float64(uint64(1) << fracBits)
	const shift64 = 1 << 64
	return float64(sgn) * (((float64(b.Hi) / shift) * shift64) + (float64(b.Lo) / shift))
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
func (a SElemC) Float64(fracBits int) float64 {
	panic("SElemC is not meant for storing fractional values")
	return 0
}
func (a SElemDS) Float64(fracBits int) float64 {
	panic("SElemDS is not meant for storing fractional values")
	return 0
}
func (a BElem) FromBytes(buf []byte) RElem {
	return BElem(binary.LittleEndian.Uint64(buf))
}
func (a LElem2NBigInt) FromBytes(buf []byte) RElem {
	return LElem2NBigInt{new(big.Int).SetBytes(buf[:a.NumBytes()])}
}
func (a LElem256) FromBytes(buf []byte) RElem {
	return LElem256{new(big.Int).SetBytes(buf[:a.NumBytes()])}
}
func (a LElem128) FromBytes(buf []byte) RElem {
	return LElem128(Uint128{binary.LittleEndian.Uint64(buf), binary.LittleEndian.Uint64(buf[8:])})
}
func (a LElem2N) FromBytes(buf []byte) RElem {
	return LElem2N(binary.LittleEndian.Uint64(buf))
}
func (a LElemP) FromBytes(buf []byte) RElem {
	return LElemP(binary.LittleEndian.Uint64(buf))
}
func (a SElemC) FromBytes(buf []byte) RElem {
	return SElemC(buf[0])
}
func (a SElemDS) FromBytes(buf []byte) RElem {
	return SElemDS(buf[0])
}
func (a BElem) ToBytes(buf []byte) {
	binary.LittleEndian.PutUint64(buf, uint64(a))
}
func (a LElem2NBigInt) ToBytes(buf []byte) {
	if a.Val.Sign() < 0 {
		a.Val.Mod(a.Val, LElem2NBigIntMod.Val)
	}
	a.Val.FillBytes(buf[:a.NumBytes()])
}
func (a LElem256) ToBytes(buf []byte) {
	if a.Val.Sign() < 0 {
		a.Val.Mod(a.Val, LElem256Mod.Val)
	}
	a.Val.FillBytes(buf[:a.NumBytes()])
}
func (a LElem128) ToBytes(buf []byte) {
	binary.LittleEndian.PutUint64(buf, a.Hi)
	binary.LittleEndian.PutUint64(buf[8:], a.Lo)
}
func (a LElem2N) ToBytes(buf []byte) {
	binary.LittleEndian.PutUint64(buf, uint64(a))
}
func (a LElemP) ToBytes(buf []byte) {
	binary.LittleEndian.PutUint64(buf, uint64(a))
}
func (a SElemC) ToBytes(buf []byte) {
	buf[0] = byte(a)
}
func (a SElemDS) ToBytes(buf []byte) {
	buf[0] = byte(a)
}
func (a BElem) Uint64() uint64 {
	return uint64(a)
}
func (a LElem2NBigInt) Uint64() uint64 {
	return a.Val.Uint64()
}
func (a LElem256) Uint64() uint64 {
	return a.Val.Uint64()
}
func (a LElem128) Uint64() uint64 {
	return a.Lo
}
func (a LElem2N) Uint64() uint64 {
	return uint64(a)
}
func (a LElemP) Uint64() uint64 {
	return uint64(a)
}
func (a SElemC) Uint64() uint64 {
	return uint64(a)
}
func (a SElemDS) Uint64() uint64 {
	return uint64(a)
}
func (a BElem) Zero() RElem {
	return BElem(0)
}
func (a LElem2NBigInt) Zero() RElem {
	return LElem2NBigInt{big.NewInt(0)}
}
func (a LElem256) Zero() RElem {
	return LElem256{big.NewInt(0)}
}
func (a LElem128) Zero() RElem {
	return LElem128(Uint128{0, 0})
}
func (a LElem2N) Zero() RElem {
	return LElem2N(0)
}
func (a LElemP) Zero() RElem {
	return LElemP(0)
}
func (a SElemC) Zero() RElem {
	return SElemC(0)
}
func (a SElemDS) Zero() RElem {
	return SElemDS(0)
}
func (a BElem) One() RElem {
	return BElem(^uint64(0))
}
func (a LElem2NBigInt) One() RElem {
	return LElem2NBigInt{big.NewInt(1)}
}
func (a LElem256) One() RElem {
	return LElem256{big.NewInt(1)}
}
func (a LElem128) One() RElem {
	return LElem128(Uint128{0, 1})
}
func (a LElem2N) One() RElem {
	return LElem2N(1)
}
func (a LElemP) One() RElem {
	return LElemP(1)
}
func (a SElemC) One() RElem {
	return SElemC(1)
}
func (a SElemDS) One() RElem {
	return SElemDS(1)
}
func (a BElem) Modulus() *big.Int {
	panic("Modulus for BElem undefined")
	return new(big.Int)
}
func (a LElem2NBigInt) Modulus() *big.Int {
	return new(big.Int).Set(LElem2NBigIntMod.Val)
}
func (a LElem256) Modulus() *big.Int {
	return new(big.Int).Set(LElem256Mod.Val)
}
func (a LElem128) Modulus() *big.Int {
	return new(big.Int).Set(LElem128ModBig)
}
func (a LElem2N) Modulus() *big.Int {
	return new(big.Int).Set(LElem2NModBig)
}
func (a LElemP) Modulus() *big.Int {
	return new(big.Int).Set(LElemPModBig)
}
func (a SElemC) Modulus() *big.Int {
	return new(big.Int).Set(SElemCModBig)
}
func (a SElemDS) Modulus() *big.Int {
	return new(big.Int).Set(SElemDSModBig)
}
func (a BElem) NumBytes() uint32 {
	return BElemBytes
}
func (a LElem2NBigInt) NumBytes() uint32 {
	return 1 + ((uint32(LElem2NBigIntModBitLen) - 2) >> 3)
}
func (a LElem256) NumBytes() uint32 {
	return LElem256Bytes
}
func (a LElem128) NumBytes() uint32 {
	return LElem128Bytes
}
func (a LElem2N) NumBytes() uint32 {
	return LElem2NBytes
}
func (a LElemP) NumBytes() uint32 {
	return LElemPBytes
}
func (a SElemC) NumBytes() uint32 {
	return SElemCBytes
}
func (a SElemDS) NumBytes() uint32 {
	return SElemDSBytes
}
func (a BElem) ModBitLength() int {
	panic("ModBitLength for BElem undefined")
	return 0
}
func (a LElem2NBigInt) ModBitLength() int {
	return LElem2NBigIntModBitLen
}
func (a LElem256) ModBitLength() int {
	return LElem256ModBitLen
}
func (a LElem128) ModBitLength() int {
	return LElem128ModBitLen
}
func (a LElem2N) ModBitLength() int {
	return LElem2NModBitLen
}
func (a LElemP) ModBitLength() int {
	return LElemPModBitLen
}
func (a SElemC) ModBitLength() int {
	return SElemCModBitLen
}
func (a SElemDS) ModBitLength() int {
	return SElemDSModBitLen
}
func (a BElem) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB >= 8*int(a.NumBytes()) {
		panic("Invalid bit position")
	}
	mask := uint64(1) << posFromLSB
	return boolToUint((uint64(a) & mask) > 0)
}
func (a LElem2NBigInt) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB >= a.ModBitLength() {
		panic("Invalid bit position")
	}
	return a.Val.Bit(posFromLSB)
}
func (a LElem256) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB >= a.ModBitLength() {
		panic("Invalid bit position")
	}
	return a.Val.Bit(posFromLSB)
}
func (a LElem128) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB >= a.ModBitLength() {
		panic("Invalid bit position")
	}
	if posFromLSB < 64 {
		return boolToUint((a.Lo & (uint64(1) << posFromLSB)) > 0)
	}
	return boolToUint((a.Hi & (uint64(1) << (posFromLSB - 64))) > 0)
}
func (a LElem2N) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB >= a.ModBitLength() {
		panic("Invalid bit position")
	}
	mask := LElem2N(1) << posFromLSB
	return boolToUint((a & mask) > 0)
}
func (a LElemP) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB > a.ModBitLength() {
		panic("Invalid bit position")
	}
	return boolToUint((a & (LElemP(1) << posFromLSB)) > 0)
}
func (a SElemC) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB > a.ModBitLength() {
		panic("Invalid bit position")
	}
	return boolToUint((a & (SElemC(1) << posFromLSB)) > 0)
}
func (a SElemDS) GetBit(posFromLSB int) uint {
	if posFromLSB < 0 || posFromLSB > a.ModBitLength() {
		panic("Invalid bit position")
	}
	return boolToUint((a & (SElemDS(1) << posFromLSB)) > 0)
}
func boolToUint(b bool) uint {
	if b {
		return 1
	}
	return 0
}
func (a BElem) Trunc(nBits int) RElem {
	if nBits < 0 || nBits >= 8*int(a.NumBytes()) {
		panic("Invalid number of bits")
	} else if nBits == 8*int(a.NumBytes()) {
		return a
	}
	return BElem(uint64(a) % (uint64(1) << nBits))
}
func (a LElem2NBigInt) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return LElem2NBigInt{new(big.Int).Set(a.Val)}
	}

	return LElem256{new(big.Int).Rem(a.Val, new(big.Int).Lsh(big.NewInt(1), uint(nBits)))}
}
func (a LElem256) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return LElem256{new(big.Int).Set(a.Val)}
	}

	return LElem256{new(big.Int).Rem(a.Val, new(big.Int).Lsh(big.NewInt(1), uint(nBits)))}
}
func (a LElem128) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return a
	}

	if nBits > 64 {
		return LElem128{a.Hi % (uint64(1) << (nBits - 64)), a.Lo}
	} else if nBits == 64 {
		return LElem128{0, a.Lo}
	}
	return LElem128{0, a.Lo % (uint64(1) << nBits)}
}
func (a LElem2N) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return a
	}
	return a % (LElem2N(1) << nBits)
}
func (a LElemP) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return a
	}
	return a % (LElemP(1) << nBits)
}
func (a SElemC) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return a
	}
	return a % (SElemC(1) << nBits)
}
func (a SElemDS) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return a
	}
	return a % (SElemDS(1) << nBits)
}
func (a BElem) TypeID() uint8 {
	return BElemUniqueID
}
func (a LElem2NBigInt) TypeID() uint8 {
	return LElem2NBigIntUniqueID
}
func (a LElem256) TypeID() uint8 {
	return LElem256UniqueID
}
func (a LElem128) TypeID() uint8 {
	return LElem128UniqueID
}
func (a LElem2N) TypeID() uint8 {
	return LElem2NUniqueID
}
func (a LElemP) TypeID() uint8 {
	return LElemPUniqueID
}
func (a SElemC) TypeID() uint8 {
	return SElemCUniqueID
}
func (a SElemDS) TypeID() uint8 {
	return SElemDSUniqueID
}

func randBytes(rtype RElem, buf []byte, prg *frand.RNG) RElem {
	prg.Read(buf)
	return rtype.FromBytes(buf)
}
func (a BElem) Rand(prg *frand.RNG) RElem {
	buf := make([]byte, a.NumBytes())
	return randBytes(a, buf, prg)
}
func (a LElem2N) Rand(prg *frand.RNG) RElem {
	buf := make([]byte, a.NumBytes())
	return randBytes(a, buf, prg)
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
func (a LElem128) Rand(prg *frand.RNG) RElem {
	buf := make([]byte, a.NumBytes())
again:
	r := randBytes(a, buf, prg).(LElem128)
	if !lessThan128(Uint128(r), Uint128(LElem128RandBnd)) {
		goto again
	}
	_, rem := div128(Uint128{0, 0}, Uint128(r), Uint128(LElem128Mod))
	return LElem128(rem)
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
func (a SElemC) Rand(prg *frand.RNG) RElem {
	buf := make([]byte, a.NumBytes())
again:
	r := randBytes(a, buf, prg).(SElemC)
	if r >= SElemCRandBnd {
		goto again
	}
	return r % SElemCMod
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
func (a BElem) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits > 8*int(a.NumBytes()) {
		panic("Requested bit length is larger than BElem size")
	}
	return a.Rand(prg).Trunc(nbits)
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
func (a LElem2NBigInt) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits >= a.ModBitLength() {
		panic("Requested bit length is larger than modulus")
	}

	buf := make([]byte, 1+((nbits-1)/8))
	prg.Read(buf)

	return LElem2NBigInt{new(big.Int).Rem(new(big.Int).SetBytes(buf), new(big.Int).Lsh(big.NewInt(1), uint(nbits)))}
}
func (a LElem256) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits >= a.ModBitLength() {
		panic("Requested bit length is larger than modulus")
	}

	buf := make([]byte, 1+((nbits-1)/8))
	prg.Read(buf)

	return LElem256{new(big.Int).Rem(new(big.Int).SetBytes(buf), new(big.Int).Lsh(big.NewInt(1), uint(nbits)))}
}
func (a LElem128) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits >= a.ModBitLength() {
		panic("Requested bit length is larger than modulus")
	}
	buf := make([]byte, a.NumBytes())
	out := randBytes(a, buf, prg).(LElem128)
	if nbits > 64 {
		out.Hi = out.Hi % (uint64(1) << (nbits - 64))
	} else if nbits == 64 {
		out.Hi = 0
	} else {
		out.Hi = 0
		out.Lo = out.Lo % (uint64(1) << nbits)
	}
	return out
}
func (a LElemP) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits >= a.ModBitLength() {
		panic("Requested bit length is larger than modulus")
	}
	buf := make([]byte, a.NumBytes())
	return randBytes(a, buf, prg).(LElemP) % (LElemP(1) << nbits)
}
func (a SElemC) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits >= a.ModBitLength() {
		panic("Requested bit length is larger than modulus")
	}
	buf := make([]byte, a.NumBytes())
	return randBytes(a, buf, prg).(SElemC) % (SElemC(1) << nbits)
}
func (a SElemDS) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits >= a.ModBitLength() {
		panic("Requested bit length is larger than modulus")
	}
	buf := make([]byte, a.NumBytes())
	return randBytes(a, buf, prg).(SElemDS) % (SElemDS(1) << nbits)
}

func (a BElem) FromBigInt(n *big.Int) RElem {
	panic("Not implemented")
	return nil
}
func (a LElem2N) FromBigInt(n *big.Int) RElem {
	panic("Not implemented")
	return nil
}
func (a LElemP) FromBigInt(n *big.Int) RElem {
	panic("Not implemented")
	return nil
}
func (a SElemC) FromBigInt(n *big.Int) RElem {
	panic("Not implemented")
	return nil
}
func (a SElemDS) FromBigInt(n *big.Int) RElem {
	panic("Not implemented")
	return nil
}
func (a BElem) Copy() RElem {
	return a
}
func (a LElem2NBigInt) Copy() RElem {
	return LElem2NBigInt{new(big.Int).Set(a.Val)}
}
func (a LElem256) Copy() RElem {
	return LElem256{new(big.Int).Set(a.Val)}
}
func (a LElem128) Copy() RElem {
	return a
}
func (a LElem2N) Copy() RElem {
	return a
}
func (a LElemP) Copy() RElem {
	return a
}
func (a SElemC) Copy() RElem {
	return a
}
func (a SElemDS) Copy() RElem {
	return a
}

//Initalizes a RElemMatrix of size n RElems
func InitRVec(val RElem, n int) RVec {
	res := make(RVec, n)
	for i := range res {
		res[i] = val.Copy()
	}
	return res
}

//Initalizes a RElemMatrix with nrows = len(RMat) and ncols = len(RMat[0])
func InitRMat(val RElem, nrow, ncol int) RMat {
	res := make(RMat, nrow)
	for i := range res {
		res[i] = InitRVec(val, ncol)
	}
	return res
}

func RMultMatVec(a RMat, b RVec) RVec {
	rtype := a.Type()
	r, c := a.Dims()
	n := len(b)
	if c != n {
		panic("Dimensions do not match")
	}

	out := InitRVec(rtype.Zero(), r)
	for i := range out {
		for j := range b {
			out[i] = out[i].Add(a[i][j].Mul(b[j]))
		}
	}
	return out
}

func RMultMat(a, b RMat) RMat {
	rtype := a.Type()
	r1, c1 := a.Dims()
	r2, c2 := b.Dims()
	if c1 != r2 {
		panic("Dimensions do not match")
	}

	out := InitRMat(rtype.Zero(), r1, c2)
	for i := range out {
		for j := range out[i] {
			for k := range a[i] {
				out[i][j] = out[i][j].Add(a[i][k].Mul(b[k][j]))
			}
		}
	}

	return out
}

func RMultElemMat(a, b RMat) RMat {
	rtype := a.Type()
	r1, c1 := a.Dims()
	r2, c2 := b.Dims()
	if r1 != r2 || c1 != c2 {
		panic("Dimensions do not match")
	}

	out := InitRMat(rtype.Zero(), r1, c1)
	for i := range out {
		for j := range out[i] {
			out[i][j] = a[i][j].Mul(b[i][j])
		}
	}
	return out
}

//MultVec computes elem-wise product of a and b
func RMultVec(a, b RVec) RVec {
	c := make(RVec, len(a))
	for i := range c {
		c[i] = a[i].Mul(b[i])
	}
	return c
}

//MultMat computes the elem-wise product of a and b

//MultConstVec multiplies c with each element in v (and performs appr modulus)
func RMultConstVec(c RElem, v RVec) RVec {
	res := make(RVec, len(v))
	for i := range res {
		res[i] = c.Mul(v[i])
	}
	return res

}

//MultConstMat multiplies c with each element in m (and performs appr modulus)
func RMultConstMat(c RElem, m RMat) RMat {
	res := make(RMat, (len(m)))
	for i := range m {
		res[i] = make(RVec, len(m[0]))
		for j := range m[i] {
			res[i][j] = c.Mul(m[i][j])

		}
	}
	return res
}

func (a RMat) Clear() {
	for i := range a {
		for j := range a[i] {
			a[i][j] = a[i][j].Zero()
		}
	}
}

func (a RVec) Clear() {
	for i := range a {
		a[i] = a[i].Zero()
	}
}

func (a RVec) Add(b RVec) {
	if len(a) != len(b) {
		panic("Inconsistent vector lengths")
	}

	for i := range a {
		a[i] = a[i].Add(b[i])
	}
}

func (a RMat) Add(b RMat) {
	if len(a) != len(b) || len(a[0]) != len(b[0]) {
		panic("Inconsistent dimensions")
	}

	for i := range a {
		for j := range a[i] {
			a[i][j] = a[i][j].Add(b[i][j])
		}
	}
}

func (a RVec) AddScalar(b RElem) {
	for i := range a {
		a[i] = a[i].Add(b)
	}
}

func (a RMat) AddScalar(b RElem) {
	for i := range a {
		for j := range a[i] {
			a[i][j] = a[i][j].Add(b)
		}
	}
}

func (a RVec) MulElem(b RVec) {
	if len(a) != len(b) {
		panic("Inconsistent vector lengths")
	}

	for i := range a {
		a[i] = a[i].Mul(b[i])
	}
}

func (a RMat) MulElem(b RMat) {
	if len(a) != len(b) || len(a[0]) != len(b[0]) {
		panic("Inconsistent dimensions")
	}

	for i := range a {
		for j := range a[i] {
			a[i][j] = a[i][j].Mul(b[i][j])
		}
	}
}

func (a RVec) MulScalar(b RElem) {
	for i := range a {
		a[i] = a[i].Mul(b)
	}
}

func (a RMat) MulScalar(b RElem) {
	for i := range a {
		for j := range a[i] {
			a[i][j] = a[i][j].Mul(b)
		}
	}
}

func (a RVec) Sub(b interface{}) {
	c := b.(RVec)
	if len(a) != len(c) {
		panic("Inconsistent vector lengths")
	}

	for i := range a {
		a[i] = a[i].Sub(c[i])
	}
}

func (a RMat) Sub(b RMat) {
	if len(a) != len(b) || len(a[0]) != len(b[0]) {
		panic("Inconsistent dimensions")
	}

	for i := range a {
		for j := range a[i] {
			a[i][j] = a[i][j].Sub(b[i][j])
		}
	}
}

// axis: 0 (row sum) or 1 (column sum)
func (a RMat) Sum(axis int) RVec {
	if axis != 0 && axis != 1 {
		panic("axis parameter must be 0 (row sum) or 1 (column sum)")
	}

	r, c := a.Dims()
	rtype := a.Type()

	var out RVec
	if axis == 0 {
		out = InitRVec(rtype.Zero(), r)
	} else {
		out = InitRVec(rtype.Zero(), c)
	}

	for i := range a {
		for j := range a[i] {
			if axis == 0 {
				out[i] = out[i].Add(a[i][j])
			} else {
				out[j] = out[j].Add(a[i][j])
			}
		}
	}

	return out
}

func (a RMat) Dims() (int, int) {
	return len(a), len(a[0])
}

func (a RMat) Type() RElem {
	return a[0][0]
}

func (a RVec) Type() RElem {
	return a[0]
}

func (a RMat) Copy() RMat {
	r, c := a.Dims()
	b := make(RMat, r)
	for i := range a {
		b[i] = make(RVec, c)
		for j := range a[i] {
			b[i][j] = a[i][j].Copy()
		}
	}
	return b
}

func (a RVec) Copy() RVec {
	n := len(a)
	b := make(RVec, n)
	for i := range a {
		b[i] = a[i].Copy()
	}
	return b
}

func (a RMat) Transpose() RMat {
	r, c := a.Dims()
	b := make(RMat, c)
	for i := range b {
		b[i] = make(RVec, r)
		for j := range b[i] {
			b[i][j] = a[j][i].Copy()
		}
	}
	return b
}

func (a RMat) NumBytes() uint32 {
	r, c := a.Dims()
	return a.Type().NumBytes() * uint32(r) * uint32(c)
}

func (a RVec) NumBytes() uint32 {
	return a.Type().NumBytes() * uint32(len(a))
}

func (a RMat) MarshalBinary() ([]byte, error) {
	buf := make([]byte, a.NumBytes())
	offset := uint64(0)
	byteSize := uint64(a.Type().NumBytes())
	for i := range a {
		for j := range a[i] {
			a[i][j].ToBytes(buf[offset:])
			offset += byteSize
		}
	}
	return buf, nil
}

func (a *RMat) UnmarshalBinary(data []byte) error {
	rtype := a.Type()
	offset := uint64(0)
	byteSize := uint64(a.Type().NumBytes())
	r, c := a.Dims()
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			(*a)[i][j] = rtype.FromBytes(data[offset:])
			offset += byteSize
		}
	}
	return nil
}

func (a RVec) MarshalBinary() ([]byte, error) {
	buf := make([]byte, a.NumBytes())
	offset := uint64(0)
	byteSize := uint64(a.Type().NumBytes())
	for i := range a {
		a[i].ToBytes(buf[offset:])
		offset += byteSize
	}
	return buf, nil
}

func (a *RVec) UnmarshalBinary(data []byte) error {
	rtype := a.Type()
	offset := uint64(0)
	byteSize := uint64(a.Type().NumBytes())
	n := len(*a)
	for i := 0; i < n; i++ {
		(*a)[i] = rtype.FromBytes(data[offset:])
		offset += byteSize
	}
	return nil
}

func (a RVec) Trunc(nBits int) {
	for i := range a {
		a[i] = a[i].Trunc(nBits)
	}
}

func (a RMat) Trunc(nBits int) {
	for i := range a {
		for j := range a[i] {
			a[i][j] = a[i][j].Trunc(nBits)
		}
	}
}

func (a RVec) ToFloat(fracBits int) []float64 {
	out := make([]float64, len(a))
	for i := range out {
		out[i] = a[i].Float64(fracBits)
	}
	return out
}

func (a RMat) ToFloat(fracBits int) [][]float64 {
	out := make([][]float64, len(a))
	for i := range out {
		out[i] = make([]float64, len(a[i]))
		for j := range out[i] {
			out[i][j] = a[i][j].Float64(fracBits)
		}
	}
	return out
}

func FloatToRVec(rtype RElem, a []float64, fracBits int) RVec {
	out := make(RVec, len(a))
	for i := range out {
		out[i] = rtype.FromFloat64(a[i], fracBits)
	}
	return out
}

func FloatToRMat(rtype RElem, a [][]float64, fracBits int) RMat {
	out := make(RMat, len(a))
	for i := range out {
		out[i] = make(RVec, len(a[i]))
		for j := range out[i] {
			out[i][j] = rtype.FromFloat64(a[i][j], fracBits)
		}
	}
	return out
}

func IntToRVec(rtype RElem, a []int) RVec {
	out := make(RVec, len(a))
	for i := range out {
		out[i] = rtype.FromInt(a[i])
	}
	return out
}

func IntToRMat(rtype RElem, a [][]int) RMat {
	out := make(RMat, len(a))
	for i := range out {
		out[i] = make(RVec, len(a[i]))
		for j := range out[i] {
			out[i][j] = rtype.FromInt(a[i][j])
		}
	}
	return out
}

func (a RVec) ToInt() []int {
	res := make([]int, len(a))
	for i := range res {
		res[i] = int(a[i].Uint64())
	}

	return res
}

func (a RMat) ToInt() [][]int {
	out := make([][]int, len(a))
	for i := range out {
		out[i] = make([]int, len(a[i]))
		for j := range out[i] {
			out[i][j] = int(a[i][j].Uint64())
		}
	}
	return out

}

func (a RVec) Filter(filt []bool) RVec {
	if len(a) != len(filt) {
		panic("Filter length does not match vector")
	}

	count := 0
	for i := range filt {
		if filt[i] {
			count++
		}
	}

	out := make(RVec, count)
	index := 0
	for i := range filt {
		if filt[i] {
			out[index] = a[i]
			index++
		}
	}
	return out
}
