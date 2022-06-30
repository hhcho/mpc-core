package mpc_core

import (
	"encoding/binary"
	"math"
	"math/big"
	"math/bits"

	"github.com/hhcho/frand"
)

type Uint128 struct {
	Hi, Lo uint64
}

type LElem128 Uint128

const LElem128Bytes uint32 = 16

func initLElem128RandBnd() LElem128 {
	mOne := Uint128{^uint64(0), ^uint64(0)}
	_, rem := div128(Uint128{0, 0}, mOne, Uint128(LElem128Mod))
	out, _ := sub128(mOne, rem, 0)
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

func equal128(a, b Uint128) bool {
	return (a.Hi == b.Hi) && (a.Lo == b.Lo)
}

func add128(a, b Uint128, carryin uint64) (c Uint128, carry uint64) {
	c.Lo, carry = bits.Add64(a.Lo, b.Lo, carryin)
	c.Hi, carry = bits.Add64(a.Hi, b.Hi, carry)
	return
}

func sub128(a, b Uint128, borrowin uint64) (c Uint128, borrow uint64) {
	c.Lo, borrow = bits.Sub64(a.Lo, b.Lo, borrowin)
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
	rhat, _ := sub128(un32, t, 0) //rhat := un32 - q1*yn1

	_, t = mul128(q1, yn0)
	for q1.Hi != 0 || lessThan128(Uint128{rhat.Lo, un1.Lo}, t) {
		//for q1 >= two32 || q1*yn0 > two32*rhat+un1 {
		q1, _ = sub128(q1, Uint128{0, 1}, 0) //	q1--
		rhat, _ = add128(rhat, yn1, 0)       //rhat += yn1
		if rhat.Hi != 0 {                    //if rhat >= two32 {
			break
		}
	}

	_, q1y := mul128(q1, y)
	un21, _ := sub128(Uint128{un32.Lo, un1.Lo}, q1y, 0) //un21 := un32*two32 + un1 - q1*y
	q0hi, r := bits.Div64(0, un21.Hi, yn1.Lo)
	q0lo, _ := bits.Div64(r, un21.Lo, yn1.Lo) //q0 := un21 / yn1
	q0 := Uint128{q0hi, q0lo}
	_, q0yn1 := mul128(q0, yn1)
	rhat, _ = sub128(un21, q0yn1, 0) //rhat = un21 - q0*yn1

	_, t = mul128(q0, yn0)
	for q0.Hi != 0 || lessThan128(Uint128{rhat.Lo, un0.Lo}, t) {
		//for q0 >= two32 || q0*yn0 > two32*rhat+un0 {
		q0, _ = sub128(q0, Uint128{0, 1}, 0) //	q0--
		rhat, _ = add128(rhat, yn1, 0)       //rhat += yn1
		if rhat.Hi != 0 {                    //if rhat >= two32 {
			break
		}
	}

	quo = Uint128{q1.Lo, q0.Lo}

	_, q0y := mul128(q0, y)
	rem, _ = sub128(Uint128{un21.Lo, un0.Lo}, q0y, 0)
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
	_, carry := add128(mNLo, lo, 0)
	t, carry1 := add128(mNHi, hi, 0)
	t, carry2 := add128(t, Uint128{0, carry}, 0)

	if carry1+carry2 > 0 || !lessThan128(t, Uint128(LElem128Mod)) {
		res, _ := sub128(t, Uint128(LElem128Mod), 0)
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
	out, carry := add128(Uint128(a), Uint128(b.(LElem128)), 0)
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
	out, _ := sub128(Uint128(LElem128Mod), Uint128(a), 0)
	return LElem128(out)
}

func (a LElem128) Inv() RElem {
	bInv := big.NewInt(0).ModInverse(a.ToBigInt(), LElem128ModBig)
	if bInv == nil {
		panic("ModInverse does not exist")
	}
	return a.FromBigInt(bInv)
}

func (a LElem128) AssertTypeFor(n RElem) RElem {
	return n.(LElem128)
}
func (a LElem128) FromInt(n int) RElem {
	if n >= 0 {
		return LElem128(Uint128{0, 0}).Add(LElem128(Uint128{0, uint64(n)}))
	}
	return LElem128(Uint128{0, 0}).Sub(LElem128(Uint128{0, uint64(-n)}))
}
func (a LElem128) FromUint64(n uint64) RElem {
	if LElem128Mod.Hi == 0 {
		return LElem128(Uint128{0, n % LElem128Mod.Lo})
	}
	return LElem128(Uint128{0, n})
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

func (a LElem128) FromBytes(buf []byte) RElem {
	return LElem128(Uint128{binary.LittleEndian.Uint64(buf), binary.LittleEndian.Uint64(buf[8:])})
}

func frombytes128(buf []byte) Uint128 {
	return Uint128{binary.LittleEndian.Uint64(buf), binary.LittleEndian.Uint64(buf[8:])}
}

func (a LElem128) ToBytes(buf []byte) {
	binary.LittleEndian.PutUint64(buf, a.Hi)
	binary.LittleEndian.PutUint64(buf[8:], a.Lo)
}
func (a LElem128) Zero() RElem {
	return LElem128(Uint128{0, 0})
}
func (a LElem128) NumBytes() uint32 {
	return LElem128Bytes
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
func (a LElem128) TypeID() uint8 {
	return LElem128UniqueID
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
func (a LElem128) Copy() RElem {
	return a
}
func (a LElem128) Uint64() uint64 {
	return a.Lo
}
func (a LElem128) One() RElem {
	return LElem128(Uint128{0, 1})
}
func (a LElem128) Modulus() *big.Int {
	return new(big.Int).Set(LElem128ModBig)
}
func (a LElem128) ModBitLength() int {
	return LElem128ModBitLen
}
