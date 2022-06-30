package mpc_core

import (
	"math/big"
)

type Uint256 struct {
	Hi, Lo Uint128
}

type LElem256 Uint256

const LElem256Bytes uint32 = 32

func initLElem256RandBnd() LElem256 {
	//TODO

	return LElem256{}
}

var LElem256RandBnd LElem256 = initLElem256RandBnd()
var LElem256Zero LElem256 = LElem256{Uint128{0, 0}, Uint128{0, 0}}

//TODO
var LElem256Mod = LElem256{} // 2^256 - 189
//TODO
var LElem256ModHalf = LElem256{} // 2^255 - 94
//TODO
var LElem256ModBitLen int = 0

const LElem256UniqueID uint8 = 5

//TODO
var LElem256ModBig *big.Int = nil

//TODO
var LElem256ModHalfBig *big.Int = nil

func mul256(a, b Uint256) (Hi, Lo Uint256) {
	h00, r0 := mul128(a.Lo, b.Lo)
	h01, l01 := mul128(a.Lo, b.Hi)
	h10, l10 := mul128(a.Hi, b.Lo)
	h11, l11 := mul128(a.Hi, b.Hi)

	t, carry1 := add128(h00, l01, 0)
	r1, carry2 := add128(t, l10, 0)
	t, carry1 = add128(l11, h01, carry1)
	r2, carry2 := add128(t, h10, carry2)
	r3, _ := add128(h11, Uint128{0, 0}, carry1+carry2)

	Hi = Uint256{r3, r2}
	Lo = Uint256{r1, r0}
	return
}

func leadingZeros256(x Uint256) int {
	return leadingZeros128(x.Hi) + leadingZeros128(x.Lo)
}

// Lsh returns u<<n.
func lsh256(u Uint256, n uint) (s Uint256) {
	if n > 128 {
		s.Lo = Uint128{0, 0}
		s.Hi = lsh128(u.Lo, n-128)
	} else {
		s.Lo = lsh128(u.Lo, n)
		s.Hi = or128(lsh128(u.Hi, n), rsh128(u.Lo, 128-n))
	}
	return
}

// Rsh returns u>>n.
func rsh256(u Uint256, n uint) (s Uint256) {
	if n > 128 {
		s.Lo = rsh128(u.Hi, n-128)
		s.Hi = Uint128{0, 0}
	} else {
		s.Lo = or128(rsh128(u.Lo, n), lsh128(u.Hi, 64-n))
		s.Hi = rsh128(u.Hi, n)

	}
	return
}

func or256(a, b Uint256) (c Uint256) {
	c.Hi = or128(a.Hi, b.Hi)
	c.Lo = or128(a.Lo, b.Lo)
	return
}

func lessThan256(a, b Uint256) bool {
	if equal128(a.Hi, b.Hi) {
		return lessThan128(a.Lo, b.Lo)
	}
	return lessThan128(a.Hi, b.Hi)
}

func add256(a, b Uint256, carryin uint64) (c Uint256, carry uint64) {
	c.Lo, carry = add128(a.Lo, b.Lo, carryin)
	c.Hi, carry = add128(a.Hi, b.Hi, carry)
	return
}

func sub256(a, b Uint256, borrowin uint64) (c Uint256, borrow uint64) {
	c.Lo, borrow = sub128(a.Lo, b.Lo, borrowin)
	c.Hi, borrow = sub128(a.Hi, b.Hi, borrow)
	return
}

func div256(Hi, Lo, y Uint256) (quo, rem Uint256) {
	if equal128(y.Hi, Uint128{0, 0}) && equal128(y.Lo, Uint128{0, 0}) {
		panic("Divide by zero error")
	}
	if !lessThan256(Hi, y) {
		panic("Overflow error")
	}
	s := uint(leadingZeros256(y))
	y = lsh256(y, s)
	yn1 := Uint256{Uint128{0, 0}, y.Hi}
	yn0 := Uint256{Uint128{0, 0}, y.Lo}
	un32 := or256(lsh256(Hi, s), rsh256(Lo, 256-s))
	un10 := lsh256(Lo, s)
	un1 := Uint256{Uint128{0, 0}, un10.Hi}
	un0 := Uint256{Uint128{0, 0}, un10.Lo}
	q1hi, r := div128(Uint128{0, 0}, un32.Hi, yn1.Lo)
	q1lo, _ := div128(r, un32.Lo, yn1.Lo)
	q1 := Uint256{q1hi, q1lo}
	_, t := mul256(q1, yn1)
	rhat, _ := sub256(un32, t, 0) //rhat := un32 - q1*yn1

	_, t = mul256(q1, yn0)
	for !equal128(q1.Hi, Uint128{0, 0}) || lessThan256(Uint256{rhat.Lo, un1.Lo}, t) {
		//for q1 >= two32 || q1*yn0 > two32*rhat+un1 {
		q1, _ = sub256(q1, uint64To256(1), 0)  //	q1--
		rhat, _ = add256(rhat, yn1, 0)         //rhat += yn1
		if !equal128(rhat.Hi, Uint128{0, 0}) { //if rhat >= two32 {
			break
		}
	}
	_, q1y := mul256(q1, y)
	un21, _ := sub256(Uint256{un32.Lo, un1.Lo}, q1y, 0) //un21 := un32*two32 + un1 - q1*y
	q0hi, r := div128(Uint128{0, 0}, un21.Hi, yn1.Lo)
	q0lo, _ := div128(r, un21.Lo, yn1.Lo) //q0 := un21 / yn1
	q0 := Uint256{q0hi, q0lo}
	_, q0yn1 := mul256(q0, yn1)
	rhat, _ = sub256(un21, q0yn1, 0) //rhat = un21 - q0*yn1

	_, t = mul256(q0, yn0)
	for !equal128(q0.Hi, Uint128{0, 0}) || lessThan256(Uint256{rhat.Lo, un0.Lo}, t) {
		//for q0 >= two32 || q0*yn0 > two32*rhat+un0 {
		q0, _ = sub256(q0, uint64To256(1), 0)  //	q0--
		rhat, _ = add256(rhat, yn1, 0)         //rhat += yn1
		if !equal128(rhat.Hi, Uint128{0, 0}) { //if rhat >= two32 {
			break
		}
	}
	quo = Uint256{q1.Lo, q0.Lo}

	_, q0y := mul256(q0, y)
	rem, _ = sub256(Uint256{un21.Lo, un0.Lo}, q0y, 0)
	rem = rsh256(rem, s)
	return quo, rem
	//return q1*two32 + q0, (un21*two32 + un0 - q0*y) >> s
}

func (a LElem256) Mul(b interface{}) RElem {
	Hi, Lo := mul256(Uint256(a), Uint256(b.(LElem256)))
	if lessThan256(Hi, Uint256(LElem256Mod)) {
		_, r := div256(Hi, Lo, Uint256(LElem256Mod))
		return LElem256(r)
	}
	_, r := div256(Uint256{Uint128{0, 0}, Uint128{0, 0}}, Hi, Uint256(LElem256Mod))
	_, r = div256(r, Lo, Uint256(LElem256Mod))
	return Lelem256(r)
}
func uint64To256(x uint64) Uint256 {
	return Uint256{Uint128{0, 0}, Uint128{0, x}}
}

func (a LElem256) Add(b interface{}) RElem {
	out, carry := add256(Uint256(a), Uint256(b.(LElem256)), 0)
	_, rem := div256(uint64To256(carry), out, Uint256(LElem256Mod))
	return LElem256(rem)
}

func (a LElem256) Sub(b interface{}) RElem {
	return a.Add(b.(LElem256).Neg())
}

func (a LElem256) Neg() RElem {
	if equal128(a.Hi, Uint128{0, 0}) && equal128(a.Lo, Uint128{0, 0}) {
		return a
	}
	out, _ := sub256(Uint256(LElem256Mod), Uint256(a), 0)
	return LElem256(out)
}

func (a LElem256) Copy() RElem {
	return a
}
func (a LElem256) Uint64() uint64 {
	return a.Lo.Lo
}
func (a LElem256) One() RElem {
	return LElem256(Uint256{Uint128{0, 0}, Uint128{0, 1}})
}
