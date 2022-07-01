package mpc_core

import (
	"encoding/binary"
	"math/big"

	"github.com/hhcho/frand"
)

type Uint256 struct {
	Hi, Lo Uint128
}

type LElem256 Uint256

const LElem256Bytes uint32 = 32

func initLElem256RandBnd() LElem256 {
	mOne := Uint256{Uint128{^uint64(0), ^uint64(0)}, Uint128{^uint64(0), ^uint64(0)}}
	_, rem := div256(uint64To256(0), mOne, Uint256(LElem256Mod))
	out, _ := sub256(mOne, rem, 0)
	return LElem256(out)
}

var LElem256RandBnd LElem256 = initLElem256RandBnd()
var LElem256Zero LElem256 = LElem256{Uint128{0, 0}, Uint128{0, 0}}

var LElem256Mod = Uint256{Uint128{18446744073709551615, 18446744073709551615}, Uint128{18446744073709551615, 18446744073709551427}}    // 2^256 - 189
var LElem256ModHalf = Uint256{Uint128{9223372036854775807, 18446744073709551615}, Uint128{18446744073709551615, 18446744073709551522}} // 2^255 - 94
var LElem256ModBitLen int = 256

const LElem256UniqueID uint8 = 5

var LElem256ModBig *big.Int = new(big.Int).Sub(big.NewInt(0).Lsh(big.NewInt(1), 256), big.NewInt(189))
var LElem256ModHalfBig *big.Int = new(big.Int).Sub(big.NewInt(0).Lsh(big.NewInt(1), 255), big.NewInt(94))

func (a LElem256) ToBigInt() *big.Int {
	hii := new(big.Int).Lsh(new(big.Int).SetUint64(a.Hi.Hi), 256-64)
	hilo := new(big.Int).Lsh(new(big.Int).SetUint64(a.Hi.Lo), 256-128)
	lohi := new(big.Int).Lsh(new(big.Int).SetUint64(a.Lo.Hi), 64)
	lo := new(big.Int).SetUint64(a.Lo.Lo)

	sumHi := new(big.Int).Add(hii, hilo)
	sumLo := new(big.Int).Add(lohi, lo)
	sum := new(big.Int).Add(sumHi, sumLo)
	r := new(big.Int)
	new(big.Int).QuoRem(sum, LElem256ModBig, r)
	return r
}

//sus
func (a LElem256) FromBigInt(n *big.Int) RElem {
	hii := big.NewInt(0).Rsh(n, 256-64)
	hilo := big.NewInt(0).Rsh(big.NewInt(0).Mod(n, big.NewInt(0).Lsh(big.NewInt(1), 256-64)), 256-128)
	lohi := big.NewInt(0).Rsh(big.NewInt(0).Mod(n, big.NewInt(0).Lsh(big.NewInt(1), 256-128)), 64)
	lo := big.NewInt(0).Mod(n, big.NewInt(0).Lsh(big.NewInt(1), 64))
	return LElem256(Uint256{Uint128{hii.Uint64(), hilo.Uint64()}, Uint128{lohi.Uint64(), lo.Uint64()}})
}

func (a LElem256) ToBigFloat(fracBits int) *big.Float {
	out := new(big.Float).SetInt(a.ToBigInt())
	out.Mul(out, new(big.Float).SetMantExp(big.NewFloat(1), -fracBits))
	return out
}

func (a LElem256) ToSignedBigInt() *big.Int {
	tmp := a.ToBigInt()
	half := new(big.Int).Rsh(LElem256ModBig, 1)
	if tmp.Cmp(half) >= 0 {
		tmp.Sub(tmp, LElem256ModBig)
	}
	return tmp
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

func mul256(a, b Uint256) (Hi, Lo Uint256) {
	h00, r0 := mul128(a.Lo, b.Lo) //r0
	h01, l01 := mul128(a.Lo, b.Hi)
	h10, l10 := mul128(a.Hi, b.Lo)
	h11, l11 := mul128(a.Hi, b.Hi)

	t, carry1 := add128(h00, l01, 0)
	r1, carry2 := add128(t, l10, 0) //r1
	t, carry1 = add128(l11, h01, carry1)
	r2, carry2 := add128(t, h10, carry2)
	r3, _ := add128(h11, Uint128{0, carry1}, carry2)

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
		s.Lo = or128(rsh128(u.Lo, n), lsh128(u.Hi, 128-n))
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
	un32 := or256(lsh256(Hi, s), rsh256(Lo, 256-s)) //correct
	un10 := lsh256(Lo, s)
	un1 := Uint256{Uint128{0, 0}, un10.Hi}
	un0 := Uint256{Uint128{0, 0}, un10.Lo}
	q1hi, r := div128(Uint128{0, 0}, un32.Hi, yn1.Lo)
	q1lo, _ := div128(r, un32.Lo, yn1.Lo)
	q1 := Uint256{q1hi, q1lo}     //should be and is 0
	_, t := mul256(q1, yn1)       //should be and is 0
	rhat, _ := sub256(un32, t, 0) //rhat := un32 - q1*yn1 //should be and is 0
	_, t = mul256(q1, yn0)        //should be and is 0
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

	_, r := div256(uint64To256(0), Hi, Uint256(LElem256Mod))
	_, r = div256(r, Lo, Uint256(LElem256Mod))
	return LElem256(r)
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

func (a LElem256) Inv() RElem {
	bInv := big.NewInt(0).ModInverse(a.ToBigInt(), LElem256ModBig)
	if bInv == nil {
		panic("ModInverse does not exist")
	}
	return a.FromBigInt(bInv)
}

func (a LElem256) AssertTypeFor(n RElem) RElem {
	return n.(LElem256)
}

func (a LElem256) FromInt(n int) RElem {
	if n >= 0 {
		return LElem256(uint64To256(0)).Add(uint64To256(uint64(n)))
	}
	return LElem256(uint64To256(0)).Sub(uint64To256(uint64(-n)))
}

func (a LElem256) FromUint64(n uint64) RElem {
	return LElem256(uint64To256(n))
}

func (a LElem256) FromFloat64(n float64, fracBits int) RElem {
	if n < 0 {
		return LElem256{Hi: Uint128{0, 0}, Lo: Uint128(LElem128Zero.FromFloat64(-n, fracBits).(LElem128))}.Neg()
	} else {
		return LElem256{Hi: Uint128{0, 0}, Lo: Uint128(LElem128Zero.FromFloat64(n, fracBits).(LElem128))}
	}
}

func (a LElem256) Float64(fracBits int) float64 {
	var sgn int
	var b LElem256
	if lessThan256(Uint256(a), LElem256ModHalf) {
		sgn, b = 1, a
	} else {
		sgn, b = -1, a.Neg().(LElem256)
	}

	shift := float64(uint64(1) << fracBits)
	const shift64 = 1 << 64
	const shift128 = 1 << 128
	const shift192 = 1 << 192

	val := float64(b.Lo.Lo) / shift
	val += float64(b.Lo.Hi) / shift * shift64
	val += float64(b.Hi.Lo) / shift * shift128
	val += float64(b.Hi.Hi) / shift * shift192

	return float64(sgn) * val
}

func (a LElem256) FromBytes(buf []byte) RElem {
	return LElem256(Uint256{frombytes128(buf[:LElem128Bytes]), frombytes128(buf[LElem128Bytes:])})
}

func (a LElem256) ToBytes(buf []byte) {
	binary.LittleEndian.PutUint64(buf, a.Hi.Hi)
	binary.LittleEndian.PutUint64(buf[8:LElem128Bytes], a.Hi.Lo)
	binary.LittleEndian.PutUint64(buf[LElem128Bytes:LElem128Bytes+8], a.Lo.Hi)
	binary.LittleEndian.PutUint64(buf[LElem128Bytes+8:], a.Lo.Lo)
	return
}
func (a LElem256) Zero() RElem {
	return LElem256(Uint256{Uint128{0, 0}, Uint128{0, 0}})
}

func (a LElem256) Modulus() *big.Int {
	r := new(big.Int).Set(LElem128ModBig)
	return r
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
	if posFromLSB < 64 {
		return boolToUint((a.Lo.Lo & (uint64(1) << posFromLSB)) > 0)
	} else if posFromLSB < 128 {
		return boolToUint((a.Lo.Hi & (uint64(1) << (posFromLSB - 64))) > 0)
	} else if posFromLSB < 192 {
		return boolToUint((a.Hi.Lo & (uint64(1) << (posFromLSB - 128))) > 0)
	}
	return boolToUint((a.Hi.Hi & (uint64(1) << (posFromLSB - 192))) > 0)
}

func (a LElem256) Trunc(nBits int) RElem {
	if nBits < 0 || nBits > a.ModBitLength() {
		panic("Invalid number of bits")
	} else if nBits == a.ModBitLength() {
		return a
	}
	if nBits > 256-64 {
		tmp := a.Hi.Hi % (uint64(1) << (nBits - (256 - 64)))
		return LElem256{Uint128{tmp, a.Hi.Lo}, a.Lo}
	} else if nBits == 256-64 {
		return LElem256{Uint128{0, a.Hi.Lo}, a.Lo}
	} else if nBits > 256-128 {
		tmp := a.Hi.Lo % (uint64(1) << (nBits - (256 - 128)))
		return LElem256{Uint128{0, tmp}, a.Lo}
	} else if nBits == 256-128 {
		return LElem256{Uint128{0, 0}, a.Lo}
	} else if nBits > 64 {
		tmp := a.Lo.Hi % (uint64(1) << (nBits - (64)))
		return LElem256{Uint128{0, 0}, Uint128{tmp, a.Lo.Lo}}
	} else if nBits == 64 {
		return LElem256{Uint128{0, 0}, Uint128{0, a.Lo.Lo}}
	}
	tmp := a.Lo.Lo % (uint64(1) << nBits)
	return LElem256{Uint128{0, 0}, Uint128{0, tmp}}
}
func (a LElem256) TypeID() uint8 {
	return LElem128UniqueID
}

func (a LElem256) Rand(prg *frand.RNG) RElem {
	buf := make([]byte, a.NumBytes())
again:
	r := randBytes(a, buf, prg).(LElem256)
	if !lessThan256(Uint256(r), Uint256(LElem256RandBnd)) {
		goto again
	}
	_, rem := div256(uint64To256(0), Uint256(r), Uint256(LElem256Mod))
	return LElem256(rem)
}

// check for length
func (a LElem256) RandBits(prg *frand.RNG, nbits int) RElem {
	if nbits >= a.ModBitLength() {
		panic("Requested bit length is larger than modulus")
	}
	buf := make([]byte, a.NumBytes())
	out := randBytes(a, buf, prg).(LElem256)
	if nbits > 256-64 {
		out.Hi.Hi = out.Hi.Hi % (uint64(1) << (nbits - (256 - 64)))
	} else if nbits == 256-64 {
		out.Hi.Hi = 0
	} else if nbits > 256-128 {
		out.Hi.Hi = 0
		out.Hi.Lo = out.Hi.Lo % (uint64(1) << (nbits - (256 - 128)))
	} else if nbits == 256-128 {
		out.Hi.Hi = 0
		out.Hi.Lo = 0
	} else if nbits > 64 {
		out.Hi.Hi = 0
		out.Hi.Lo = 0
		out.Lo.Hi = out.Lo.Hi % (uint64(1) << (nbits - (64)))
	} else if nbits == 64 {
		out.Hi.Hi = 0
		out.Hi.Lo = 0
		out.Lo.Hi = 0
	} else {
		out.Hi.Hi = 0
		out.Hi.Lo = 0
		out.Lo.Hi = 0
		out.Lo.Lo = out.Lo.Lo % (uint64(1) << nbits)
	}
	return out
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
