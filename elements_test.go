package mpc_core

import (
	"fmt"
	"math/big"
	"math/rand"
	"runtime"
	"sync"
	"testing"
	"time"
)

func TestAddSub256(t *testing.T) {
	fmt.Println("start")
	source := rand.NewSource(0)
	r := rand.New(source)
	for i := 0; i < 10000; i++ {

		aint := new(big.Int).Rand(r, LElem256ModBig)
		bint := new(big.Int).Rand(r, LElem256ModBig)

		a := new(LElem256Big).FromBigInt(aint)
		b := new(LElem256Big).FromBigInt(bint)

		atest := new(LElem256).FromBigInt(aint)
		btest := new(LElem256).FromBigInt(bint)

		c := a.Add(b).(LElem256Big)
		ctest := atest.Add(btest).(LElem256)
		d := a.Sub(b).(LElem256Big)
		dtest := atest.Sub(btest).(LElem256)

		dif := new(big.Int).Sub(c.ToBigInt(), ctest.ToBigInt())
		dif2 := new(big.Int).Sub(d.ToBigInt(), dtest.ToBigInt())
		if !dif.IsUint64() {
			fmt.Println(c.ToBigInt(), "answer, add+div")
			fmt.Println(ctest.ToBigInt(), "test, add+div")

			out, carry := Add256(Uint256(atest.(LElem256)), Uint256(btest.(LElem256)), 0)

			fmt.Println(new(big.Int).Add(a.(LElem256Big).ToBigInt(), b.(LElem256Big).ToBigInt()), "answer, add")
			fmt.Println(carry, LElem256(out).ToBigInt(), "test, add")

			_, rem := Div256(uint64To256(carry), out, Uint256(LElem256Mod))
			fmt.Println(LElem256(rem).ToBigInt(), "test, div after add")

			t.Fatalf(fmt.Sprint("can't be a uint64:", dif))
		} else {
			if dif.Uint64() != 0 {
				fmt.Println(c.ToBigInt(), "answer, add+div")
				fmt.Println(ctest.ToBigInt(), "test, add+div")

				out, carry := Add256(Uint256(atest.(LElem256)), Uint256(btest.(LElem256)), 0)

				fmt.Println(new(big.Int).Add(a.(LElem256Big).ToBigInt(), b.(LElem256Big).ToBigInt()), "answer, add")
				fmt.Println(carry, LElem256(out).ToBigInt(), "test, add")

				_, rem := Div256(uint64To256(carry), out, Uint256(LElem256Mod))
				fmt.Println(LElem256(rem).ToBigInt(), "test, div after add")
				rem, _ = Sub256(out, Uint256(LElem256Mod), 0)
				fmt.Println(LElem256(rem).ToBigInt(), "test, sub after add")
				t.Fatalf(fmt.Sprint("not equal to zero:", dif))
			}
		}

		if !dif2.IsUint64() {
			fmt.Println(c.ToBigInt(), "answer, sub+div")
			fmt.Println(ctest.ToBigInt(), "test, sub+div")

			out, carry := Add256(Uint256(atest.(LElem256)), Uint256(btest.(LElem256)), 0)

			fmt.Println(new(big.Int).Sub(a.(LElem256Big).ToBigInt(), b.(LElem256Big).ToBigInt()), "answer, sub")
			fmt.Println(carry, LElem256(out).ToBigInt(), "test, sub")

			_, rem := Div256(uint64To256(carry), out, Uint256(LElem256Mod))
			fmt.Println(LElem256(rem).ToBigInt(), "test, div after sub")

			t.Fatalf(fmt.Sprint("can't be a uint64:", dif))
		} else {
			if dif.Uint64() != 0 {
				t.Fatalf(fmt.Sprint("not equal to zero:", dif))
			}
		}
	}

}

func TestEqualMod(t *testing.T) {
	a := LElem256(LElem256ModHalf).ToBigInt()
	b := LElem256Big(LElem256BigModHalf).ToBigInt()
	fmt.Println(a, b)
	dif := new(big.Int).Sub(a, b)
	fmt.Println(dif)
	if dif.Uint64() != 0 {
		t.Fatalf("not equal to zero")
	}
}

func TestFromToBigInt(t *testing.T) {
	source := rand.NewSource(0)
	r := rand.New(source)
	for i := 0; i < 10000; i++ {
		aint := new(big.Int).Rand(r, LElem256ModBig)
		a := new(LElem256Big).FromBigInt(aint).(LElem256Big)

		atest := new(LElem256).FromBigInt(aint).(LElem256)
		dif := new(big.Int).Sub(a.ToBigInt(), atest.ToBigInt())
		if !dif.IsUint64() {
			fmt.Println(a.ToBigInt(), atest.ToBigInt())
			t.Fatalf("can't be a uint64")
		} else {
			if dif.Uint64() != 0 {
				t.Fatalf("not equal to zero")
			}
		}
	}
}
func TestFromToBytes(t *testing.T) {
	source := rand.NewSource(0)
	r := rand.New(source)
	for i := 0; i < 10000; i++ {
		aint := new(big.Int).Rand(r, LElem256ModBig)
		a := new(LElem256Big).FromBigInt(aint).(LElem256Big)
		atest := new(LElem256).FromBigInt(aint).(LElem256)
		buf := make([]byte, 32)
		buftest := make([]byte, 32)
		a.ToBytes(buf)
		atest.ToBytes(buftest)

		c := a.FromBytes(buftest).(LElem256Big)
		ctest := atest.FromBytes(buf).(LElem256)

		dif := new(big.Int).Sub(c.ToBigInt(), ctest.ToBigInt())
		if !dif.IsUint64() {
			fmt.Println(a.ToBigInt(), atest.ToBigInt())
			t.Fatalf("can't be a uint64")
		} else {
			if dif.Uint64() != 0 {
				t.Fatalf("not equal to zero")
			}
		}
	}
}

// func TestRandBits(t *testing.T) {
// 	source := rand.NewSource(5)
// 	r := rand.New(source)
// 	counts := make(map[int]int)
// 	for i := 0; i < 100000; i++ {
// 		aint := new(big.Int).Rand(r, LElem256ModBig)
// 		n := rand.Intn(255)
// 		a := new(LElem256Big).FromBigInt(aint)
// 		atest := new(LElem256).FromBigInt(aint)
// 		seed := make([]byte, chacha.KeySize)
// 		prg := frand.NewCustom(seed, 1024, 20)

// 		c := a.RandBits(prg, n).(LElem256Big)
// 		ctest := atest.RandBits(prg, n).(LElem256)
// 		res := c.ToBigInt().BitLen()
// 		res2 := ctest.ToBigInt().BitLen()
// 		_ = res2
// 		dif := res - n
// 		if val, ok := counts[dif]; ok {
// 			counts[dif] = val + 1
// 		} else {
// 			counts[dif] = 1
// 		}
// 	}
// 	fmt.Println(counts)
// 	t.Fatalf("can't be right")
// }

func TestFromToBigFloat(t *testing.T) {
	source := rand.NewSource(0)
	r := rand.New(source)
	for i := 0; i < 10000; i++ {
		aint := new(big.Int).Rand(r, LElem256ModBig)
		a := new(LElem256Big).FromBigInt(aint).(LElem256Big)

		atest := new(LElem256).FromBigInt(aint).(LElem256)
		dif := new(big.Float).Sub(a.ToBigFloat(50), atest.ToBigFloat(50))
		ans, _ := dif.Float64()
		if ans != 0.0 {
			fmt.Println(a.ToBigInt(), atest.ToBigInt())
			t.Fatalf("can't be a uint64")
		}
	}
}
func TestMul256(t *testing.T) {
	fmt.Println()
	fmt.Println("start")
	source := rand.NewSource(0)
	r := rand.New(source)
	for i := 0; i < 50000; i++ {

		aint := new(big.Int).Rand(r, LElem256ModBig)
		bint := new(big.Int).Rand(r, LElem256ModBig)
		// aint, _ := new(big.Int).SetString("106850440303375662358867864312180327228148047350501935339945988082313438748850", 10)
		// bint, _ := new(big.Int).SetString("57039647654281895328592137303582190442763477193021006335380483862547864133732", 10)

		a := new(LElem256Big).FromBigInt(aint)
		b := new(LElem256Big).FromBigInt(bint)

		atest := new(LElem256).FromBigInt(aint)
		btest := new(LElem256).FromBigInt(bint)

		c := a.Mul(b).(LElem256Big)
		ctest := atest.Mul(btest).(LElem256)

		dif := new(big.Int).Sub(c.ToBigInt(), ctest.ToBigInt())
		if !dif.IsUint64() {
			fmt.Println(aint, "a")
			fmt.Println(bint, "b")
			fmt.Println()
			fmt.Println(c.ToBigInt(), "answer")
			fmt.Println(ctest.ToBigInt(), "test")

			// out, carry := Add256(Uint256(atest.(LElem256)), Uint256(btest.(LElem256)), 0)

			// fmt.Println(new(big.Int).Add(a.(LElem256).ToBigInt(), b.(LElem256).ToBigInt()), "answer, add")
			// fmt.Println(carry, LElem256(out).ToBigInt(), "test, add")

			// _, rem := Div256(uint64To256(carry), out, Uint256(LElem256Mod))
			// fmt.Println(LElem256(rem).ToBigInt(), "test, div after add")

			t.Fatalf(fmt.Sprint("can't be a uint64:", dif))
		} else {
			if dif.Uint64() != 0 {
				t.Fatalf(fmt.Sprint("not equal to zero:", dif))
			}
		}
	}
}
func TestNegInv(t *testing.T) {
	source := rand.NewSource(0)
	r := rand.New(source)
	for i := 0; i < 50000; i++ {
		aint := new(big.Int).Rand(r, LElem256ModBig)
		a := new(LElem256Big).FromBigInt(aint)

		atest := new(LElem256).FromBigInt(aint)

		c := a.Neg().(LElem256Big)
		ctest := atest.Neg().(LElem256)

		dif := new(big.Int).Sub(c.ToBigInt(), ctest.ToBigInt())
		if !dif.IsUint64() {
			t.Fatalf(fmt.Sprint("can't be a uint64:", dif))
		} else {
			if dif.Uint64() != 0 {
				t.Fatalf(fmt.Sprint("not equal to zero:", dif))
			}
		}
		c = a.Inv().(LElem256Big)
		ctest = atest.Inv().(LElem256)

		dif = new(big.Int).Sub(c.ToBigInt(), ctest.ToBigInt())
		if !dif.IsUint64() {
			t.Fatalf(fmt.Sprint("can't be a uint64:", dif))
		} else {
			if dif.Uint64() != 0 {
				t.Fatalf(fmt.Sprint("not equal to zero:", dif))
			}
		}
	}
}

func TestGetBit256(t *testing.T) {
	source := rand.NewSource(0)
	r := rand.New(source)
	for i := 0; i < 50000; i++ {
		aint := new(big.Int).Rand(r, LElem256ModBig)
		n := rand.Intn(255)
		a := new(LElem256Big).FromBigInt(aint)

		atest := new(LElem256).FromBigInt(aint)

		c := a.GetBit(n)
		ctest := atest.GetBit(n)

		dif := c - ctest
		if dif != 0 {
			fmt.Println(c, ctest)
			t.Fatalf("can't be a uint64")
		}
	}
}

func TestTrunc256(t *testing.T) {
	source := rand.NewSource(0)
	r := rand.New(source)
	for i := 0; i < 50000; i++ {
		aint := new(big.Int).Rand(r, LElem256ModBig)
		n := rand.Intn(255)
		a := new(LElem256Big).FromBigInt(aint)

		atest := new(LElem256).FromBigInt(aint)

		c := a.Trunc(n).(LElem256Big)
		ctest := atest.Trunc(n).(LElem256)

		dif := new(big.Int).Sub(c.ToBigInt(), ctest.ToBigInt())
		if !dif.IsUint64() {
			fmt.Println(c.ToBigInt(), ctest.ToBigInt())
			t.Fatalf("can't be a uint64")
		} else {
			if dif.Uint64() != 0 {
				t.Fatalf("not equal to zero")
			}
		}
	}
}

func TestFromToInt(t *testing.T) {
	source := rand.NewSource(0)
	r := rand.New(source)
	for i := 0; i < 10000; i++ {
		aint := r.Int()
		a := new(LElem256Big).FromInt(aint).(LElem256Big)
		atest := new(LElem256).FromInt(aint).(LElem256)
		dif := new(big.Int).Sub(a.ToBigInt(), atest.ToBigInt())
		ans := dif.Uint64()
		if ans != 0.0 {
			fmt.Println(a.ToBigInt(), atest.ToBigInt())
			t.Fatalf("can't be a uint64")
		}
	}
}
func TestFromToUInt64(t *testing.T) {
	source := rand.NewSource(0)
	r := rand.New(source)
	for i := 0; i < 10000; i++ {
		aint := r.Uint64()
		a := new(LElem256Big).FromUint64(aint).(LElem256Big)
		atest := new(LElem256).FromUint64(aint).(LElem256)
		dif := new(big.Int).Sub(a.ToBigInt(), atest.ToBigInt())
		ans := dif.Uint64()
		if ans != 0.0 {
			fmt.Println(a.ToBigInt(), atest.ToBigInt())
			t.Fatalf("can't be a uint64")
		}
	}
}
func TestFromToFloat64(t *testing.T) {
	// source := rand.NewSource(0)
	// r := rand.New(source)
	for i := 0; i < 10000; i++ {
		aint := 1.203410921
		a := new(LElem256Big).FromFloat64(aint, 50).(LElem256Big)
		atest := new(LElem256).FromFloat64(aint, 50).(LElem256)
		dif := new(big.Int).Sub(a.ToBigInt(), atest.ToBigInt())
		ans := dif.Uint64()
		if ans != 0.0 {
			fmt.Println(a.ToBigInt(), atest.ToBigInt())
			t.Fatalf("can't be a uint64")
		}
	}
}

func TestChangeInLoop(t *testing.T) {
	runtime.GOMAXPROCS(15)
	nblocks := []int{1, 5, 15}
	for _, nb := range nblocks {
		startBlock := time.Now()
		var wg sync.WaitGroup
		for i := 0; i < nb; i++ {
			wg.Add(1)
			// go func(i int, mpcObj *mpc.MPC, mod *big.Int, m mpc_core.RVec) {
			go func(i int) {
				defer wg.Done()
				a := new(big.Int).SetUint64(132)
				b := new(big.Int).SetUint64(1123232)
				c := new(big.Int).SetUint64(1323232)
				start := time.Now()
				for j := 0; j < 5*50; j++ {
					for k := 0; k < 1000; k++ {
						// add128(a, b)
						c.Add(a, b)
						b.Mod(c, a)
						// d := new(big.Int).Add(c, a)
						// _ = d
					}
					// m.add(m)

				}
				fmt.Printf("%d-%d:", nb, i)
				fmt.Println(time.Since(start))
			}(i)
			// }(i, prot.GetMpc()[i], new(big.Int).Set(mpc_core.LElem256Mod.Val), mpc_core.InitRVec(rtype.One(), 1000))
		}
		wg.Wait()
		fmt.Printf("[Total] %d blocks:", nb)
		fmt.Println(time.Since(startBlock))
	}
}
