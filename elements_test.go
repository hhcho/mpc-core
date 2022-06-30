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
	source := rand.NewSource(0)
	r := rand.New(source)
	for i := 0; i < 10; i++ {
		aint := new(big.Int).Rand(r, LElem256testModBig)
		bint := new(big.Int).Rand(r, LElem256testModBig)

		a := new(LElem256).FromBigInt(aint)
		b := new(LElem256).FromBigInt(bint)

		atest := new(LElem256test).FromBigInt(aint)
		btest := new(LElem256test).FromBigInt(bint)

		c := a.Add(b).(LElem256)
		ctest := atest.Add(btest).(LElem256test)

		dif := new(big.Int).Sub(c.ToBigInt(), ctest.ToBigInt())
		if !dif.IsUint64() {
			fmt.Println(LElem256testModBig)
			fmt.Println(LElem256Mod)
			fmt.Println(a.(LElem256).ToBigInt())
			fmt.Println(atest.(LElem256test).ToBigInt())
			fmt.Println(b.(LElem256).ToBigInt())
			fmt.Println(btest.(LElem256test).ToBigInt())
			fmt.Println(c.ToBigInt(), "answer, add+div")
			fmt.Println(ctest.ToBigInt(), "test, add+div")

			out, carry := add256(Uint256(atest.(LElem256test)), Uint256(btest.(LElem256test)), 0)

			fmt.Println(new(big.Int).Add(a.(LElem256).ToBigInt(), b.(LElem256).ToBigInt()), "answer, add")
			fmt.Println(carry, LElem256test(out).ToBigInt(), "test, add")

			_, rem := div256(uint64To256(carry), out, Uint256(LElem256testMod))
			fmt.Println(LElem256test(rem).ToBigInt(), "test, div after add")

			t.Fatalf(fmt.Sprint("can't be a uint64:", dif))
		} else {
			if dif.Uint64() != 0 {
				t.Fatalf(fmt.Sprint("not equal to zero:", dif))
			}
		}
	}

}

func TestFromToBigInt(t *testing.T) {
	source := rand.NewSource(0)
	r := rand.New(source)
	for i := 0; i < 10; i++ {
		aint := new(big.Int).Rand(r, LElem256testModBig)
		a := new(LElem256).FromBigInt(aint).(LElem256)

		atest := new(LElem256test).FromBigInt(aint).(LElem256test)
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
func TestMul256(t *testing.T) {

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
