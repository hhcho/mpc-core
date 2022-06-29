package mpc_core

import (
	"math/big"

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

func boolToUint(b bool) uint {
	if b {
		return 1
	}
	return 0
}

func randBytes(rtype RElem, buf []byte, prg *frand.RNG) RElem {
	prg.Read(buf)
	return rtype.FromBytes(buf)
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
