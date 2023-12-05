package main

import (
	"fmt"
	"math"
	"testing"

	"github.com/wonderstone/GEP-MOD-2ND/functions"
	"github.com/wonderstone/GEP-MOD-2ND/gene"
	"github.com/wonderstone/GEP-MOD-2ND/genome"
	"github.com/wonderstone/GEP-MOD-2ND/genomeset"

	// "github.com/wonderstone/GEP-MOD-2ND/grammars"
	"github.com/wonderstone/GEP-MOD-2ND/model"
)

// func TestGEP(t *testing.T) {
// 	m, _ := NewModelConfig("/Users/alexxiong/GolandProjects/quant/strategy/S20230824-090000-000/.vqt/train/T20230824-110000-000/input/model.yaml")

// 	model.NewHandler(
// 		model.ModelTypeGenome,
// 		model.WithModelConfig(*m.Gep),
// 		model.WithKarvaExpressionFile("/Users/alexxiong/GolandProjects/quant/strategy/S20230824-090000-000/.vqt/bt/B20230824-110000-000/input/expression.yaml"),
// 		model.WithGenome(
// 			func(g *genome.Genome) {
// 				fmt.Println(
// 					g.EvalMath(
// 						[]float64{
// 							7.68, 7.27, 7.33, 7.3, 24222613, 176749046, 7.32, 7.31, 7.31, 7.31,
// 							// 15.3, 14.31, 15.01, 14.56, 30688983, 445475572, math.NaN(), math.NaN(), math.NaN(),
// 						},
// 					),
// 				)
// 				gr, _ := grammars.LoadGoMathGrammar()

// 				fmt.Println(g.WriteExpsStr(gr, []string{"Open", "Close", "High", "Low", "MA3", "MA5", "MA8", "MA10"}))
// 			},
// 		),

// 		model.WithNumTerminal(8),
// 	).RunOnce()

// }

// srTests is a random sample of inputs and outputs for the function "a^4 + a^3 + a^2 + a"
var srTests = []struct {
	in  []float64
	out float64
}{
	// {[]float64{0}, 0},
	{[]float64{2.81}, 95.2425},
	{[]float64{6}, 1554},
	{[]float64{7.043}, 2866.55},
	{[]float64{8}, 4680},
	{[]float64{10}, 11110},
	{[]float64{11.38}, 18386},
	{[]float64{12}, 22620},
	{[]float64{14}, 41370},
	{[]float64{15}, 54240},
	{[]float64{20}, 168420},
	{[]float64{100}, 101010100},
	{[]float64{-100}, 99009900},
}

func validateFunc(g *genome.Genome) float64 {
	result := 0.0
	for _, n := range srTests {
		r := g.EvalMath(n.in)
		// fmt.Printf("r=%v, n.in=%v, n.out=%v, g=%v\n", r, n.in, n.out, g)
		if math.IsInf(r, 0) {
			return 0.0
		}
		fitness := math.Abs(r - n.out)
		fitness = 1000.0 / (1.0 + fitness) // fitness is normalized and max value is 1000
		// fmt.Printf("r=%v, n.in=%v, n.out=%v, fitness=%v, g=%v\n", r, n.in, n.out, fitness, g)
		result += fitness
	}
	return result / float64(len(srTests))
}

func validateFuncGS(gs *genomeset.GenomeSet) float64 {
	result := 0.0
	for _, n := range srTests {
		fitness := 0.0

		for _, g := range gs.Genomes {
			r := g.EvalMath(n.in)
			// fmt.Printf("r=%v, n.in=%v, n.out=%v, g=%v\n", r, n.in, n.out, g)
			if math.IsInf(r, 0) {
				return 0.0
			}
			fitness += r
			// fmt.Printf("r=%v, n.in=%v, n.out=%v, fitness=%v, g=%v\n", r, n.in, n.out, fitness, g)

		}

		fitness = math.Abs(fitness - n.out)

		fitness = 1000.0 / (1.0 + fitness) // fitness is normalized and max value is 1000

		result += fitness
	}
	return result / float64(len(srTests))
}


// 测试Genome和符号+-*/
func TestGenomeSymbol(t *testing.T) {
	funcs := []gene.FuncWeight{
		{Symbol: "+", Weight: 1},
		{Symbol: "-", Weight: 1},
		{Symbol: "*", Weight: 1},
		{Symbol: "/", Weight: 1},
	}
	numIn := len(srTests[0].in)
	e := model.New(funcs, functions.Float64, 50, 4, 2, numIn, 0, "+", validateFunc)
	s := e.Evolve(1000, 50000, 0.8, 0.5, 3, 0.5, 3, 0.5, 0.01, 0.01, 0.01)
	
	fmt.Printf("// (a^4 + a^3 + a^2 + a) solution karva expression:\n// %q, score=%v\n", s, validateFunc(s))

	fmt.Println(s.Genes, validateFunc(s))

	gene1 := gene.New("*.+.*.*.d0.d0.d0.d0.d0",functions.Float64)
	gene2 := gene.New("*.+.d0./.d0.d0.d0.d0.d0",functions.Float64)
	genome1 := genome.New([]*gene.Gene{gene1, gene2}, "+")
	fmt.Println(validateFunc(genome1))
}
// 测试GenomeSet和符号+-*/
func TestGenomeSetSymbol(t *testing.T) {
	funcs := []gene.FuncWeight{
		{Symbol: "+", Weight: 1},
		{Symbol: "-", Weight: 1},
		{Symbol: "*", Weight: 1},
		{Symbol: "/", Weight: 1},
	}
	numIn := len(srTests[0].in)
	es := model.NewGS(funcs, functions.Float64, 50, 4, 2,2, numIn, 0, "+", validateFuncGS)
	ss := es.EvolveGS(1000, 50000, 0.8, 0.5, 3, 0.5, 3, 0.5, 0.01, 0.01, 0.01)


	fmt.Println(ss.Genomes, validateFuncGS(ss))

	gene1 := gene.New("*.*.d0.d0.d0.d0.d0.d0.d0",functions.Float64)
	gene2 := gene.New("*.*.d0.*.d0.d0.d0.d0.d0",functions.Float64)
	genome1 := genome.New([]*gene.Gene{gene1, gene2}, "+")
	gene3 := gene.New("-.*.d0.d0.d0.d0.d0.d0.d0",functions.Float64)
	gene4 := gene.New("+.d0.d0.d0.d0.d0.d0.d0.d0",functions.Float64)
	genome2 := genome.New([]*gene.Gene{gene3, gene4}, "+")

	gs := genomeset.New([]*genome.Genome{genome1, genome2}, "+")
	fmt.Println(validateFuncGS(gs))

}
