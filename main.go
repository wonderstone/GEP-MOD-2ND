package main

import (
	// "fmt"
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"time"

	// "github.com/jinzhu/copier"
	// "github.com/wonderstone/GEP-MOD-2ND/functions"
	"github.com/wonderstone/GEP-MOD-2ND/functions"
	"github.com/wonderstone/GEP-MOD-2ND/gene"
	"github.com/wonderstone/GEP-MOD-2ND/genome"

	// "github.com/wonderstone/GEP-MOD-2ND/grammars"
	"github.com/wonderstone/GEP-MOD-2ND/model"
)

func main() {

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
	// 1. 设置GepModel实例的参数如下
	// iterations , expectFitness , pm , pis , glis , pris , glris , pgene , p1p , p2p , pr
	// 1000, 50000, 0.8, 0.5, 3, 0.5, 3, 0.5, 0.01, 0.01, 0.01
	// fs []gene.FuncWeight, funcType functions.FuncType,
	// numGenomes, headSize, numGenesPerGenome, numTerminals, numConstants int,
	// linkFunc string, sf genome.ScoringFunc,

	FW := []gene.FuncWeight{
		{Symbol: "+", Weight: 1},
		{Symbol: "-", Weight: 1},
		{Symbol: "*", Weight: 1},
		{Symbol: "/", Weight: 1},
	}

	GM := model.GepModel{
		Iteration:              1000,
		Glis:                   3,
		Glris:                  3,
		NumGenomeSet:           200,
		NumGenomes:             200,
		HeadSize:               5,
		NumGenomesPerGenomeSet: 2,
		NumGenesPerGenome:      2,
		NumConstants:           0,
		PMutate:                0.8,
		Pis:                    0.5,
		Pris:                   0.5,
		PGene:                  0.5,
		P1p:                    0.01,
		P2p:                    0.01,
		Pr:                     0.01,
		FunctionWeights:        FW,
		LinkFunc:               "+",
		Mode:                   "Genome",
	}
	FitnessFunction := func(g *genome.Genome) float64 {
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
	// give an expectation
	ExpectFitness := 1000.0

	// 2. 设置PerformanceFunc
	PF := func(iteration int, genomes []*genome.Genome) (g *genome.Genome, accomplished bool) {
		// 1. 从genomes中选择最优的一个
		bestScore := 0.0
		bestGenome := genomes[0]
		c := make(chan *genome.Genome)
		for i := 0; i < len(genomes); i++ { // Evaluate genomes concurrently
			go genomes[i].EvaluateWithScore(FitnessFunction, c)
		}
		for i := 0; i < len(genomes); i++ { // Collect and return the highest scoring Genome
			gn := <-c
			if gn.Score > bestScore {
				bestGenome = gn
				bestScore = gn.Score
			}
		}
		// 2. 判断是否完成，完成包含达到预期与外部停止
		// 这里不对外部停止进行处理，可通过结构体属性来实现

		if bestScore >= ExpectFitness {
			return bestGenome, true
		}
		return bestGenome, false
	}
	// 3. 设置GenomeFunc
	GF := func(g *genome.Genome) {
		fmt.Println("Genome: ", g.StringSlice())
	}



	Modetype := model.ModelTypeGenome
	Opls := []model.WithOption{
		model.WithModelConfig(GM),
		model.WithPerformance(PF),
		model.WithGenome(GF),
		model.WithNumTerminal(1),
		model.WithFuncType(functions.Float64),
	}

	HL := model.NewHandler(Modetype, Opls...)
	fmt.Println(HL.Evolve())
	HL.RunOnce()

}

func init() {
	// check the compiler version, if it is less than 1.20,
	// otherwise use the old rand.Seed() function
	ver := runtime.Version()
	if ver <= "go1.20" {
		rand.Seed(time.Now().UnixNano())
	}
}
