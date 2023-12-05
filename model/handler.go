package model

import (
	"fmt"
	"os"

	"github.com/wonderstone/GEP-MOD-2ND/functions"
	"github.com/wonderstone/GEP-MOD-2ND/gene"
	"github.com/wonderstone/GEP-MOD-2ND/genome"
	"github.com/wonderstone/GEP-MOD-2ND/genomeset"
	"gopkg.in/yaml.v3"
)

type InputValues []float64
type OutputValues []float64

type EvaluateFunc func(InputValues) OutputValues
type FitnessFunc func(OutputValues) float64
// ModelType 模型类型
type ModelType string

const (
	ModelTypeGenome    ModelType = "Genome"
	ModelTypeGenomeSet ModelType = "GenomeSet"
)

// GetEvaluateFunc 用于单个基因组的评估
func GetEvaluateFunc(g *genome.Genome) EvaluateFunc {
	return func(input InputValues) OutputValues {
		return OutputValues{g.EvalMath(input)}
	}
}

// GetEvaluateFunc2 用于多个基因组的评估
func GetEvaluateFunc2(g *genomeset.GenomeSet) EvaluateFunc {
	return func(input InputValues) OutputValues {
		o := OutputValues{}
		for _, g := range g.Genomes {
			o = append(o, g.EvalMath(input))
		}
		return o
	}
}

type Handler interface {
	Init(option ...WithOption) error

	Evolve() *Record  // 用于多次评估
	RunOnce() *Record // 用于单次运行
}

type TrainInfo struct {
	Genome    *genome.Genome       // 当前基因组
	GenomeSet *genomeset.GenomeSet // 当前基因组集合
}

type Record struct {
	Mode  string     `yaml:"mode"`
	Score float64    `yaml:"score"`
	KES   [][]string `yaml:"kes"`
}

type GepModel struct {
	Iteration              int               `yaml:"iteration"`                // 迭代次数
	Glis                   int               `yaml:"glis"`                     // 插入长度
	Glris                  int               `yaml:"glris"`                    // 重复长度
	NumGenomeSet           int               `yaml:"num-genomeset"`            // 基因组集合数量
	NumGenomes             int               `yaml:"num-genome"`               // 基因组数量
	HeadSize               int               `yaml:"head-size"`                // 头部长度
	NumGenomesPerGenomeSet int               `yaml:"num-genome-per-genomeset"` // 每个基因组集合中基因组数量
	NumGenesPerGenome      int               `yaml:"num-gene-per-genome"`      // 每个基因组中基因数量
	NumConstants           int               `yaml:"num-constants"`            // 常量数量
	PMutate                float64           `yaml:"pmutate"`                  // 变异概率
	Pis                    float64           `yaml:"pis"`                      // 插入概率
	Pris                   float64           `yaml:"pris"`                     // 重复概率
	PGene                  float64           `yaml:"pgene"`                    // 基因突变概率
	P1p                    float64           `yaml:"p1p"`                      // 一点交叉概率
	P2p                    float64           `yaml:"p2p"`                      // 两点交叉概率
	Pr                     float64           `yaml:"pr"`                       // 交换概率
	FunctionWeights        []gene.FuncWeight `yaml:"input-function"`           // [func][weight] 注入的函数权重
	LinkFunc               string            `yaml:"link-func"`                // 连接函数
	Mode                   string            `yaml:"mode"`                     // 模式
}

// *  from here, Use GEP to solve the specific problem
type Op struct {
	Conf                   *GepModel
	Perf                   PerformanceFunc
	Perf2                  PerformanceFunc2
	GenomeF                GenomeFunc
	GenomeSetF             GenomeSetFunc
	Record                 Record
	NumTerminal            int
	FuncType               functions.FuncType
	Indicator2FormulaIndex map[string]int // 指标名称到公式索引的映射
}

type WithOption func(op *Op)

// PerformanceFunc 优化函数，返回最优基因组，以及是否达到预期
type PerformanceFunc func(int, []*genome.Genome) (g *genome.Genome, accomplished bool)
type PerformanceFunc2 func(int, []*genomeset.GenomeSet) (g *genomeset.GenomeSet, accomplished bool)

type GenomeFunc func(*genome.Genome)
type GenomeSetFunc func(*genomeset.GenomeSet)

// * WithOption mode
func WithOp(opt *Op) WithOption {
	return func(op *Op) {
		*op = *opt
	}
}

func WithModelConfig(model GepModel) WithOption {
	return func(op *Op) {
		op.Conf = &model
	}
}
func WithPerformance(perf PerformanceFunc) WithOption {
	return func(op *Op) {
		op.Perf = perf
	}
}
func WithPerformanceSet(perf PerformanceFunc2) WithOption {
	return func(op *Op) {
		op.Perf2 = perf
	}
}
func WithGenome(f GenomeFunc) WithOption {
	return func(op *Op) {
		op.GenomeF = f
	}
}

func WithGenomeSet(f GenomeSetFunc) WithOption {
	return func(op *Op) {
		op.GenomeSetF = f
	}
}

func WithIndicator2FormulaIndex(m map[string]int) WithOption {
	return func(op *Op) {
		op.Indicator2FormulaIndex = m
	}
}

func WithKarvaExpressionFile(f string) WithOption {
	return func(op *Op) {
		rec := Record{
			KES: [][]string{},
		}
		content, err := os.ReadFile(f)
		if err != nil {
			fmt.Printf("读取karva表达式文件失败: %v", err)
		}

		err = yaml.Unmarshal(content, &rec)
		if err != nil {
			fmt.Printf("解析karva表达式文件失败: %v", err)
		}
		op.Record = rec
	}
}

func WithNumTerminal(num int) WithOption {
	return func(op *Op) {
		op.NumTerminal = num
	}
}
func WithFuncType(funcType functions.FuncType) WithOption {
	return func(op *Op) {
		op.FuncType = funcType
	}
}
func NewOp(option ...WithOption) *Op {
	op := &Op{
		FuncType: functions.Float64,
	}
	for _, opt := range option {
		opt(op)
	}
	if op.Conf == nil {
		fmt.Print("未指定模型配置")
	}
	if op.Perf == nil && op.Perf2 == nil && op.GenomeF == nil && op.GenomeSetF == nil {
		fmt.Print("未指定筛选函数 或者 评估函数")
	}
	if op.NumTerminal == 0 {
		fmt.Print("未指定模型参数数量，此参数由指标数量决定")
	}
	return op
}

func NewHandler(mode ModelType, option ...WithOption) Handler {
	switch mode {
	case ModelTypeGenome:
		handler := &GenomeHandler{}
		if err := handler.Init(option...); err != nil {
			fmt.Printf("初始化Genome模型失败: %v", err)
		}

		return handler
	case ModelTypeGenomeSet:
		handler := &GenomeSetHandler{}
		if err := handler.Init(option...); err != nil {
			fmt.Printf("初始化GenomeSet模型失败: %v", err)
		}
		return handler
	default:
		fmt.Printf("未知模型类型: %s", mode)
	}

	return nil
}
