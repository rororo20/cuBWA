# cuBWA

cuBWA is a GPU-accelerated implementation toolkit of the BWA algorithm, currently under active development.

## Methodology

### SMEMS Optimization
During sequence extension phase, the algorithm requires SA interval calculations involving k and l parameters. To address the limitation of low thread utilization (only 8 concurrent threads for 4-base forward/backward processing), we propose the following optimizations:

**BWT Array Restructuring**
- Implement 8bp sampling with restructured BWT array indexing: [0,4,1,5,2,6,3,7]
- Thread allocation scheme:
  - Thread 0: [0,4]
  - Thread 1: [1,5] 
  - Thread 2: [2,6]
  - Thread 3: [3,7]

**Masking Strategy**
Developed a dynamic mask allocation protocol based on k/l offsets:

| Offset | Thread 0 | Thread 1 | Thread 2 | Thread 3 |
|--------|----------|----------|----------|----------|
| 0      | 0x0      | 0x0      | 0x0      | 0x0      |
| 1      | 0x02     | 0x0      | 0x0      | 0x0      |
| 2      | 0x02     | 0x02     | 0x0      | 0x0      |
| 3      | 0x02     | 0x02     | 0x02     | 0x0      |
| 4      | 0x02     | 0x02     | 0x02     | 0x02     |
| 5      | 0x03     | 0x02     | 0x02     | 0x02     |
| 6      | 0x03     | 0x03     | 0x02     | 0x02     |
| 7      | 0x03     | 0x03     | 0x03     | 0x02     |

**Performance Optimization**
- Achieves 32-thread/warp utilization through warp shuffle operations
- Conducting random I/O efficiency tests on sample datasets
- Estimating extension iterations for 30X WGS data to determine module throughput limits

### SAI2Coordinate Implementation
**Architecture Options**
1. *CPU Implementation (Preferred)*
   - Batch processing of multiple suffix arrays
   - Direct SA access without sampling

2. GPU Implementation (Alternative)
   - Parallel computation capability
   - Requires memory access optimization

### Chaining Module
**Implementation Strategy**
- Primary development on CPU architecture
- Potential GPU adaptation using established algorithms from recent research
- Maintaining flexibility for hybrid computing approaches

This technical documentation presents optimized GPU utilization strategies for BWA-algorithm components, with particular focus on warp-level parallelism and memory access patterns. Current development emphasizes balancing computational efficiency with implementation complexity.