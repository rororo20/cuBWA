#  cuBWA

## Method

### SMEMS
extend过程中需要SA interval，也就是需要继续k，l两个值。foward和backwards需要算4个碱基。这样并发度才8个，这样占不满wrap，需要按如下索引需要调整一下  
假如每隔8bp采样，bwt数组调整为 [0,4,1,5,2,6,3,7]  
thread 0 被分配成 [0,4]  
thread 1 被分配成 [1,5]  
thread 2 被分配成 [2,6]  
thread 3 被分配成 [3,7]  
最终需要wrap shuffle来统计个数，以及需要计算mask  

*Mask策略如下:*
假如k,l的采样偏移量为0，thread0-thread3都有0x0;  
假如k,l的采样偏移量为1，thread0为0x01， 其余0x0;  
假如k,l的采样偏移量为2，thread0-thread1都有0x01，其余0;  
假如k,l的采样偏移量为3，thread0-thread2都有0x01，其余0;  
假如k,l的采样偏移量为4，thread0都有0x02，其余0x01;  
...
假如k,l的采样偏移量为7，thread0-thread3都有0x02;  
...
这样4*8=32线程，这样wrap利用率会比较高 
sample文件需要测试一下随机读写的效率，估算30X WGS数据要Extension多少次，来推导该模块的上限； 

### SAI2coordinate
该部分可以在CPU实现（倾向）,也可以GPU实现，同时可以计算多个suffix Array，如果不采样，则可以一起并发，算完后把直接可以访问到SA;  

### Chaining
该部分可以在CPU实现，如果不行，则可以抄论文里面的；  