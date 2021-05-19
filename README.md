# Definite Non-Ancestral Relations and Structure Learning
This is an implemetation of the structure learning algorithm using definite non-ancestral (DNA) relations as intermediate instrument. A modified version of sparsest permutation algorithm (SP) is implemneted in this package. 

## Installation
```R
  # install.packages("devtools")
  devtools::install_github("WY-Chen/DAGDNA")
  ```

## DNA Learning Methods
Learning DNA from sample.
```R
LearnDNAforward(
  dat,  #dat$X hold the data matrix. 
  k,    #learning level
  CIFUN_skel, #Conditional independence test used for skeleton learning. 
              #CIFUN(x,y,S) returns 0 (dependence) or 1 (independence)
  CIFUN_DNA   #Conditional independence test used for DNA learning. 
              #CIFUN(x,y,S) returns 0 (dependence) or 1 (independence)
  )
```
Deduce all DNA from known DAG.
```R
allDNA(
  g #DAG matrix (g_ij=1 if and only if i->j)
  )
```
Given an arbitrary set of DNA, deduce an order-constraining DNA set and layering of DAG.
```R
orderConstraining(
  D, #DNA output
  verbose=F
  )
```

## Modified Structure Learning Methods

```R
  SparsestPermutation(
    X,      #observed data
    DNA=F,  #use DNA learning or not
    k=0,    #if DNA=T, set level of DNA learning
    alpha=0.05,  #significance level in CI test (gaussCItest)
    r=5,    #number of restarts
    d=5,    #depth of SP search
    verbose=F #if true, prints the learning path with scores
    )
  ```
