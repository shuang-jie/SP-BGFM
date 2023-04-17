# SP-BGFM
SP-BGFM (Sparse Bayesian Group Factor Model for Feature Interactions in Multiple Count Tables Data) is a Bayesian semi-parametric model for multiple domians of next-generation sequencing microbiome abundance data. It models cross-domain interaction between microbial features directly using observed count tables, and provides a flexible DP mixture structure for interpretattion.

For more information, read the paper: Sparse Bayesian Group Factor Model for Feature Interactions in Multiple Count Tables Data (to be submitted). An presentation will be on JSM 2023, Toronto. 

Contact: Shuangjie Zhang (szhan209 AT ucsc DOT edu)


## Installation

SP-BGFM requires Rcpp and RcppArmadillo. Download and install R from [https://www.r-project.org/](https://www.r-project.org/). Once installed, open R from the terminal with R and run the following command:

```
install.packages(c("Rcpp", "RcppArmadillo"))
```

Please make sure the C++ compiler is correctly installed. 

All the parameters update are in update_github.Rcpp file. To load the update functions, use the command:

```
sourceCpp("~update_github.cpp")
```



## Simulation 

