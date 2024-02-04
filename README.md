# SP-BGFM
SP-BGFM (Sparse Bayesian Group Factor Model for Feature Interactions in Multiple Count Tables Data) is a Bayesian semi-parametric model for multiple domains of next-generation sequencing microbiome abundance data. It models cross-domain interaction between microbial features directly using observed count tables, and provides a flexible DP mixture structure for interpretation.

Code to reproduce simulation results, figures, and real data analysis results from the paper "Sparse Bayesian Group Factor Model for Feature Interactions in Multiple Count Tables Data" by Shuangjie Zhang, Yuning Shen, Irene A. Chen and Juhee Lee.

Contact: Shuangjie Zhang (szhan209 AT ucsc DOT edu)

## Environment Setup

SP-BGFM requires R packages, especially Rcpp and RcppArmadillo for Rc++ functions: 

```
install.packages(c("Rcpp", "RcppArmadillo", "statmod", "GIGrvg", "extraDistr", "abind", "mvnfast", "mvnfast", "statmod", "extraDistr"))
```

For a comparison of the SPIEC-EASI method, please install

```
install.packages("SpiecEasi")
```

Please make sure the C++ compiler is correctly installed. For Mac user, please install Xcode from command line tools. Execute the command ```xcode-select --install``` on Terminal.

## Organization

### simulation code

Running our_type1.R with this.sim.id between 1-6 produces the results displayed in Figure 4.

Running cond_power.R with this.sim.id between 1-3 produces the results displayed in Figure 5.

Running our_type1_est.R with this.sim.id between 1-4 produces the results displayed in Figure S1.

Running power.R with this.sim.id between 1-9 produces the results displayed in Figures S2-S3.

