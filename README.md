# SP-BGFM
SP-BGFM (Sparse Bayesian Group Factor Model for Feature Interactions in Multiple Count Tables Data) is a Bayesian semi-parametric model for multiple domains of next-generation sequencing microbiome abundance data. It models cross-domain interaction between microbial features directly using observed count tables, and provides a flexible DP mixture structure for interpretation.

Code to reproduce simulation results, figures, and real data analysis results from the paper "Sparse Bayesian Group Factor Model for Feature Interactions in Multiple Count Tables Data" by Shuangjie Zhang, Yuning Shen, Irene A. Chen and Juhee Lee.

Contact: Shuangjie Zhang (szhan209 AT tamu DOT edu)

## Environment Setup

SP-BGFM requires R packages, especially Rcpp and RcppArmadillo for Rc++ functions: 

```
install.packages(c("Rcpp", "RcppArmadillo", "statmod", "GIGrvg", "extraDistr", "abind", "mvnfast", "LaplacesDemon", "latex2exp", "tikzDevice"))
```



For a comparison of the SPIEC-EASI method, please install

```
library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
```

For a comparison of the metagenomeSeq method, please install

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

Then, we are able to install metagenomeSeq package and ComplexHeatmap for plotting:

```
BiocManager::install("ComplexHeatmap")
BiocManager::install("metagenomeSeq")
```

Please make sure the C++ compiler is correctly installed. For Mac user, please install Xcode from command line tools. Execute the command ```xcode-select --install``` on Terminal.


## R project reproducing instruction

An user friendly way to reproduce the results of the main text and the algorithm, one can create a R project by following the steps: Open Rstudio - File - New Project - Version Control - Git - enter the github link: https://github.com/shuang-jie/SP-BGFM.

Then, All simulation codes and real data analysis codes are in each sub-folder. All figures in the maintext can be reproduced by the figures-code.

## Organization

### simulation-code

Running Sim 1-5.R produces the results displayed in Sim 1-5 in the paper and please save the result as RData with this.sim.id.

### real-data

Filtered7539OTUs.RData is the multi-domain skin microbiome data from the paper. 

In the real data, it contains:

``` Y1 ``` : bacterial microbiome count table. 60 samples $\times$ 75 OTUs. Each row is a sample, and each column is a bacterial OTU. 

``` Y2 ``` : viral microbiome count table. 60 samples $\times$ 39 OTUs. Each row is a sample, and each column is a viral OTU.   

``` Y ``` : combined multi-domain skin microbiome data.  60 samples $\times$ 114(75+39) OTUs.

``` X ``` : a categorical covariate representing experimental conditions. (1,0,0) pre-treatment & (0,1,0) post-treatment &  (0,0,1) healthy condition.

``` J ``` : number of OTUs in each domain. (75, 39)

``` Jsum ``` : number of total OTUs. 114

``` n ``` : number of samples. 60

``` S ``` : number of subjects. 20

### real-data-code

Real Data.R produces the results displayed in real data in the paper and please save the result as RData with Real Data.RData

### figures-code

Produces Figures 1-10 in the main text.

### figures

Contains the results from calling the code in the figures-code folder.








