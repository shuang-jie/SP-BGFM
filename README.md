# SP-BGFM
SP-BGFM (Sparse Bayesian Group Factor Model for Feature Interactions in Multiple Count Tables Data) is a Bayesian semi-parametric model for multiple domains of next-generation sequencing microbiome abundance data. It models cross-domain interaction between microbial features directly using observed count tables, and provides a flexible DP mixture structure for interpretation.

For more information, read the paper: Sparse Bayesian Group Factor Model for Feature Interactions in Multiple Count Tables Data (to appear). 

Contact: Shuangjie Zhang (szhan209 AT ucsc DOT edu)


## Installation

SP-BGFM requires some packages, especially Rcpp and RcppArmadillo for Rc++ functions: 

```
install.packages(c("Rcpp", "RcppArmadillo", "statmod", "GIGrvg", "extraDistr", "abind", "mvnfast", "mvnfast", "statmod", "extraDistr"))
```

Please make sure the C++ compiler is correctly installed. 

## Simulation 

We present the code of Sim 1 in the papaer. We need another a few packages installed to simulate the synthentic dataset. 

```
install.packages(c("mvnfast", "statmod", "GIGrvg", "extraDistr", "abind"))
```


```
rm(list = ls())
m = 2
n = 20
s = 20
S = rep(1:s, each = 1)
J = c(150, 50)
M = rep(c(1:m), J)
Jsum = sum(J)
K.true = 5
seed = 1

########################## truth Omega setting ########################

set.seed(seed)
Lambda.true = matrix(0, nrow = Jsum, ncol =K.true)
sigma2.true = rep(0.25, m)

zero.rate=0.8
set.seed(seed)
for(j in 1:25)
{
  num=sample(1:5, size = 1, replace = F)
  ttt = extraDistr::rtnorm(1, 0, 1, -1, 1)
  ttt = ttt + (ttt>0) * 1 + (ttt<0) * (-1) 
  Lambda.true[j, num] =  ttt
}

set.seed(seed+1)
for(j in 1:25+50)
{
  num=sample(1:5, size = 1, replace = F)
  ttt = extraDistr::rtnorm(1, 0, 1, -1, 1)
  ttt = ttt + (ttt>0) * 1 + (ttt<0) * (-1) 
  Lambda.true[j, num] =  ttt
}

set.seed(seed+2)
for(j in 1:25+150)
{
  num=sample(1:5, size = 1, replace = F)
  ttt = extraDistr::rtnorm(1, 0, 1, -1, 1)
  ttt = ttt + (ttt>0) * 1 + (ttt<0) * (-1) 
  Lambda.true[j, num] =  ttt
}

Omega.true = Lambda.true %*% t(Lambda.true) + diag(rep(sigma2.true, J))
true.cor <- cov2cor(Omega.true)

########################## truth ri alphai setting ########################

library(mvnfast)
location.mean.alpha.true = matrix(0, Jsum, 3)
set.seed(seed)
Weights = matrix(0, Jsum, 3)
for(j in 1:Jsum){
  Weights[j,] = extraDistr::rdirichlet(1, alpha=c(30, 40, 30))
}

location.mean.alpha.true[,1] = -5
set.seed(seed)
for(j in 1:Jsum){
  location.mean.alpha.true[j,2] = rnorm(1, 4, 1)
  location.mean.alpha.true[j,3] = rnorm(1, 10, 1)
}

alphaSij.true = matrix(0, n, Jsum)

for(ii in 1:n){
  for(jj in 1:Jsum){
    try = sample(1:3, 1, prob = Weights[jj,])
    alphaSij.true[ii, jj] = location.mean.alpha.true[jj, try]    
  }
}

alphaij.true = alphaSij.true[S, ]
ri.true = matrix(runif(m*n, 0, 2), m, n)

########################## truth mu setting ########################

mu.true = matrix(0, n, Jsum)
for(mi in 1:m){
  Jcol = which(M==mi)
  for(i in 1:n){
    for(j in Jcol){
      mu.true[i, j] = ri.true[mi, i] + alphaij.true[i, j]
    }
  }
}

########################## simulate data #######################

y.star.true = matrix(NA, nrow = n, ncol = Jsum)
set.seed(seed)

for(i in 1:n){
  y.star.true[i, ] = rmvn(1, mu.true[i, ], Omega.true)
}

Y = floor(exp(y.star.true))
sum(Y==0)/(n*Jsum)
max(log(Y))

normalized.count = Y

for(i in 1:n){
  normalized.count[i, ] = normalized.count[i, ]/exp(ri.true[,i])
}

######################### prior hyper-parameter #########################

library(statmod)
library(GIGrvg)
library(extraDistr)

K = 10

a.sig = 3; b.sig = 3
a.phi = 1/20 
a.tau = 0.1; b.tau = 1/Jsum

sig.pro = matrix(1, Jsum, K)
nacc = true.nacc = matrix(0, Jsum, K)
acc = matrix(0, Jsum, K)
acc.tar = 0.234

######################### prior mean hyper-parameter #########################

hat.ri = matrix(NA, nrow = m, ncol = n)
for(mi in 1:m){
  hat.ri[mi,] = rowSums(log(Y[,(1:J[mi])+sum(J[1:(mi-1)])* (mi!=1) ]+0.01)/J[mi])
}

hat.alpha = c()
for(mi in 1:m){
  hat.alpha[(1:J[mi])+sum(J[1:(mi-1)])* (mi!=1)] = 
    colMeans( log(Y[,(1:J[mi])+sum(J[1:(mi-1)])* (mi!=1)]+0.01) - matrix(hat.ri[mi,], n, J[mi]) )
}

hat.alphasij = matrix(0, s, Jsum)

for(i in 1:s){
  for(j in 1:Jsum){
    hat.alphasij[i, j] = mean(log(Y[which(S==i),j]+0.01) - hat.ri[M[j],i])
  }
}


### ri hyper-parameter ###
nu.r = rowMeans(hat.ri)
Lr = 30; a.psi.r = 1; a.w = b.w = 1;
ur2 = 1
a.xi = nu.r
ri = hat.ri
Si1 = matrix(NA, m, n)
for(mi in 1:m){
  Si1[mi, ] = sample(1:Lr, n, replace = T)
}
Si2 = matrix(sample(0:1, m * n, replace = T), nrow = m, ncol = n)
set.seed(seed)
w.l.r = lapply(1:m, function(x) rbeta(Lr, a.w, b.w))
V.r = lapply(1:m, function(x){
  res <- rbeta(Lr-1, 1, a.psi.r)
  return(res)
})
V.recover.psi = function(x){
  log_x = log(x)
  log_x[is.infinite(log_x)] = -.Machine$double.xmax
  res <- c()
  res[1] <- x[1]
  
  for(i in 2:(length(x))){
    res[i] = exp( log_x[i] + sum(log(1- x[1:(i-1)])))
  }
  res[length(x)+1] = exp(sum(log(1-x)))
  return(res)
}
psi.r =  lapply(1:m, function(x) V.recover.psi(V.r[[x]]))
psi.r = matrix(unlist(psi.r), m, Lr[1], byrow = T)
w.l.r = matrix(unlist(w.l.r), m, Lr[1], byrow = T)
xi = matrix(0, m, Lr[1])

### alphai hyper-parameter ###
nu.alpha = hat.alpha
L.alpha = 35; a.psi.alpha = 1; a.w.alpha = b.w.alpha = 1; 
u2.alpha = Rfast::colVars(hat.alphasij)
a.xi.alpha = nu.alpha
set.seed(seed)
alphasij = hat.alphasij
alphaij = alphasij[S, ]
Sij1 = matrix(NA, s, Jsum)
set.seed(seed)
for(i in 1:s){
  for(j in 1:Jsum){
    Sij1[i, j] = sample(1:L.alpha, 1)
  }
}
set.seed(seed)
Sij2 = matrix(1, s, Jsum)
xi.alpha = matrix(rnorm(Jsum * L.alpha, nu.alpha, sqrt(u2.alpha)), Jsum, L.alpha)
w.alpha = lapply(1:m, function(x) rbeta(L.alpha, a.w.alpha, b.w.alpha))
V.alpha <- lapply(1:m, function(x){rbeta(L.alpha-1, 1, a.psi.alpha)})
psi.alpha = lapply(1:m, function(x) V.recover.psi(V.alpha[[x]]))

w.alpha = matrix(unlist(w.alpha), m, L.alpha, byrow = T)
psi.alpha = matrix(unlist(psi.alpha), m, L.alpha, byrow = T)

######################### start point #########################
library(mvnfast)
ls = list()

phi.m = matrix(0, Jsum, K)
til.phi.m = til.til.phi.m = matrix(0, Jsum, K)
for(ki in 1:K){
  for(mi in 1:m){
    Jcol = which(M==mi)
    try = rgamma(J[mi], a.phi, 1)
    phi.m[Jcol, ki] = try/sum(try)
    til.phi.m[Jcol, ki] = try
    til.til.phi.m[Jcol, ki] = log(try)
  }
}

tau = rep(1, K)

set.seed(seed)
sig2 = rep(1, m)
sig2.m = matrix(rep(sig2, J), n, Jsum, byrow = T)
Dj = diag(rep(1/sig2, J))

eta = rmvn(n, rep(0, K), diag(1, K))
Lambda = matrix(0, nrow = Jsum, ncol =K)

ri = hat.ri
ri.muij = t(ri)[,rep(1:m, J)]
y.star = log(Y+0.01)
Z = matrix(1/rgamma(Jsum *K, 1/2, 1), Jsum, K)
zeta = matrix(1/rgamma(Jsum *K, 1/2, 1/Z), Jsum, K)

######################### likelihood function #########################

pos.phi = function(t.phi, wj){
  ttt = sqrt(zeta[, ki] * t.phi/sum(t.phi) * tau[ki])
  ttt[ttt==0] = sort(ttt[ttt!=0])[1]
  dgamma(t.phi[wj], a.phi, 1, log = T) + sum(dnorm(Lambda[, ki], 0, ttt, log=T)) + log(t.phi[wj])
} 

######################### run MCMC #########################
niter = 100000

sig2.xi.r = rep(1, m)
log.Y = log(Y)
log.Y1 = log(Y+1)
J.ls = lapply(1:m, function(x) which(M==x))
nsamp = niter/10

count.st = 0
ri.st <- array(NA, dim=c(m, n, nsamp))
#alphaij.st <- array(NA, dim=c(n, Jsum, nsamp))
alphasij.st <- array(NA, dim=c(s, Jsum, nsamp))

xi.alpha.st <- array(NA, dim=c(Jsum, L.alpha[1], nsamp))
psi.alpha.st <- array(NA, dim=c(m, L.alpha[1], nsamp))
w.alpha.st <- array(NA, dim=c(m, L.alpha[1], nsamp))

Lambda.st <- array(NA, dim=c(Jsum, K, nsamp))
tau.st <- matrix(NA, K, nsamp)
phi.m.st <- array(NA, dim=c(Jsum, K, nsamp))
sig2.st <- matrix(NA, m, nsamp)
Sij1.st <- array(NA, dim=c(s, Jsum, nsamp))
Sij2.st <- array(NA, dim=c(s, Jsum, nsamp))
eta.st <- array(NA, dim=c(n, K, nsamp))
RSS = y.star - ri.muij - alphaij - tcrossprod(eta, Lambda)

library(Rcpp)
library(RcppArmadillo)
library(abind)

J_bound = J
J_bound[1] = J_bound[1] -1
J_bound = c(0, cumsum(J_bound))

sourceCpp("~/update_github.cpp")

time.stamp.0 <- proc.time()
for(ni in 1:niter){
  # impute continuous latent variable
  RSS = y.star - RSS
  y.star = update_ystar(log.Y, log.Y1, n, Jsum, RSS, sig2, M-1)
  RSS = y.star - RSS
  
  ## update sig2
  sig2 = update_sig2(m, n, RSS, J_bound, J, a.sig, b.sig)
  Dj = diag(rep(1/sig2, J))
  
  # update eta_i
  RSS = RSS + tcrossprod(eta, Lambda)
  eta = update_eta(K, Lambda, Dj, n, RSS)
  
  ## update lambda_j^(m)
  Lambda = update_Lambda(Jsum, K, zeta, Z, tau, phi.m, eta, RSS, sig2, M-1)
  RSS = RSS - tcrossprod(eta, Lambda)
  
  ## update zeta
  try = update_Z_Zeta(Jsum, K, Lambda, phi.m, tau, Z, zeta)
  zeta = try[,,1]; Z = try[,,2]
  
  ## M-H: phi.m
  delta = min(0.01, 1/sqrt(ni))
  try = update_phi(K, Jsum, til.phi.m, sig.pro, a.phi, Lambda, zeta, tau, acc, phi.m, nacc, true.nacc, ni, acc.tar, delta)
  til.phi.m = try[,,1]
  phi.m = try[,,2]
  nacc = try[,,3]
  true.nacc = try[,,4]
  sig.pro = try[,,5]
  
  ##  tau.k
  tau = update_tau_GIG(K, a.tau, b.tau, Jsum, Lambda, phi.m, zeta)
  
  ## update ri related
  RSS = RSS + ri.muij
  
  try = update_psi_w_r(m, Lr, Si1-1, Si2, a.w, b.w, a.psi.r)
  psi_r_m = adrop(try[,,1,drop = F], drop = 3);
  w_l_r_m = adrop(try[,,2,drop = F], drop = 3);
  
  try = update_Si12(Lr, n, m, ri, xi, ur2, nu.r,w_l_r_m, psi_r_m)
  Si1 = adrop(try[,,1,drop = F], drop = 3) + 1
  Si2 = adrop(try[,,2,drop = F], drop = 3)
  
  xi = update_xi(m, Lr, ri, Si1-1, Si2, a.xi, sig2.xi.r, w_l_r_m, psi_r_m, nu.r, ur2)
  ri = update_ri(m, n, ur2, J, sig2, RSS, J_bound, Si1-1, Si2, nu.r, xi, w_l_r_m)
  
  ri.muij = t(ri)[,rep(1:m, J)]
  RSS = RSS - ri.muij
  
  ## update alphai
  try = update_psi_w(m, L.alpha, J_bound = J_bound, Sij1-1, Sij2, a.w.alpha, 
                     b.w.alpha, a.psi.alpha)
  psi.alpha = adrop(try[,,1,drop = F], drop = 3)
  w.alpha = adrop(try[,,2,drop = F], drop = 3)
  
  RSS = RSS + alphaij
  try = update_Sij12(RSS, s, Jsum, L.alpha, xi.alpha, nu.alpha, w.alpha, psi.alpha, M-1,
                     sig2, S-1)
  Sij1 = try[,,1]+1
  Sij2 = try[,,2]
  # RSS = RSS - alphaij
  
  # RSS = RSS + alphaij
  IIJ1 = Sij1[S, ]
  IIJ2 = Sij2[S, ]
  xi.alpha = update_xi_alpha(L.alpha, Jsum, RSS, IIJ1-1, IIJ2, u2.alpha, a.xi.alpha, M-1, w.alpha, sig2, nu.alpha)  
  for(ii in 1:s){
    for(j in 1:Jsum){
      alphasij[ii, j] = Sij2[ii, j] * xi.alpha[j, Sij1[ii, j]] +
        (1-Sij2[ii, j]) * ((nu.alpha[j] - w.alpha[M[j], Sij1[ii, j]] * xi.alpha[j, Sij1[ii, j]]) / (1-w.alpha[M[j], Sij1[ii, j]] ))
    }
  }
  alphaij = alphasij[S, ]
  RSS = RSS - alphaij
  
  ## save res
  if((ni > niter/2) & (ni%%5 ==0)){
    count.st = count.st + 1
    ri.st [,, count.st] <- ri
    alphasij.st[,, count.st] <- alphasij
    xi.alpha.st[,, count.st] <- xi.alpha
    psi.alpha.st[, ,count.st] <- psi.alpha
    w.alpha.st[, ,count.st] <- w.alpha
    Sij1.st[,, count.st] <- Sij1
    Sij2.st[,, count.st] <- Sij2
    sig2.st[, count.st] <- sig2
    
    Lambda.st[,, count.st] <- Lambda
    tau.st[,count.st] <- tau
    phi.m.st[,, count.st] <- phi.m
    eta.st[,, count.st] <- eta
    print(ni)
  }
}

run.time <- proc.time() - time.stamp.0
run.time


```
