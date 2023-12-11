######################### generate data #########################


#### load real data
rm(list = ls())
load("~/Desktop/Project 2/23-1-25/Subject simulation/75-39/Filtered7539OTUs.RData")
library(SpiecEasi)
seed = 3
set.seed(seed)
colnames(Y) <- paste0('OTU', 1:Jsum)

#### generate group 1 bacterial-like data 
Y1 = Y[, 1:(J[1])]
# generate cov 1
set.seed(seed)
true.cor <- diag(rep(1, ncol(Y1)))
vineBeta <- function(d, b){
  P = matrix(0, d, d)
  S = diag(d)
  
  for(k in 1:(d-1)){
    for(i in (k+1):d){
      P[k, i] = rbeta(1, b, b) 
      P[k, i] = (P[k, i] - 0.5) * 2
      if(abs(P[k ,i]) <0.8) {P[k ,i]=0}
      p = P[k, i]
      
      if(k == 1){
        p = p
      } else {
        for(l in (k-1):1 ){
          p = p * sqrt((1-P[l,i]^2)*(1-P[l,k]^2)) + P[l,i]*P[l,k]
        }
      }
      
      S[k, i] = p
      S[i, k] = p
    }
  }
  permutation = sample(1:d, d, replace = F)
  return(S)
}
set.seed(seed)
ttt = vineBeta(J[1], 1)
matrixcalc::is.positive.definite(ttt)
ttt = as.matrix(Matrix::nearPD(ttt)$mat)
corrplot::corrplot(ttt)
matrixcalc::is.positive.definite(ttt)
factoextra::fviz_eig(prcomp(ttt), addlabels = T, ncp = 50)
true.cor = ttt
# generate count table 1
depths <- rowSums(Y1)
Y.n  <- t(apply(Y1, 1, norm_to_total))
Y.cs <- round(Y.n * min(depths)) ### use min depth
n = nrow(Y1)
true.cor1 = true.cor
hist(true.cor1[true.cor1!=0 & true.cor1!=1])
# synth_comm_from_counts function is from SpiecEasi package
Y1 <- synth_comm_from_counts(Y.cs, mar=2, distr='zinegbin', Sigma=true.cor1, n=n, retParams = T)
mean(Y1==0)

#### generate group 2 bacterial-like data 
Y2 = Y[, 1:(J[2])+(J[1])]
# generate cov 2
set.seed(seed)
ttt = vineBeta(J[2], 1)
matrixcalc::is.positive.definite(ttt)
ttt = as.matrix(Matrix::nearPD(ttt)$mat)
corrplot::corrplot(ttt)
matrixcalc::is.positive.definite(ttt)
factoextra::fviz_eig(prcomp(ttt), addlabels = T, ncp = 50)
true.cor2 = ttt
hist(true.cor2[true.cor2!=0 & true.cor2!=1])
# generate count table 2
depths <- rowSums(Y2)
Y.n  <- t(apply(Y2, 1, norm_to_total))
Y.cs <- round(Y.n * min(depths))
Y2 <- synth_comm_from_counts(Y.cs, mar=2, distr='zinegbin', Sigma=true.cor2, n=n)
mean(Y2==0)

#### combine group 1,2 count data, cov matrix
Y <- cbind(Y1, Y2)
true.cor = cbind(true.cor1, matrix(0, J[1], J[2]))
true.cor = rbind(true.cor, cbind(matrix(0, J[2], J[1]), true.cor2))
corrplot::corrplot(true.cor)
factoextra::fviz_eig(prcomp(true.cor), addlabels = T, ncp = 50)
rm(list=setdiff(ls(), c("Y","true.cor", "J")))

######################### prior hyper-parameter #########################
seed = 3
s = 60
Jsum = ncol(Y)
n = nrow(Y)
S = rep(1:s, each = 1)
library(statmod)
library(GIGrvg)
library(extraDistr)

K = 25
a.sig = 3; b.sig = 3
a.phi = 1/20 
a.tau = 0.1; b.tau = 1/Jsum

sig.pro = matrix(1, Jsum, K)
nacc = true.nacc = matrix(0, Jsum, K)
acc = matrix(0, Jsum, K)
acc.tar = 0.234

######################### prior mean hyper-parameter #########################

m = 2
M = rep(c(1:m), J)
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
Lr = 30; a.psi.r = 1; a.w = b.w = 5;
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
L.alpha = 35; a.psi.alpha = 1; a.w.alpha = b.w.alpha = 5; 
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
niter = 150000
library(statmod)
library(GIGrvg)

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

sourceCpp("~/Desktop/Rcpp 2/update_GIG_Subject.cpp")
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
save.image("~/Desktop/All new sim/New sim 2/New sim 2 60 M2 mindepth.RData")

