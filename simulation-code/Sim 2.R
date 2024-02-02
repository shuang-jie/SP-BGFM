rm(list = ls())

######################### load real data and obtain empirical zero rate ########

load("Filtered7539OTUs.RData")
Y1 = Y[, 1:J[1]]; Y2 = Y[, 1:J[2] + J[1]]
zero.rate1 = apply(Y1, 2, function(x) sum(x==0)/n)
zero.rate2 = apply(Y2, 2, function(x) sum(x==0)/n)
rm(list=setdiff(ls(), c("zero.rate1","zero.rate2")))

########################## load packages #########################

library(statmod)
library(GIGrvg)
library(extraDistr)
library(Rcpp)
library(RcppArmadillo)
library(abind)
library(mvnfast)
library(statmod)
library(extraDistr)

########################## load update function ##################

sourceCpp("update_GIG_Subject.cpp")
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
pos.phi = function(t.phi, wj){
  ## likelihood function of Dir-HS prior for phi
  ttt = sqrt(zeta[, ki] * t.phi/sum(t.phi) * tau[ki])
  ttt[ttt==0] = sort(ttt[ttt!=0])[1]
  dgamma(t.phi[wj], a.phi, 1, log = T) + sum(dnorm(Lambda[, ki], 0, ttt, log=T)) + log(t.phi[wj])
}
vineBeta <- function(d, b){
  P = matrix(0, d, d)
  S = diag(d)
  ## vine method to generate random covariance
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

####################### set up dimension in the paper ############

m = 2; n = 40; s = 20; K.true = 5; seed = 2
J = c(150, 50); Jsum = sum(J)
S = rep(1:s, each = 2)
M = rep(c(1:m), J)

####################### set up covariance matrix truth ###########

set.seed(seed)
sigma2.true = rep(1, m)

set.seed(seed)
ttt = vineBeta(Jsum, 1)
ttt = as.matrix(Matrix::nearPD(ttt)$mat)
true.cor = ttt

sigma2.true = runif(Jsum, 1, 1.5)
Omega.true <- diag(rep(1, Jsum))
for(i in 1:Jsum){
  for(j in 1:Jsum){
    Omega.true[i,j] <- 
      true.cor[i, j] * sqrt( sigma2.true[i] * sigma2.true[j] )
  }
}

####################### set up ri alpha_Sij truth ################

ri.true = matrix(runif(m*n, 0, 2), m, n)

location.mean.alpha.true = matrix(0, Jsum, 3)
set.seed(seed)
Weights = matrix(0, Jsum, 3)
id1 = sample(1:length(zero.rate1), J[1], replace = T)
for(j in 1:J[1]){
  ww = zero.rate1[id1[j]]
  if(ww==0){ww = 0.0000001}
  Weights[j,] = extraDistr::rdirichlet(1, alpha=c(ww*100, (1-ww)*0.6*100,  (1-ww)*0.4*100))
}
set.seed(seed)
id2 = sample(1:length(zero.rate2), J[2], replace = T)
for(j in 1:J[2]){
  ww = zero.rate2[id2[j]]
  if(ww==0){ww = 0.0000001}
  Weights[j+J[1],] = extraDistr::rdirichlet(1, alpha=c(ww*100, (1-ww)*0.6*100,  (1-ww)*0.4*100))
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

####################### set up beta truth ########################

set.seed(seed)
X = cbind( rep(c(1,0), n/2), rep(c(0,1), n/2))
p = ncol(X)
beta_jp_truth <- matrix(0, nrow = Jsum, ncol = p)
for(j in 1:Jsum){
  delta = LaplacesDemon::rbern(1, 0.8)
  if(delta == 1){
    beta_jp_truth[j, 1] = 0 
  } else {
    c_v= rnorm(1, 0, sqrt(1/3)) 
    if(c_v >0){
      c_v = c_v + 1
    } else {
      c_v = c_v - 1
    }
    beta_jp_truth[j, 1] = c_v
  }
}
beta_jp_truth[,2] = 0

########################## set up mu(mean) truth ########################

mu.true = matrix(0, n, Jsum)
for(mi in 1:m){
  Jcol = which(M==mi)
  for(i in 1:n){
    for(j in Jcol){
      mu.true[i, j] = ri.true[mi, i] + alphaij.true[i, j]
    }
  }
}
mu.true = mu.true + tcrossprod(X, beta_jp_truth)

########################## simulate count data #######################

y.star.true = matrix(NA, nrow = n, ncol = Jsum)
set.seed(seed)
for(i in 1:n){
  y.star.true[i, ] = rmvn(1, mu.true[i, ], Omega.true)
}
Y = floor(exp(y.star.true))

######################### prior hyper-parameter #########################

K = 15; a.sig = 3; b.sig = 3; a.phi = 1/20 
a.tau = 0.1; b.tau = 1/Jsum; acc.tar = 0.234

######################### prior mean hyper-parameter ####################

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
Lr = 30; a.psi.r = 1; a.w = b.w = 5; ur2 = 1; a.xi = nu.r
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
psi.r =  lapply(1:m, function(x) V.recover.psi(V.r[[x]]))
psi.r = matrix(unlist(psi.r), m, Lr[1], byrow = T)
w.l.r = matrix(unlist(w.l.r), m, Lr[1], byrow = T)
xi = matrix(0, m, Lr[1])
sig2.xi.r = rep(1, m)

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

sig.pro = matrix(1, Jsum, K)
nacc = true.nacc = matrix(0, Jsum, K)
acc = matrix(0, Jsum, K)

ri = hat.ri
ri.muij = t(ri)[,rep(1:m, J)]
y.star = log(Y+0.01)
Z = matrix(1/rgamma(Jsum *K, 1/2, 1), Jsum, K)
zeta = matrix(1/rgamma(Jsum *K, 1/2, 1/Z), Jsum, K)

log.Y = log(Y)
log.Y1 = log(Y+1)
J.ls = lapply(1:m, function(x) which(M==x))

J_bound = J
J_bound[1] = J_bound[1] -1
J_bound = c(0, cumsum(J_bound))

mu_beta = matrix(0, nrow = p, ncol=1)
u2_beta = 3
beta_jp <- matrix(0, nrow = Jsum, ncol = p)

RSS = y.star - ri.muij - alphaij - tcrossprod(eta, Lambda)  - tcrossprod(X, beta_jp)

######################### create storage ###################

ls = list()
niter = 100000
nsamp = niter/10
count.st = 0

ri.st <- array(NA, dim=c(m, n, nsamp))
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
beta.st <- array(NA, dim=c(Jsum, p, nsamp))

######################### run MCMC #########################

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
  
  ## update alphai related
  try = update_psi_w(m, L.alpha, J_bound = J_bound, Sij1-1, Sij2, a.w.alpha, 
                     b.w.alpha, a.psi.alpha)
  psi.alpha = adrop(try[,,1,drop = F], drop = 3)
  w.alpha = adrop(try[,,2,drop = F], drop = 3)
  
  RSS = RSS + alphaij
  try = update_Sij12(RSS, s, Jsum, L.alpha, xi.alpha, nu.alpha, w.alpha, psi.alpha, M-1,
                     sig2, S-1)
  Sij1 = try[,,1]+1
  Sij2 = try[,,2]
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
  
  ## update beta
  RSS = RSS + tcrossprod(X, beta_jp)
  for(ji in 1:Jsum){
    U.l = chol(diag(rep(1/u2_beta, p), nrow = p, ncol = p) + crossprod(X) /sig2[M[ji]] )
    a.l = crossprod(X, RSS[, ji]) / sig2[M[ji]]
    beta_jp[ji, ] = backsolve(U.l, backsolve(U.l, a.l, transpose = T) + rnorm(p))
  }
  RSS = RSS - tcrossprod(X, beta_jp)
  
  ## save result
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
    beta.st[,, count.st] <- beta_jp
    
    Lambda.st[,, count.st] <- Lambda
    tau.st[,count.st] <- tau
    phi.m.st[,, count.st] <- phi.m
    eta.st[,, count.st] <- eta
  }
}

save.image("Sim 2.RData")


