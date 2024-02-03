######################################## 
############    !!note!!     ###########
###### Please first run Sim 2 and ###### 
###### save result as Sim 2.RData ###### 
######################################## 

################## load simulation 2 result ####################

load("./simulation-code/Sim 2.RData")

######################### load package #########################

library(ggplot2)
library(latex2exp)

######################### plot 8(a) ############################

mi=1
burn=(nsamp-5000):nsamp
Gmj = array(NA, dim = c(n, J[mi], length(burn)))
Xb.st <- array(0, dim = c(n, J[mi], length(burn)))
for(li in burn){
  Xb.st[,,li-burn[1]+1] <- tcrossprod(X, beta.st[J.ls[[mi]],,li]) 
}
Sigmali <- matrix(0, nrow = length(burn), ncol = J[mi])
for(li in burn){
  Sigmali[li-burn[1]+1,] <- diag(tcrossprod(Lambda.st[J.ls[[mi]],,li])) + sig2.st[mi, li]
}
for(li in burn){
  psi.alpha.pre = psi.alpha.st[mi, 1:L.alpha, li]
  w.alpha.pre = w.alpha.st[mi, 1:L.alpha, li]
  xi.a.pre = xi.alpha.st[1:J[1], 1:L.alpha, li]
  alpha.pre = matrix(0, s, J[mi])
  
  for(ss in 1:s){
    for(j in 1:(J[1])){
      Sij1.pre = sample(1:(L.alpha[1]), size = 1, prob = psi.alpha.pre)
      Sij2.pre = sample(c(1,0), size = 1, prob = c(w.alpha.pre[Sij1.pre], 1-w.alpha.pre[Sij1.pre]))
      alpha.pre[ss,j] = Sij2.pre * xi.a.pre[j, Sij1.pre] +
        (1-Sij2.pre) * ((nu.alpha[j] - w.alpha.pre[Sij1.pre] * xi.a.pre[j, Sij1.pre]) / (1-w.alpha.pre[Sij1.pre]))
    }
  }
  
  alphaij = alpha.pre[S, ]
  
  for(i in 1:n){
    for(j in 1:(J[1])){
      Gmj[i, j, li-burn[1]+1] = floor(exp(rnorm(1, alphaij[i,j] + Xb.st[i, j, li-burn[1]+1], sqrt(Sigmali[li-burn[1]+1, j]) ) ) )
    }
  }
}

pos.ri.mean1= rowMeans(ri.st[1,,])
Gmj.log = log(Gmj+1)
Gmj.mean = rowMeans(Gmj.log, dims = 2)
Gmj.low = apply(Gmj.log, 1:2, function(x) quantile(x,0.025))
Gmj.up = apply(Gmj.log, 1:2, function(x) quantile(x,0.975))
Gmj.median = apply(Gmj.log, 1:2, function(x) quantile(x,0.5))
Gmj.true =  log(floor(exp(log(Y[, J.ls[[mi]]]+1)- pos.ri.mean1))+1)

odds = seq(1, n, 2); even = seq(2, n, 2)
df1 = data.frame(x1 = c(Gmj.log[odds, 12,]))
df1_c2 = data.frame(x1 = c(Gmj.log[even, 12,]))
set.seed(1)
ggplot() +
  geom_density(aes(x = x1), df1,  alpha = 0.3, adjust = 1) +
  geom_density(aes(x = x1), df1_c2,  alpha = 0.3, adjust = 1, col = 'red',
               lty = 2) +
  xlim(-2, 17) +
  ylim(0, .39) +
  xlab(TeX("$\\log(Y_{i1j}+1)$")) +
  ylab("Posterior Predictive Density") +
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=25),
        text=element_text(size=25)) + 
  theme(panel.grid.minor = element_blank())+ 
  geom_point(aes(x = jitter(Gmj.true[odds, 12], 10) , y = 0.35), 
             col = "black", shape = 20, size =5)+ 
  geom_point(aes(x = jitter(Gmj.true[even, 12], 10) , y = 0.38), 
             col = "red", shape = 4, size =5)

######################### plot 8(b) ############################

odds = seq(1, n, 3); df2 = data.frame(x2 = c(Gmj.log[odds, 32,]))
df2_c2 = data.frame(x2 = c(Gmj.log[even, 32,]))
set.seed(2)
ggplot() +
  geom_density(aes(x = x2), df2,  alpha = 0.3, adjust = 1) +
  geom_density(aes(x = x2), df2_c2,  alpha = 0.3, adjust = 1, col = 'red',
               lty = 2) +
  xlim(-2, 15) +
  ylim(0, 0.7) +
  xlab(TeX("$\\log(Y_{i1j}+1)$")) +
  ylab("Posterior Predictive Density") +
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=25),
        text=element_text(size=25)) + 
  theme(panel.grid.minor = element_blank())+ 
  geom_point(aes(x = jitter(Gmj.true[odds, 32], 10) , y = 0.65), 
             col = "black", shape = 20, size =5)+ 
  geom_point(aes(x = jitter(Gmj.true[even, 32], 10) , y = 0.70), 
             col = "red", shape = 4, size =5)


######################### plot 8(c) ############################

mi=2
Gmj2 = array(NA, dim = c(n, J[mi], length(burn)))
Xb.st <- array(0, dim = c(n, J[mi], length(burn)))
for(li in burn){
  Xb.st[,,li-burn[1]+1] <- tcrossprod(X, beta.st[J.ls[[mi]],,li]) 
}
Sigmali <- matrix(0, nrow = length(burn), ncol = J[mi])
for(li in burn){
  Sigmali[li-burn[1]+1,] <- diag(tcrossprod(Lambda.st[J.ls[[mi]],,li])) + sig2.st[mi, li]
}
for(li in burn){
  psi.alpha.pre = psi.alpha.st[mi, 1:L.alpha, li]
  w.alpha.pre = w.alpha.st[mi, 1:L.alpha, li]
  xi.a.pre = xi.alpha.st[J.ls[[mi]], 1:L.alpha, li]
  alpha.pre = matrix(0, s, J[mi])
  
  for(ss in 1:s){
    for(j in 1:(J[2])){
      Sij1.pre = sample(1:(L.alpha[1]), size = 1, prob = psi.alpha.pre)
      Sij2.pre = sample(c(1,0), size = 1, prob = c(w.alpha.pre[Sij1.pre], 1-w.alpha.pre[Sij1.pre]))
      alpha.pre[ss,j] = Sij2.pre * xi.a.pre[j, Sij1.pre] +
        (1-Sij2.pre) * ((nu.alpha[j+J[1]] - w.alpha.pre[Sij1.pre] * xi.a.pre[j, Sij1.pre]) / (1-w.alpha.pre[Sij1.pre]))
    }
  }
  
  alphaij = alpha.pre[S, ]
  
  for(i in 1:n){
    for(j in 1:(J[2])){
      Gmj2[i, j, li-burn[1]+1] = floor(exp(rnorm(1, alphaij[i,j] + Xb.st[i, j, li-burn[1]+1], sqrt(Sigmali[li-burn[1]+1, j]) ) ) )
    }
  }
}

pos.ri.mean2= rowMeans(ri.st[2,,])
Gmj2.log = log(Gmj2+1)
Gmj2.mean = rowMeans(Gmj2.log, dims = 2)
Gmj2.low = apply(Gmj2.log, 1:2, function(x) quantile(x,0.025))
Gmj2.up = apply(Gmj2.log, 1:2, function(x) quantile(x,0.975))
Gmj2.median = apply(Gmj2.log, 1:2, function(x) quantile(x,0.5))
Gmj2.true =  log(floor(exp(log(Y[, J.ls[[mi]]]+1)- pos.ri.mean2))+1)

j3 = 161-J[1]
df3 = data.frame(x3 = c(Gmj2.log[odds, j3,]))
even = seq(2,n, 2)
df3_c2 = data.frame(x3 = c(Gmj2.log[even, j3,]))
set.seed(1)
ggplot() +
  geom_density(aes(x = x3), df3,  alpha = 0.3, adjust = 2) +
  geom_density(aes(x = x3), df3_c2,  alpha = 0.3, adjust = 2, col = 'red',
               lty = 2) +
  xlim(-2, 17) +
  ylim(0, 1.63) +
  xlab(TeX("$\\log(Y_{i1j}+1)$")) +
  ylab("Posterior Predictive Density") +
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=25),
        text=element_text(size=25)) + 
  theme(panel.grid.minor = element_blank())+ 
  geom_point(aes(x = jitter(Gmj2.true[odds, j3], 1) , y = 1.55), 
             col = "black", shape = 20, size =5)+ 
  geom_point(aes(x = jitter(Gmj2.true[even, j3], 5) , y = 1.63), 
             col = "red", shape = 4, size =5)



