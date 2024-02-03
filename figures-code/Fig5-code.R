######################################## 
############    !!note!!     ###########
###### Please first run Sim 1 and ###### 
###### save result as Sim 1.RData ###### 
######################################## 

################## load simulation 1 result ####################

load("./simulation-code/Sim 1.RData")

######################### load package #########################

library(ggplot2)
library(latex2exp)

######################### plot 5(a) ############################

mi=1
burn=(nsamp-5000):nsamp
Gmj = array(NA, dim = c(n, J[mi], length(burn)))
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
      Gmj[i, j, li-burn[1]+1] = floor(exp(rnorm(1, alphaij[i,j] , sqrt(Sigmali[li-burn[1]+1, j]) ) ) )
    }
  }
}

pos.ri.mean1= rowMeans(ri.st[1,,])
Gmj.log = log(Gmj+1)
Gmj.mean = rowMeans(Gmj.log, dims = 2)
Gmj.low = apply(Gmj.log, 1:2, function(x) quantile(x,0.025))
Gmj.up = apply(Gmj.log, 1:2, function(x) quantile(x,0.975))
Gmj.median = apply(Gmj.log, 1:2, function(x) quantile(x,0.5))
Gmj.true = log(Y[, J.ls[[mi]]]+1)
for(i in 1:n){
  for(j in 1:(J[1])){
    if(Y[i, j] != 0){
      Gmj.true[i, j] = log(floor(Y[i, j]/exp(pos.ri.mean1)[i])+1)
    }
  }
}

set.seed(2)
df1 = data.frame(x = c(Gmj.log[, 30,]))
ggplot() +
  geom_density(aes(x = x),df1,  alpha = 0.3) +
  xlim(-1, 12) +
  ylim(0, 0.8) +
  xlab(TeX("$\\log(Y_{i1j}+1)$")) +
  ylab("Posterior Predictive Density") +
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=25),
        text=element_text(size=25)) + 
  theme(panel.grid.minor = element_blank())+ 
  geom_point(aes(x = jitter(Gmj.true[, 30], 30) , y = 0), 
             col = "red", shape = 4, size = 5)

######################### plot 5(b) ############################

df4 = data.frame(x4 = c(Gmj.log[, 133,]))
set.seed(2)
ggplot() +
  geom_density(aes(x = x4), df4,  alpha = 0.3, adjust = 1.5) +
  xlim(-1, 12) +
  ylim(0, 0.3) +
  xlab(TeX("$\\log(Y_{i1j}+1)$")) +
  ylab("Posterior Predictive Density") +
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=25),
        text=element_text(size=25)) + 
  theme(panel.grid.minor = element_blank())+ 
  geom_point(aes(x = jitter(Gmj.true[, 133], 30) , y = 0), 
             col = "red", shape = 4, size = 5)


######################### plot 5(c) ############################

mi=2
Gmj = array(NA, dim = c(n, J[mi], length(burn)))
Sigmali <- matrix(0, nrow = length(burn), ncol = J[mi])
for(li in burn){
  Sigmali[li-burn[1]+1,] <- diag(tcrossprod(Lambda.st[J.ls[[mi]],,li])) + sig2.st[mi, li]
}
for(li in burn){
  psi.alpha.pre = psi.alpha.st[mi, 1:L.alpha, li]
  w.alpha.pre = w.alpha.st[mi, 1:L.alpha, li]
  xi.a.pre = xi.alpha.st[1:J[2]+J[1], 1:L.alpha, li]
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
      Gmj[i, j, li-burn[1]+1] = floor(exp(rnorm(1, alphaij[i, j] , sqrt(Sigmali[li-burn[1]+1, j]) ) ) )
    }
  }
}

pos.ri.mean2= rowMeans(ri.st[2,,])
Gmj.log = log(Gmj+1)
Gmj.mean = rowMeans(Gmj.log, dims = 2)
Gmj.low = apply(Gmj.log, 1:2, function(x) quantile(x,0.025))
Gmj.up = apply(Gmj.log, 1:2, function(x) quantile(x,0.975))
Gmj.median = apply(Gmj.log, 1:2, function(x) quantile(x,0.5))
Gmj.true = log(Y[, J.ls[[mi]]]+1)
for(i in 1:n){
  for(j in 1:(J[2])){
    if(Y[i, j+J[1]] != 0){
      Gmj.true[i, j] = log(floor(Y[i, j+J[1]]/exp(pos.ri.mean2)[i])+1)
    }
  }
}

df4 = data.frame(x4 = c(Gmj.log[, 31,]))
set.seed(2)
ggplot() +
  geom_density(aes(x = x4), df4,  alpha = 0.3, adjust = 2) +
  xlim(-3, 12) +
  ylim(0, 0.4) +
  xlab(TeX("$\\log(Y_{i2j}+1)$")) +
  ylab("Posterior Predictive Density") +
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=25),
        text=element_text(size=25)) + 
  theme(panel.grid.minor = element_blank())+ 
  geom_point(aes(x = jitter(Gmj.true[, 31], 80) , y = 0), 
             col = "red", shape = 4, size = 5)
