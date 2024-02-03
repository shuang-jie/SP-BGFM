rm(list = ls())

########################## load packages #########################

library(LaplacesDemon)
library(latex2exp)
library(tikzDevice)

########################## plot 2(a) ############################
seed = 1

#a.phi = 1/2
#a.phi = 1/10

a.phi = 1/20; N = 10000
set.seed(seed)
gamma1 = rep(1, N)
phij1 = rbeta(N, a.phi, a.phi); phij2 = 1-phij1

set.seed(seed)
Zetaj1 = rhalfcauchy(N, 1) ; Zetaj2 = rhalfcauchy(N, 1) 
lambdaj_DH1 = rnorm(N, 0, sqrt(Zetaj1^2 * phij1)); lambdaj_DH2 = rnorm(N, 0, sqrt(Zetaj2^2 * phij2))

grids = seq(-8, 8, length = 3000); ngrids = length(grids) 

grid.m = matrix(grids, N, ngrids, byrow = T)
sd.m = matrix(sqrt(Zetaj1^2 * phij1 * gamma1), N, ngrids, byrow = F)
dens1 = dnorm(grid.m, 0, sd.m)
dens_eval1 = apply(dens1, 2, mean)

sd.m = matrix(sqrt(Zetaj2^2 * phij2 * gamma1), N, ngrids, byrow = F)
dens1.2 = dnorm(grid.m, 0, sd.m)
dens_eval1.2 = apply(dens1.2, 2, mean)

Z = replicate(ngrids, dens_eval1) * t(replicate(ngrids, dens_eval1.2))
Z = log(Z)
plot(lambdaj_DH1, lambdaj_DH2, pch = 16, cex = 0.5, xlim = c(-8, 8), ylim = c(-8, 8),
     main = '', xlab = expression(lambda[1]), cex.lab = 2, ylab="", cex.axis = 1.5)
title(ylab=expression(lambda[2]), line=2.35, cex.lab=2)
contour(grids, grids, Z, xlim = c(-8, 8), ylim = c(-8, 8), add = T, col = 'red', labcex=1.3)

########################## plot 2(b) ############################

set.seed(seed)
lambdaj_DL1 = LaplacesDemon::rlaplace(N, location = 0, scale = phij1)
phij2[phij2==0] = min(phij2[phij2!=0])
lambdaj_DL2 = LaplacesDemon::rlaplace(N, location = 0, scale = phij2)

sd3.m = matrix(phij1 * gamma1, N, ngrids, byrow = F)
dens3 = matrix(LaplacesDemon::dlaplace(grid.m, 0, scale = sd3.m), N, ngrids)
dens_eval3 = apply(dens3, 2, mean)

sd3.m = matrix(phij2 * gamma1, N, ngrids, byrow = F)
dens3.2 = matrix(LaplacesDemon::dlaplace(grid.m, 0, scale = sd3.m), N, ngrids)
dens_eval3.2 = apply(dens3.2, 2, mean)

Z.3 = replicate(ngrids, dens_eval3) * t(replicate(ngrids, dens_eval3.2))
Z.3 = log(Z.3)

plot(lambdaj_DL1, lambdaj_DL2, pch = 16, cex = 0.5, xlim = c(-8, 8), ylim = c(-8, 8),
     main = '', xlab = expression(lambda[1]), cex.lab = 2, ylab="", cex.axis = 1.5)
title(ylab=expression(lambda[2]), line=2.35, cex.lab=2)
contour(grids, grids, Z.3, xlim = c(-8, 8), ylim = c(-8, 8), add = T, col = 'red', labcex=1.4)

########################## plot 2(c) ############################

set.seed(seed)
lambdaj_H1 = rnorm(N, 0, sqrt(Zetaj1^2/2))
lambdaj_H2 = rnorm(N, 0, sqrt(Zetaj2^2/2))

sd2.m = matrix(sqrt(Zetaj1^2 * gamma1 * 1/2), N, ngrids, byrow = F)
dens2 = dnorm(grid.m, 0, sd2.m)
dens_eval2 = apply(dens2, 2, mean)

sd2.m = matrix(sqrt(Zetaj2^2 * gamma1 * 1/2), N, ngrids, byrow = F)
dens2.2 = dnorm(grid.m, 0, sd2.m)
dens_eval2.2 = apply(dens2.2, 2, mean)

Z.2 = replicate(ngrids, dens_eval2) * t(replicate(ngrids, dens_eval2.2))
Z.2 = log(Z.2)
plot(lambdaj_H1, lambdaj_H2, pch = 16, cex = 0.5, xlim = c(-8, 8), ylim = c(-8, 8),
     main = '', xlab = expression(lambda[1]), cex.lab = 2, ylab="", cex.axis = 1.5)
title(ylab=expression(lambda[2]), line=2.35, cex.lab=2)
contour(grids, grids, Z.2, xlim = c(-8, 8), ylim = c(-8, 8), add = T, col = 'red', labcex=1.4)

