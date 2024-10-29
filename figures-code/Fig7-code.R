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
library(ComplexHeatmap)
library(circlize)
library(metagenomeSeq)
col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))

######################### Fig 7(a) #########################

burn=1:nsamp
bb = matrix(NA, Jsum, nsamp)
for(i in 1:nsamp){
  bb[, i] = beta.st[,1,i] - beta.st[,2,i]
}
df.b = c(apply(bb, 1, function(x) quantile(x, 0.025)))
df.b = cbind(df.b, c(apply(bb, 1, median)))
df.b = cbind(df.b, c(apply(bb, 1, function(x) quantile(x, 0.975))))
df.b = cbind(df.b, c( beta_jp_truth[,1] - beta_jp_truth[,2]))
df.b = data.frame(df.b)
colnames(df.b) <- c('2.5%', 'Pos.Mean', '97.5%', 'True')
ggplot(df.b[1:J[1], ], aes(x = True)) +
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.1, linewidth =.2,
                colour = "light grey") +
  geom_point(aes(y = Pos.Mean, colour = "Pos.Mean"), size = 1.5) +
  geom_line(aes(x = True, y = True, colour = "True")) + 
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(TeX("$\\widehat{beta_{1j1}-beta_{1j2}}$"))+
  xlab(TeX("$beta^{TR}_{1j1}-beta^{TR}_{1j2}$")) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=25),
        text=element_text(size=25))

######################### Fig 7(b) #########################

ggplot(df.b[1:J[2]+J[1], ], aes(x = True)) +
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.1, linewidth =.2,
                colour = "light grey") +
  geom_point(aes(y = Pos.Mean, colour = "Pos.Mean"), size = 1.5) +
  geom_line(aes(x = True, y = True, colour = "True")) + 
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(TeX("$\\widehat{beta_{2j1}-beta_{2j2}}$"))+
  xlab(TeX("$beta^{TR}_{2j1}-beta^{TR}_{2j2}$")) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=25),
        text=element_text(size=25))

######################### Fig 7(c) #########################

Ymeta = newMRexperiment(t(Y[, 1:J[1]]))
Y.Trim = cumNorm(Ymeta)
intercept = rep(1,n)
mod <- model.matrix(~ -1 + X)
settings <- zigControl(maxit = 11, verbose = T)
fit <- fitZig(obj = Y.Trim, mod = mod, useCSSoffset = T, control = settings)
df.b = fit@fit$coefficients[,1] - fit@fit$coefficients[,2]
df.b = cbind(df.b, c( beta_jp_truth[1:J[1],1] - beta_jp_truth[1:J[1],2]))
df.b = data.frame(df.b)
colnames(df.b) <- c('Pos.Mean', 'True')
ggplot(df.b, aes(x = True)) +
  geom_point(aes(y = Pos.Mean, colour = "Pos.Mean"), size = 1.5) +
  geom_line(aes(x = True, y = True, colour = "True")) + 
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(TeX("$\\widehat{beta_{1j1}-beta_{1j2}}$"))+
  xlab(TeX("$beta^{TR}_{1j1}-beta^{TR}_{1j2}$")) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=25),
        text=element_text(size=25))

######################### Fig 7(d) #########################

Ymeta = newMRexperiment(t(Y[, 1:J[2] + J[1]]))
Y.Trim = cumNorm(Ymeta)
intercept = rep(1,n)
mod <- model.matrix(~ -1 + X)
zeromod <- model.matrix(~ -1 + X)
settings <- zigControl(maxit = 20, verbose = T)
fit <- fitZig(obj = Y.Trim, mod = mod, useCSSoffset = T, control = settings)
df.b = fit@fit$coefficients[,1] - fit@fit$coefficients[,2]
df.b = cbind(df.b, c( beta_jp_truth[1:J[2]+J[1],1] - beta_jp_truth[1:J[2]+J[1],2]))
df.b = data.frame(df.b)
colnames(df.b) <- c('Pos.Mean', 'True')
ggplot(df.b, aes(x = True)) +
  geom_point(aes(y = Pos.Mean, colour = "Pos.Mean"), size = 1.5) +
  geom_line(aes(x = True, y = True, colour = "True")) + 
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(TeX("$\\widehat{beta_{1j1}-beta_{1j2}}$"))+
  xlab(TeX("$beta^{TR}_{1j1}-beta^{TR}_{1j2}$")) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=25),
        text=element_text(size=25))
