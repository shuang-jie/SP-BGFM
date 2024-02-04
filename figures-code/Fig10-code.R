######################################## 
############    !!note!!     ###########
###### Please first Real Data.R and ####
###### save result as Real Data.RData ##
######################################## 

################## load Real Data result ######################

load("./real-data/Real Data.RData")

######################### load package #########################

library(ggplot2)
library(latex2exp)

######################### plot 10(a) ############################

burn=1:nsamp
pos.mean.beta = rowMeans(beta.st[,,burn, drop=F], dims = 2)
bb1 = matrix(NA, Jsum, nsamp)
for(i in 1:nsamp){
  bb1[, i] = beta.st[,1,i] - beta.st[,3,i]
}
pos.low.bb1 = apply(bb1, 1, function(x) quantile(x, 0.025))
pos.median.bb1 = apply(bb1, 1, function(x) median(x))
pos.up.bb1 = apply(bb1, 1, function(x) quantile(x, 0.975))
df.betap =c()
df.betap = c(pos.low.bb1)
df.betap = cbind(df.betap, c(pos.median.bb1))
df.betap = cbind(df.betap, c(pos.up.bb1))
df.betap = data.frame(df.betap)
colnames(df.betap) <- c('2.5%', 'Pos.Mean', '97.5%')
import.bb1 = which(df.betap$`97.5%` <0 | df.betap$`2.5%`>0)
import.bb1.p = which((df.betap$`97.5%` <0 | df.betap$`2.5%`>0) & df.betap$Pos.Mean>0)
import.bb1.n = which((df.betap$`97.5%` <0 | df.betap$`2.5%`>0) & df.betap$Pos.Mean<0)
df.betap = cbind(df.betap, OTU = 1:Jsum)
df.betap$import = (df.betap$`97.5%` <0 | df.betap$`2.5%`>0)
df.betap$import <- as.numeric(df.betap$import)

df.betap.B <- df.betap[which(df.betap$OTU<=J[1]), ]
ttt = df.betap.B$import
ttt[ttt==0] = 0.4; ttt[ttt==1] = 1.1
ggplot(df.betap.B, aes(x = OTU, y = Pos.Mean)) + 
  geom_hline(yintercept=0) + theme_bw() + 
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`, color = factor(import)), size=ttt, width = 0.01,
                position=position_dodge(0.05))+
  geom_point() + 
  scale_color_manual("import", breaks=c(0, 1), values=c("grey", "red")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(legend.position = "none") + 
  labs(x= 'bOTU', y = 'Compare pre-treatment to healthy', title = 'Bacterial')+
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=25),
        text=element_text(size=25))

######################### plot 10(b) ############################

df.betap.V <- df.betap[which(df.betap$OTU>J[1]), ]
ttt = df.betap.V$import
ttt[ttt==0] = 0.4; ttt[ttt==1] = 1.1
ggplot(df.betap.V, aes(x = OTU-J[1], y = Pos.Mean)) + 
  geom_hline(yintercept=0) + theme_bw() + 
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`, color = factor(import)), size=ttt, width = 0.01,
                position=position_dodge(0.05))+
  geom_point() + 
  scale_color_manual("import", breaks=c(0, 1), values=c("grey", "red")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(legend.position = "none") + 
  labs(x= 'vOTU', y = 'Compare pre-treatment to healthy', title = 'Viral')+
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=25),
        text=element_text(size=25))

######################### plot 10(c) ############################

bb2 = matrix(NA, Jsum, nsamp)
for(i in 1:nsamp){
  bb2[, i] = beta.st[,2,i] - beta.st[,3,i]
}
pos.low.bb2 = apply(bb2, 1, function(x) quantile(x, 0.025))
pos.median.bb2 = apply(bb2, 1, function(x) median(x))
pos.up.bb2 = apply(bb2, 1, function(x) quantile(x, 0.975))
df.betap2 =c()
df.betap2 = c(pos.low.bb2)
df.betap2 = cbind(df.betap2, c(pos.median.bb2))
df.betap2 = cbind(df.betap2, c(pos.up.bb2))
df.betap2 = data.frame(df.betap2)
colnames(df.betap2) <- c('2.5%', 'Pos.Mean', '97.5%')
import.bb2 = which(df.betap2$`97.5%` <0 | df.betap2$`2.5%`>0)
import.bb2.p = which((df.betap2$`97.5%` <0 | df.betap2$`2.5%`>0) & df.betap2$Pos.Mean>0)
import.bb2.n = which((df.betap2$`97.5%` <0 | df.betap2$`2.5%`>0) & df.betap2$Pos.Mean<0)

df.betap2 = cbind(df.betap2, OTU = 1:Jsum)
df.betap2$import = (df.betap2$`97.5%` <0 | df.betap2$`2.5%`>0)
df.betap2$import <- as.numeric(df.betap2$import)

df.betap2.B <- df.betap2[which(df.betap2$OTU<=J[1]), ]
ttt = df.betap2.B$import; ttt[ttt==0] = 0.4; ttt[ttt==1] = 1.1
ggplot(df.betap2.B, aes(x = OTU, y = Pos.Mean)) + 
  geom_hline(yintercept=0) + theme_bw() + 
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`, color = factor(import)), size=ttt, width = 0.01,
                position=position_dodge(0.05))+
  geom_point() + 
  scale_color_manual("import", breaks=c(0, 1), values=c("grey", "red")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(legend.position = "none") + 
  labs(x= 'bOTU', y = 'Compare post-treatment to healthy', title = 'Bacterial')+
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=25),
        text=element_text(size=25))

######################### plot 10(d) ############################

df.betap2.V <- df.betap2[which(df.betap2$OTU>J[1]), ]
ttt = df.betap2.V$import; ttt[ttt==0] = 0.4; ttt[ttt==1] = 1.1
ggplot(df.betap2.V, aes(x = OTU-J[1], y = Pos.Mean)) + 
  geom_hline(yintercept=0) + theme_bw() + 
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`, color = factor(import)), size=ttt, width = 0.01,
                position=position_dodge(0.05))+
  geom_point() + 
  scale_color_manual("import", breaks=c(0, 1), values=c("grey", "red")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(legend.position = "none") + 
  labs(x= 'vOTU', y = 'Compare post-treatment to healthy', title = 'Viral')+
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=25),
        text=element_text(size=25))

######################### plot 10(e) ############################

bb3 = matrix(NA, Jsum, nsamp)
for(i in 1:nsamp){
  bb3[, i] = beta.st[,2,i] - beta.st[,1,i]
}
pos.low.bb3 = apply(bb3, 1, function(x) quantile(x, 0.025))
pos.median.bb3 = apply(bb3, 1, function(x) median(x))
pos.up.bb3 = apply(bb3, 1, function(x) quantile(x, 0.975))
df.betap3 =c()
df.betap3 = c(pos.low.bb3)
df.betap3 = cbind(df.betap3, c(pos.median.bb3))
df.betap3 = cbind(df.betap3, c(pos.up.bb3))
df.betap3 = data.frame(df.betap3)
colnames(df.betap3) <- c('2.5%', 'Pos.Mean', '97.5%')
import.bb3 = which(df.betap3$`97.5%` <0 | df.betap3$`2.5%`>0)
df.betap3 = cbind(df.betap3, OTU = 1:Jsum)
df.betap3$import = (df.betap3$`97.5%` <0 | df.betap3$`2.5%`>0)
df.betap3$import <- as.numeric(df.betap3$import)

df.betap3.B <- df.betap3[which(df.betap3$OTU<=J[1]), ]
ttt = df.betap3.B$import; ttt[ttt==0] = 0.4; ttt[ttt==1] = 1.1
ggplot(df.betap3.B, aes(x = OTU, y = Pos.Mean)) + 
  geom_point() + 
  geom_hline(yintercept=0) + theme_bw() + 
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`, color = factor(import)), size=ttt, width = 0.01,
                position=position_dodge(0.05))+
  scale_color_manual("import", breaks=c(0, 1), values=c("grey", "red")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(legend.position = "none") + 
  labs(x= 'bOTU', y = 'Compare post-treatment to pre-treatment', title = 'Bacterial')+
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=25),
        text=element_text(size=25))

######################### plot 10(f) ############################

df.betap3.V <- df.betap3[which(df.betap3$OTU>J[1]), ]
ttt = df.betap3.V$import; ttt[ttt==0] = 0.4; ttt[ttt==1] = 1.1
ggplot(df.betap3.V, aes(x = OTU-J[1], y = Pos.Mean)) + 
  geom_point() + 
  geom_hline(yintercept=0) + theme_bw() + 
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`, color = factor(import)), size=ttt, width = 0.01,
                position=position_dodge(0.05))+
  scale_color_manual("import", breaks=c(0, 1), values=c("grey", "red")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(legend.position = "none") + 
  labs(x= 'vOTU', y = 'Compare post-treatment to pre-treatment', title = 'Viral')+
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=25),
        text=element_text(size=25))
