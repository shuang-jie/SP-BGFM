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
col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))

################### Plot Fig 6(a)  ############################

burn=1:nsamp
mar.pos.cor.array = array(NA, dim=c(Jsum, Jsum, length(burn)))
for(i in burn){
  mar.pos.cor.array[,,i] = cov2cor(tcrossprod(Lambda.st[,,i]) + diag(rep(sig2.st[,i], J)))
}
mar.pos.cor = apply(mar.pos.cor.array, c(1,2), median)
Hplot = Heatmap(mar.pos.cor[J.ls[[1]], J.ls[[1]]], name = "Cor", col = col_fun, cluster_rows = T,  cluster_columns = T, heatmap_legend_param = list(
  at = c(-1, -0.5, 0, 0.5, 1),
  labels =c("-1", "-0.5", "0", "0.5","1"),
  title = "Cor",
  legend_height = unit(15, "cm")
), column_names_gp = gpar(fontsize = 8.5), row_names_gp = gpar(fontsize = 5))
new.order.1 = row_order(Hplot)
Hplot = Heatmap(mar.pos.cor[J.ls[[2]], J.ls[[2]]], name = "Cor", col = col_fun, cluster_rows = T,  cluster_columns = T, heatmap_legend_param = list(
  at = c(-1, -0.5, 0, 0.5, 1),
  labels =c("-1", "-0.5", "0", "0.5","1"),
  title = "Cor",
  legend_height = unit(15, "cm")
), column_names_gp = gpar(fontsize = 8.5), row_names_gp = gpar(fontsize = 5))
new.order.2 = row_order(Hplot)
new.order = c(new.order.1, new.order.2+J[1])
poscor.vs.trcor = true.cor[new.order, new.order]
poscor.vs.trcor[upper.tri(poscor.vs.trcor)] = (mar.pos.cor[new.order, new.order])[upper.tri(mar.pos.cor)]
poscor.vs.trcor.1 = poscor.vs.trcor[, 1:J[1]]
poscor.vs.trcor.2 = poscor.vs.trcor[, 1:J[2] + J[1]]
ht1 = Heatmap(poscor.vs.trcor.1, column_title = "Group 1", name = "Cor1", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm"), border = "black"
), border_gp = gpar(col = "black"), column_title_gp=gpar(fontsize=20))
ht2 = Heatmap(poscor.vs.trcor.2, column_title = "Group 2",  name = "Cor2", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm") , border = "black"
), show_heatmap_legend = FALSE, border_gp = gpar(col = "black"), column_title_gp=gpar(fontsize=20))
ht_list = ht1  + ht2 
draw(ht_list, row_title = "", column_title = "")
decorate_heatmap_body("Cor1",{
  grid.lines(c(0, 1), c(0.25, 0.25), gp = gpar(lty = 1, lwd = 1.5))
})
decorate_heatmap_body("Cor2",{
  grid.lines(c(0, 1), c(0.25, 0.25), gp = gpar(lty = 1, lwd = 1.5))
})

################### Plot Fig 6(b)  ############################

est.cor <- read.csv("./figures-code/Sim 2-MOFA.csv", header=FALSE)
est.cor <- as.matrix(est.cor); est.cor = cov2cor(est.cor)
est.cor <- est.cor[new.order, new.order]
poscor.vs.trcor.MOFA = true.cor[new.order, new.order]
poscor.vs.trcor.MOFA[upper.tri(poscor.vs.trcor.MOFA)] = est.cor[upper.tri(est.cor)]
poscor.vs.trcor.MOFA[upper.tri(poscor.vs.trcor.MOFA)] = est.cor[upper.tri(est.cor)]
poscor.vs.trcor.MOFA.1 = poscor.vs.trcor.MOFA[, 1:J[1]]
poscor.vs.trcor.MOFA.2 = poscor.vs.trcor.MOFA[, 1:J[2] + J[1]]
ht1 = Heatmap(poscor.vs.trcor.MOFA.1, column_title = "Group 1", name = "Cor1", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm"), border = "black"
), border_gp = gpar(col = "black"), column_title_gp=gpar(fontsize=20))
ht2 = Heatmap(poscor.vs.trcor.MOFA.2, column_title = "Group 2",  name = "Cor2", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm") , border = "black"
),show_heatmap_legend = FALSE, border_gp = gpar(col = "black"), column_title_gp=gpar(fontsize=20))
ht_list = ht1  + ht2 
draw(ht_list, row_title = "", column_title = "")
decorate_heatmap_body("Cor1",{
  grid.lines(c(0, 1), c(0.25, 0.25), gp = gpar(lty = 1, lwd = 1.5))
})
decorate_heatmap_body("Cor2",{
  grid.lines(c(0, 1), c(0.25, 0.25), gp = gpar(lty = 1, lwd = 1.5))
})

################### Plot Fig 6(c)  ############################

load("./figures-code/Sim 2-SPIEC-EASI.RData")
se.gl.cor = as.matrix(se.gl.cor)
se.gl.cor = se.gl.cor[new.order, new.order]
poscor.vs.trcor.SPIEC.EASI = true.cor[new.order, new.order]
poscor.vs.trcor.SPIEC.EASI[upper.tri(poscor.vs.trcor.SPIEC.EASI)] = se.gl.cor[upper.tri(se.gl.cor)]
poscor.vs.trcor.SPIEC.EASI[upper.tri(poscor.vs.trcor.SPIEC.EASI)] = se.gl.cor[upper.tri(se.gl.cor)]
poscor.vs.trcor.SPIEC.EASI.1 = poscor.vs.trcor.SPIEC.EASI[, 1:J[1]]
poscor.vs.trcor.SPIEC.EASI.2 = poscor.vs.trcor.SPIEC.EASI[, 1:J[2] + J[1]]
ht1 = Heatmap(poscor.vs.trcor.SPIEC.EASI.1, column_title = "Group 1", name = "Cor1", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm"), border = "black"
), border_gp = gpar(col = "black"), column_title_gp=gpar(fontsize=20))
ht2 = Heatmap(poscor.vs.trcor.SPIEC.EASI.2, column_title = "Group 2",  name = "Cor2", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm") , border = "black"
),show_heatmap_legend = FALSE, border_gp = gpar(col = "black"), column_title_gp=gpar(fontsize=20))
ht_list = ht1  + ht2 
draw(ht_list, row_title = "", column_title = "")
decorate_heatmap_body("Cor1",{
  grid.lines(c(0, 1), c(0.25, 0.25), gp = gpar(lty = 1, lwd = 1.5))
})
decorate_heatmap_body("Cor2",{
  grid.lines(c(0, 1), c(0.25, 0.25), gp = gpar(lty = 1, lwd = 1.5))
})
