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
library(circlize)
library(ComplexHeatmap)

######################### plot 9(a) ############################

burn=1:nsamp
col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
mar.pos.cor.array = array(NA, dim=c(Jsum, Jsum, length(burn)))
for(i in burn){
  mar.pos.cor.array[,,i] = cov2cor(tcrossprod(Lambda.st[,,i]) + diag(rep(sig2.st[,i], J)))
}
mar.pos.cor = apply(mar.pos.cor.array, c(1,2), median)
emp.ri = rbind( log(rowMeans(Y1)), log(rowMeans(Y2)))
empirical.cov = cov(log(Y+0.01) - t(emp.ri)[,rep(1:m, J)])
emp.cor = cov2cor(empirical.cov)
poscor.vs.trcor = emp.cor
poscor.vs.trcor[upper.tri(poscor.vs.trcor)] = mar.pos.cor[upper.tri(mar.pos.cor)]
poscor.vs.trcor[upper.tri(poscor.vs.trcor)] = mar.pos.cor[upper.tri(mar.pos.cor)]
poscor.vs.trcor.1 = poscor.vs.trcor[, 1:J[1]]
poscor.vs.trcor.2 = poscor.vs.trcor[, 1:J[2] + J[1]]
ht1 = Heatmap(poscor.vs.trcor.1, column_title = "Bacteria", name = "Cor1", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm"), border = "black"
), column_title_gp = grid::gpar(fontsize = 20), row_names_gp = gpar(fontsize = 1),border_gp = gpar(col = "black"))
ht2 = Heatmap(poscor.vs.trcor.2, column_title = "Viruses",  name = "Cor2", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm") , border = "black"
), show_heatmap_legend = F,column_title_gp = gpar(fontsize = 20), row_names_gp = gpar(fontsize = 1),border_gp = gpar(col = "black"))
ht_list = ht1  + ht2 
draw(ht_list, row_title = "", column_title = "")
decorate_heatmap_body("Cor1", {
  grid.lines(c(0, 1), c(0.342, 0.342), gp = gpar(lty = 1, lwd = 1.5))
})
decorate_heatmap_body("Cor2", {
  grid.lines(c(0, 1), c(0.342, 0.342), gp = gpar(lty = 1, lwd = 1.5))
})

######################### plot 9(b) ############################

realMOFAres <- read.csv("./figures-code/Real Data-MOFA.csv", header=FALSE)
realMOFAres <- as.matrix(realMOFAres)
realMOFAres <- cov2cor(realMOFAres)

poscor.vs.trcor = emp.cor
poscor.vs.trcor[upper.tri(poscor.vs.trcor)] = realMOFAres[upper.tri(realMOFAres)]
poscor.vs.trcor[upper.tri(poscor.vs.trcor)] = realMOFAres[upper.tri(realMOFAres)]
poscor.vs.trcor.1 = poscor.vs.trcor[, 1:J[1]]
poscor.vs.trcor.2 = poscor.vs.trcor[, 1:J[2] + J[1]]
ht1 = Heatmap(poscor.vs.trcor.1, column_title = "Bacteria", name = "Cor1", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm"), border = "black"
), column_title_gp = grid::gpar(fontsize = 20), row_names_gp = gpar(fontsize = 1),border_gp = gpar(col = "black"))
ht2 = Heatmap(poscor.vs.trcor.2, column_title = "Viruses",  name = "Cor2", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm") , border = "black"
), show_heatmap_legend = F,column_title_gp = gpar(fontsize = 20), row_names_gp = gpar(fontsize = 1),border_gp = gpar(col = "black"))
ht_list = ht1  + ht2 
draw(ht_list, row_title = "", column_title = "")
decorate_heatmap_body("Cor1", {
  grid.lines(c(0, 1), c(0.341, 0.341), gp = gpar(lty = 1, lwd = 1.5))
})
decorate_heatmap_body("Cor2", {
  grid.lines(c(0, 1), c(0.341, 0.341), gp = gpar(lty = 1, lwd = 1.5))
})

######################### plot 9(c) ############################

load("./figures-code/Real Data-SPIEC-EASI.RData")
poscor.vs.trcor = emp.cor
rownames(poscor.vs.trcor) = colnames(poscor.vs.trcor) = NULL
poscor.vs.trcor[upper.tri(poscor.vs.trcor)] = se.gl.cor[upper.tri(se.gl.cor)]
poscor.vs.trcor[upper.tri(poscor.vs.trcor)] = se.gl.cor[upper.tri(se.gl.cor)]
poscor.vs.trcor.1 = poscor.vs.trcor[, 1:J[1]]
poscor.vs.trcor.2 = poscor.vs.trcor[, 1:J[2] + J[1]]
ht1 = Heatmap(poscor.vs.trcor.1, column_title = "Bacteria", name = "Cor1", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm"), border = "black"
), column_title_gp = grid::gpar(fontsize = 20), row_names_gp = gpar(fontsize = 1),border_gp = gpar(col = "black"))
ht2 = Heatmap(poscor.vs.trcor.2, column_title = "Viruses",  name = "Cor2", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm") , border = "black"
), show_heatmap_legend = F,column_title_gp = gpar(fontsize = 20), row_names_gp = gpar(fontsize = 1),border_gp = gpar(col = "black"))
ht_list = ht1  + ht2 
draw(ht_list, row_title = "", column_title = "")
decorate_heatmap_body("Cor1", {
  grid.lines(c(0, 1), c(0.341, 0.341), gp = gpar(lty = 1, lwd = 1.5))
})
decorate_heatmap_body("Cor2", {
  grid.lines(c(0, 1), c(0.341, 0.341), gp = gpar(lty = 1, lwd = 1.5))
})
