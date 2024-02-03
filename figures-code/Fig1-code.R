rm(list = ls())

######################### load real data #########################

load("./real-data/Filtered7539OTUs.RData")

######################### load package #########################

library(circlize)
library(ComplexHeatmap)
library(MetaLonDA)

######################### plot Fig 1(a) #########################

norm.Y1 = MetaLonDA::normalize(Y1); norm.Y2 = MetaLonDA::normalize(Y2)
normalized.count = cbind(norm.Y1, norm.Y2)
pre_trt = seq(1, n, 3)
post_trt = seq(2, n, 3)
normal = seq(3, n, 3)
normalized.count = log(normalized.count[c(pre_trt, post_trt, normal), ]+0.01)

col_fun = colorRamp2(c(-4, 0, 15), c("green", "white", "red"))
colnames(normalized.count) <- c(paste0('b', 1:J[1]), paste0('v', 1:J[2]))
normalized.count = t(normalized.count)
dim(normalized.count)
rownames(normalized.count) <- NULL
colnames(normalized.count) <- NULL
normalized.count.pre = normalized.count[, 1:20]
normalized.count.post = normalized.count[, 1:20+20]
normalized.count.norm = normalized.count[, 1:20+40]
ht1 = Heatmap(normalized.count.pre, column_title = "Pre-trt Samples", name = "mat1", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = " ",
  legend_height = unit(15, "cm"), border = "black", grid_width = unit(1, "cm")
), column_title_gp = gpar(fontsize = 20), row_names_gp = gpar(fontsize = 1),border_gp = gpar(col = "black"))
ht2 = Heatmap(normalized.count.post, column_title = "Post-trt Samples",  name = "mat2", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = " ",
  legend_height = unit(15, "cm") , border = "black"
), show_heatmap_legend = FALSE, column_title_gp = gpar(fontsize = 20), row_names_gp = gpar(fontsize = 1),border_gp = gpar(col = "black"))
ht3 = Heatmap(normalized.count.norm, column_title = "Healthy Samples",  name = "mat3", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = " ",
  legend_height = unit(15, "cm") , border = "black"
), show_heatmap_legend = FALSE,  column_title_gp = gpar(fontsize = 20), row_names_gp = gpar(fontsize = 1),border_gp = gpar(col = "black"))

ht_list = ht1  + ht2 + ht3
draw(ht_list, row_title = "", column_title = "")
decorate_heatmap_body("mat1", {
  grid.lines(c(0.005, .995), c(0.334, 0.334), gp = gpar(lty = 1, lwd = 3))
})
decorate_heatmap_body("mat2", {
  grid.lines(c(0.005, .995), c(0.334, 0.334), gp = gpar(lty = 1, lwd = 3))
})
decorate_heatmap_body("mat3", {
  grid.lines(c(0.005, .995), c(0.334, 0.334), gp = gpar(lty = 1, lwd = 3))
  grid.text("bOTU", 1.11, 0.9, default.units = "npc", gp = gpar(fontsize = 19))
  grid.text("vOTU", 1.11, 0.09, default.units = "npc", gp = gpar(fontsize = 19))
})

######################### plot Fig 1(b) #########################

normalized.count = cbind(norm.Y1, norm.Y2)
normalized.count = log(normalized.count[c(pre_trt, post_trt, normal), ]+0.01)
emp.cor = cor(normalized.count)
colnames(emp.cor) = rownames(emp.cor) = NULL
m = 2; M = rep(c(1:m), J)
J.ls = lapply(1:m, function(x) which(M==x))
Hplot = Heatmap(emp.cor[J.ls[[1]], J.ls[[1]]], name = "Cor", col = col_fun, cluster_rows = T,  cluster_columns = T, heatmap_legend_param = list(
  at = c(-1, -0.5, 0, 0.5, 1),
  labels =c("-1", "-0.5", "0", "0.5","1"),
  title = "Cor",
  legend_height = unit(15, "cm")
), column_names_gp = gpar(fontsize = 8.5), row_names_gp = gpar(fontsize = 5))
new.order.1 = row_order(Hplot)
Hplot = Heatmap(emp.cor[J.ls[[2]], J.ls[[2]]], name = "Cor", col = col_fun, cluster_rows = T,  cluster_columns = T, heatmap_legend_param = list(
  at = c(-1, -0.5, 0, 0.5, 1),
  labels =c("-1", "-0.5", "0", "0.5","1"),
  title = "Cor",
  legend_height = unit(15, "cm")
), column_names_gp = gpar(fontsize = 8.5), row_names_gp = gpar(fontsize = 5))
new.order.2 = row_order(Hplot)
new.order = c(new.order.1, new.order.2+J[1])

emp.cor.reorder = emp.cor[new.order, new.order]
col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
emp.cor.reorder.1 = emp.cor.reorder[, 1:J[1]]
emp.cor.reorder.2 = emp.cor.reorder[, 1:J[2] + J[1]]
ht1 = Heatmap(emp.cor.reorder.1, column_title = "Bacteria", name = "Cor1", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm"), border = "black"
), column_title_gp = gpar(fontsize = 20), row_names_gp = gpar(fontsize = 5),border_gp = gpar(col = "black"))
ht2 = Heatmap(emp.cor.reorder.2, column_title = "Viruses",  name = "Cor2", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm") , border = "black"
), show_heatmap_legend = FALSE,   column_title_gp = gpar(fontsize = 20), row_names_gp = gpar(fontsize = 5),border_gp = gpar(col = "black"))
ht_list = ht1  + ht2 
draw(ht_list, row_title = "", column_title = "")
decorate_heatmap_body("Cor1", {
  grid.lines(c(0.00025, 1), c(0.343, 0.343), gp = gpar(lty = 1, lwd = 1.5))
})
decorate_heatmap_body("Cor2", {
  grid.lines(c(0.00025, 1), c(0.343, 0.343), gp = gpar(lty = 1, lwd = 1.5))
})



