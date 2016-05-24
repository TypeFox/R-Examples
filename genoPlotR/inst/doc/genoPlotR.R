### R code from vignette source 'genoPlotR.Rnw'

###################################################
### code chunk number 1: genoPlotR.Rnw:37-38
###################################################
library(genoPlotR)


###################################################
### code chunk number 2: genoPlotR.Rnw:45-56
###################################################
data(three_genes)
comparisons[[1]]$col <- apply_color_scheme(c(0.6, 0.4, 0.5), "grey")
names <- c("Huey", "Dewey", "Louie")
names(dna_segs) <- names
tree <- newick2phylog("(((Huey:4.2,Dewey:3.9):3.1,Louie:7.3):1);")
mid_pos <- middle(dna_segs[[1]])
xlims <- list(c(Inf, -Inf), c(-Inf, Inf), c(1850, 2800))
annot <- annotation(x1=c(mid_pos[1], dna_segs[[1]]$end[2]),
                     x2=c(NA, dna_segs[[1]]$end[3]),
                     text=c(dna_segs[[1]]$name[1], "region1"),
                     rot=c(30, 0), col=c("blue", "black"))


###################################################
### code chunk number 3: quick_plot
###################################################
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              annotations=annot, annotation_height=1.3,
              tree=tree, tree_width=2,
              xlims=xlims,
              main="Comparison of Huey, Dewey and Louie")



###################################################
### code chunk number 4: genoPlotR.Rnw:70-71
###################################################
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              annotations=annot, annotation_height=1.3,
              tree=tree, tree_width=2,
              xlims=xlims,
              main="Comparison of Huey, Dewey and Louie")



###################################################
### code chunk number 5: genoPlotR.Rnw:82-83 (eval = FALSE)
###################################################
## help.start()


###################################################
### code chunk number 6: genoPlotR.Rnw:86-87 (eval = FALSE)
###################################################
## library(help=genoPlotR)


###################################################
### code chunk number 7: genoPlotR.Rnw:93-95 (eval = FALSE)
###################################################
## help("read_functions")
## help("plot_gene_map")


###################################################
### code chunk number 8: genoPlotR.Rnw:117-127
###################################################
names1 <- c("feat1", "feat2", "feat3")
starts1 <- c(2, 1000, 1050)
ends1 <- c(600, 800, 1345)
strands1 <- c("-", -1, 1)
cols1 <- c("blue", "grey", "red")
df1 <- data.frame(name=names1, start=starts1, end=ends1,
                  strand=strands1, col=cols1)

dna_seg1 <- dna_seg(df1)
str(dna_seg1)


###################################################
### code chunk number 9: genoPlotR.Rnw:141-148
###################################################
starts1 <- c(2, 1000, 1050)
ends1 <- c(600, 800, 1345)
starts2 <- c(50, 800, 1200)
ends2 <- c(900, 1100, 1322)
comparison1 <- as.comparison(data.frame(start1=starts1, end1=ends1,
                                        start2=starts2, end2=ends2))
str(comparison1)


###################################################
### code chunk number 10: genoPlotR.Rnw:157-160
###################################################
mid_pos <- middle(dna_segs[[1]])
annot1 <- annotation(x1=mid_pos, text=dna_segs[[1]]$name)
str(annot1)


###################################################
### code chunk number 11: genoPlotR.Rnw:166-168
###################################################
tree <- newick2phylog("(((A_aaa:4.2,B_bbb:3.9):3.1,C_ccc:7.3):1);")
str(tree$leaves)


###################################################
### code chunk number 12: loadLib
###################################################
library(genoPlotR)


###################################################
### code chunk number 13: ex1_create_dna_segs
###################################################
df1 <- data.frame(name=c("feat1", "feat2", "feat3"),
                  start=c(2, 1000, 1050),
                  end=c(600, 800, 1345),
                  strand=c(-1, -1, 1),
                  col=c("blue", "grey", "red"))
dna_seg1 <- dna_seg(df1)
df2 <- data.frame(name=c("feat1", "feat2", "feat3"),
                  start=c(50, 800, 1200),
                  end=c(900, 1100, 1322),
                  strand=c(-1, 1, 1),
                  col=c("blue", "grey", "red"))
dna_seg2 <- dna_seg(df2)
df3 <- data.frame(name=c("feat1", "feat2", "feat3"),
                  start=c(1899, 2108, 2803),
                  end=c(2034, 2732, 3620),
                  strand=c(-1, -1, 1),
                  col=rep("blue", 3))
dna_seg3 <- dna_seg(df3)
dna_segs <- list(dna_seg1, dna_seg2, dna_seg3)


###################################################
### code chunk number 14: ex1_create_comps
###################################################
df4 <- data.frame(start1=dna_seg1$start,
                  end1=dna_seg1$end,
                  start2=dna_seg2$start,
                  end2=dna_seg2$end)
comparison1 <- comparison(df4)
df5 <- data.frame(start1=c(50, 800),
                  end1=c(500, 1100),
                  start2=c(1899, 2732),
                  end2=c(2034, 2508),
                  col=c("#67000D", "#08306B"))
comparison2 <- comparison(df5)
comparisons <- list(comparison1, comparison2)


###################################################
### code chunk number 15: ex1_raw_plot
###################################################
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons)


###################################################
### code chunk number 16: genoPlotR.Rnw:320-321
###################################################
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons)


###################################################
### code chunk number 17: ex1_col_scheme
###################################################
comparisons[[1]]$col <- apply_color_scheme(c(0.6, 0.4, 0.5), "grey")


###################################################
### code chunk number 18: ex1_tree
###################################################
names <- c("Huey", "Dewey", "Louie")
names(dna_segs) <- names
tree_HDL <- newick2phylog("(((Huey:4.2,Dewey:3.9):3.1,Louie:7.3):1);")


###################################################
### code chunk number 19: ex1_annots
###################################################
mid_pos <- middle(dna_segs[[1]])
annot <- annotation(x1=c(mid_pos[1], dna_segs[[1]]$end[2]),
                    x2=c(NA, dna_segs[[1]]$end[3]),
                    text=c(dna_segs[[1]]$name[1], "region1"),
                    rot=c(30, 0), col=c("grey", "black"))


###################################################
### code chunk number 20: ex1_adv_plot
###################################################
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              annotations=annot, annotation_height=1.3,
              tree=tree_HDL, tree_width=2,
              main="Comparison of Huey, Dewey and Louie")


###################################################
### code chunk number 21: genoPlotR.Rnw:380-381
###################################################
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              annotations=annot, annotation_height=1.3,
              tree=tree_HDL, tree_width=2,
              main="Comparison of Huey, Dewey and Louie")


###################################################
### code chunk number 22: ex_2_load_fake (eval = FALSE)
###################################################
## BH <- try(read_dna_seg_from_file("NC_005956.gbk"))
## BQ <- try(read_dna_seg_from_file("NC_005955.gbk"))
## BH_vs_BQ <- try(read_comparison_from_blast("NC_005956_vs_NC_005955.blast"))


###################################################
### code chunk number 23: ex2_load_silent
###################################################
data(barto)
BH <- barto$dna_segs[[3]]
BQ <- barto$dna_segs[[4]]
BH_vs_BQ <- barto$comparisons[[3]]


###################################################
### code chunk number 24: ex2_plot
###################################################
xlims <- list(c(1,50000), c(1,50000))
plot_gene_map(dna_segs=list(BH, BQ),
              comparisons=list(BH_vs_BQ),
              xlims=xlims,
              main="BH vs BQ, comparison of the first 50 kb",
              gene_type="side_blocks",
              dna_seg_scale=TRUE, scale=FALSE)



###################################################
### code chunk number 25: genoPlotR.Rnw:472-473
###################################################
xlims <- list(c(1,50000), c(1,50000))
plot_gene_map(dna_segs=list(BH, BQ),
              comparisons=list(BH_vs_BQ),
              xlims=xlims,
              main="BH vs BQ, comparison of the first 50 kb",
              gene_type="side_blocks",
              dna_seg_scale=TRUE, scale=FALSE)



###################################################
### code chunk number 26: ex3_read_data
###################################################
bbone_file <- system.file('extdata/barto.backbone', package = 'genoPlotR')
bbone <- read_mauve_backbone(bbone_file, ref=2, filter_low=10000)
names <- c("B_bacilliformis", "B_grahamii", "B_henselae", "B_quintana")
names(bbone$dna_segs) <- names


###################################################
### code chunk number 27: ex3_length
###################################################
for (i in 1:length(bbone$comparisons)){
  cmp <- bbone$comparisons[[i]]
  bbone$comparisons[[i]]$length <- 
    abs(cmp$end1 - cmp$start1) + abs(cmp$end2 - cmp$start2)
}


###################################################
### code chunk number 28: ex3_plot
###################################################
plot_gene_map(dna_segs=bbone$dna_segs, 
              comparisons=bbone$comparisons,
              global_color_scheme=c("length", "increasing", "red_blue", 0.7),
              override_color_schemes=TRUE)


###################################################
### code chunk number 29: genoPlotR.Rnw:524-525
###################################################
plot_gene_map(dna_segs=bbone$dna_segs, 
              comparisons=bbone$comparisons,
              global_color_scheme=c("length", "increasing", "red_blue", 0.7),
              override_color_schemes=TRUE)


###################################################
### code chunk number 30: ex4_load_tree
###################################################
data(barto)
tree_barto <- newick2phylog("(BB:2.5,(BG:1.8,(BH:1,BQ:0.8):1.9):3);")


###################################################
### code chunk number 31: ex4_xlims
###################################################
xlims <- list(c(1445000, 1415000, 1380000, 1412000),
              c(  10000,   45000,   50000,   83000, 90000, 120000),
              c(  15000,   36000,   90000,  120000, 74000,  98000),
              c(   5000,   82000))


###################################################
### code chunk number 32: ex4_annots
###################################################
annots <- lapply(barto$dna_segs, function(x){
  mid <- middle(x)
  annot <- annotation(x1=mid, text=x$name, rot=30)
  idx <- grep("^[^B]", annot$text, perl=TRUE)
  annot[idx[idx %% 4 == 0],] 
})


###################################################
### code chunk number 33: ex4_plot
###################################################
plot_gene_map(barto$dna_segs, barto$comparisons, tree=tree_barto,
              annotations=annots,
              xlims=xlims,
              limit_to_longest_dna_seg=FALSE,
              dna_seg_scale=TRUE, scale=FALSE,
              main="Comparison of homologous segments in 4 Bartonella genomes")


###################################################
### code chunk number 34: genoPlotR.Rnw:589-590
###################################################
plot_gene_map(barto$dna_segs, barto$comparisons, tree=tree_barto,
              annotations=annots,
              xlims=xlims,
              limit_to_longest_dna_seg=FALSE,
              dna_seg_scale=TRUE, scale=FALSE,
              main="Comparison of homologous segments in 4 Bartonella genomes")


###################################################
### code chunk number 35: ex5_load_annot
###################################################
data(chrY_subseg)
genes_homo <- unique(chrY_subseg$dna_segs[[1]]$gene)
x_homo <- sapply(genes_homo, function(x)
                 range(chrY_subseg$dna_segs[[1]]
                       [chrY_subseg$dna_segs[[1]]$gene == x,])
                 )
annot_homo <- annotation(x1=x_homo[1,], x2=x_homo[2,],
                         text=dimnames(x_homo)[[2]])
genes_pan <- unique(chrY_subseg$dna_segs[[2]]$gene)
x_pan <- sapply(genes_pan, function(x)
                range(chrY_subseg$dna_segs[[2]]
                      [chrY_subseg$dna_segs[[2]]$gene == x,])
                 )
annot_pan <- annotation(x1=x_pan[1,], x2=x_pan[2,],
                        text=dimnames(x_pan)[[2]])


###################################################
### code chunk number 36: ex5_plot
###################################################
main <- "Comparison of two subsegments in H. sapiens and P. troglodytes"
plot_gene_map(chrY_subseg$dna_segs, chrY_subseg$comparison,
              annotations=list(annot_homo, annot_pan),
              dna_seg_scale=TRUE,
              main=main,
              scale=FALSE)


###################################################
### code chunk number 37: genoPlotR.Rnw:632-633
###################################################
main <- "Comparison of two subsegments in H. sapiens and P. troglodytes"
plot_gene_map(chrY_subseg$dna_segs, chrY_subseg$comparison,
              annotations=list(annot_homo, annot_pan),
              dna_seg_scale=TRUE,
              main=main,
              scale=FALSE)


###################################################
### code chunk number 38: ex6_startdevice
###################################################
sysname <- Sys.info()["sysname"]
if (sysname == "Windows"){
  pdf(file="ex6.pdf", width=6, height=8, onefile=TRUE)
} else {
  cairo_pdf(file="ex6.pdf", width=6, height=8, onefile=TRUE)
}


###################################################
### code chunk number 39: ex6_vpext
###################################################
pushViewport(viewport(layout=grid.layout(3,1,
                        heights=unit(c(1,1.3,0.8), rep("null", 3))),
                      name="overall_vp"))


###################################################
### code chunk number 40: ex6_plot
###################################################
## Panel A
pushViewport(viewport(layout.pos.row=1, name="panelA"))
plot_gene_map(dna_segs=bbone$dna_segs, comparisons=bbone$comparisons,
              dna_seg_scale=c(FALSE, FALSE, FALSE, TRUE),
              scale=FALSE, main="A", main_pos="left", plot_new=FALSE)
upViewport()
## Panel B
pushViewport(viewport(layout.pos.row=2, name="panelB"))
plot_gene_map(barto$dna_segs, barto$comparisons, 
              annotations=annots, tree=tree_barto, xlims=xlims,
              limit_to_longest_dna_seg=FALSE, scale=FALSE,
              dna_seg_scale=TRUE, main="B", main_pos="left",
              annotation_height=0.6, annotation_cex=0.5, 
              plot_new=FALSE)
upViewport()
## Panel C
pushViewport(viewport(layout.pos.row=3, name="panelC"))
plot_gene_map(chrY_subseg$dna_segs, chrY_subseg$comparison,
              annotations=list(annot_homo, annot_pan),
              dna_seg_scale=TRUE, scale=FALSE, main="C", main_pos="left",
              plot_new=FALSE)
upViewport(0)


###################################################
### code chunk number 41: ex6_currentvp
###################################################
grid_list <- grid.ls(grob=TRUE, viewports=TRUE, print=FALSE)
str(grid_list)
current.vpTree()


###################################################
### code chunk number 42: ex6_mod_panela
###################################################
downViewport("panelA")
for (i in 1:length(names)){
  new_label <- sub("_", ". ", names[[i]])
  grid.edit(paste("label", i, sep="."), label=new_label, 
            gp=gpar(fontface="italic"))
}
grid.remove("label.2")
upViewport(0)


###################################################
### code chunk number 43: ex6_mod_panelb
###################################################
downViewport("panelB")
downViewport("dna_seg.3.2")
grid.rect(height = unit(2.2, "npc"), gp=gpar(col="red", lwd=2, fill=0))
upViewport(0)


###################################################
### code chunk number 44: genoPlotR.Rnw:734-735
###################################################
dev.off()


###################################################
### code chunk number 45: ex7_starGrob
###################################################
starGrob <- function(gene, ...){
  ## Coordinates for the star
  x <- sin(((0:5)/2.5)*pi)*(gene$end-gene$start)/2 + (gene$end+gene$start)/2
  y <- cos(((0:5)/2.5)*pi)*gene$strand*0.5 + 0.5
  idx <- c(1, 3, 5, 2, 4, 1)
  ## Attribute line_col only if present in the gene
  line_col <- if (!is.null(gene$line_col)) gene$line_col else gene$col
  ## Having a conditional transparency, depending on a length cut-off
  ## passed via dots
  length_cutoff <- list(...)$length_cutoff
  if (!is.null(length_cutoff)){
    alpha <- if ((gene$end-gene$start) < length_cutoff)  0.3 else  0.8
  } else alpha <- 1
  
  ## Grobs
  g <- polygonGrob(x[idx], y[idx], gp=gpar(fill=gene$col, col=line_col,
                                     lty=gene$lty, lwd=gene$lwd, alpha=alpha),
                   default.units="native")
  t <- textGrob(label="***", x=(gene$end+gene$start)/2, y=0.5,
                default.units="native")
  gList(g, t)
}


###################################################
### code chunk number 46: ex7_replace_geneType
###################################################
barto$dna_segs[[2]]$gene_type <- "starGrob"
barto$dna_segs[[4]]$gene_type <- "starGrob"
line_col <- rep(1:20, (nrow(barto$dna_segs[[3]]) %% 20)+1)
barto$dna_segs[[2]]$line_col <- line_col[1:nrow(barto$dna_segs[[2]])]
plot_gene_map(barto$dna_segs, barto$comparisons, tree = tree_barto,
              annotations = annots,
              xlims = xlims, 
              dna_seg_scale = TRUE,
              length_cutoff = 600)


###################################################
### code chunk number 47: genoPlotR.Rnw:793-794
###################################################
barto$dna_segs[[2]]$gene_type <- "starGrob"
barto$dna_segs[[4]]$gene_type <- "starGrob"
line_col <- rep(1:20, (nrow(barto$dna_segs[[3]]) %% 20)+1)
barto$dna_segs[[2]]$line_col <- line_col[1:nrow(barto$dna_segs[[2]])]
plot_gene_map(barto$dna_segs, barto$comparisons, tree = tree_barto,
              annotations = annots,
              xlims = xlims, 
              dna_seg_scale = TRUE,
              length_cutoff = 600)


###################################################
### code chunk number 48: ex7_definesegplots
###################################################
seg_plots <- lapply(barto$dna_segs, function(ds){
  x <- seq(1, range(ds)[2], by=1000)
  y <- jitter(seq(100, 300, length=length(x)), amount=50)
  seg_plot(func=linesGrob, args=list(x=x, y=y, gp=gpar(col=grey(0.3), lty=2)))
})
str(seg_plots[[1]])


###################################################
### code chunk number 49: ex7_plotwithsegplots
###################################################
plot_gene_map(barto$dna_segs, barto$comparisons, tree = tree_barto,
              annotations = annots,
              xlims = xlims,
              seg_plots = seg_plots,
              seg_plot_height = 0.5,
              seg_plot_height_unit = "null",
              seg_plot_yaxis = 2,
              seg_plot_yaxis_cex = 0.7)


###################################################
### code chunk number 50: genoPlotR.Rnw:833-834
###################################################
plot_gene_map(barto$dna_segs, barto$comparisons, tree = tree_barto,
              annotations = annots,
              xlims = xlims,
              seg_plots = seg_plots,
              seg_plot_height = 0.5,
              seg_plot_height_unit = "null",
              seg_plot_yaxis = 2,
              seg_plot_yaxis_cex = 0.7)


