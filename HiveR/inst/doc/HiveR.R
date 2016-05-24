## ----getVersion, echo = FALSE--------------------------------------------
desc <- packageDescription("HiveR")
vers <- paste("version", desc$Version)

## ----SetUp, echo = FALSE, results = "hide"-------------------------------
set.seed(123)
suppressMessages(library(grid))
suppressMessages(library(FuncMap))
suppressMessages(library(HiveR))
suppressMessages(library(sna))
suppressMessages(library(xtable))
suppressMessages(library(bipartite))
suppressMessages(library(reshape))
suppressMessages(library(lattice))
suppressMessages(library(ggplot2))
suppressMessages(library(mvbutils))
suppressMessages(library(rgl))
#suppressMessages(library(animation))
suppressMessages(library(knitr))

if (!file.exists("graphics")) dir.create("graphics")

# Stuff specifically for knitr:

opts_chunk$set(out.width = "0.9\\textwidth", fig.align = "center", fig.width = 7, fig.height = 7, cache = FALSE, echo = FALSE, crop = TRUE)
knit_hooks$set(rgl.static = hook_rgl) # use to capture a single rgl graphic
knit_hooks$set(rgl.mov = hook_plot_custom) # use to capture a movie
if (Sys.which("pdfcrop") != "") knit_hooks$set(crop = hook_pdfcrop) # use pdfcrop if it exists


# Note: defaults are eval = TRUE, echo = TRUE


## ----PPNdata-------------------------------------------------------------
data(Safariland)

## ----PPNA, fig.cap = "Safariland data set using visweb"------------------
visweb(Safariland)

## ----PPN4, fig.cap = "Safariland data set using plotweb"-----------------
plotweb(Safariland)

## ----PPN5, fig.cap = "Safariland data set using gplot (mode = circle)", warning = FALSE----
gplot(Safariland, gmode = "graph", edge.lwd = 0.05,
	vertex.col = c(rep("green", 9), rep("red", 27)),
	mode = "circle")

## ----PPN6, fig.cap = "Safariland data set using gplot (mode = Fruchterman-Reingold)", warning = FALSE----
gplot(Safariland, gmode = "graph", edge.lwd = 0.05,
	vertex.col = c(rep("green", 9), rep("red", 27)))

## ----PPN2, fig.cap = "Hive Panel comparing Safari with Arroyo", eval = FALSE----
#  data(Safari)
#  Safari$nodes$size <- 0.5
#  data(Arroyo)
#  Arroyo$nodes$size <- 0.5
#  
#  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#  #
#  grid.newpage()
#  pushViewport(viewport(layout = grid.layout(2, 1)))
#  #
#  pushViewport(vplayout(1, 1)) # upper plot
#  plotHive(Safari, ch = 0.1, axLabs = c("plants", "pollinators"), axLab.pos = c(0.15, 0.15), rot = c(-90, 90), np = FALSE, axLab.gpar = gpar(fontsize = 16, col = "white"))
#  grid.text("Safari (undisturbed)", x = 0.5, y = 0.95, default.units = "npc", gp = gpar(fontsize = 20, col = "white"))
#  popViewport(2)
#  #
#  pushViewport(vplayout(2, 1)) # lower plot
#  plotHive(Arroyo, ch = 0.1, axLabs = c("plants", "pollinators"), axLab.pos = c(0.15, 0.15), rot = c(-90, 90), np = FALSE, axLab.gpar = gpar(fontsize = 16, col = "white"))
#  grid.text("Arroyo (disturbed)", x = 0.5, y = 0.95, default.units = "npc", gp = gpar(fontsize = 20, col = "white"))
#  

## ----E_coli_1aa, results = "asis"----------------------------------------
tmp <- readLines("network_tf_gene.parsed.dot")[1595:1605]
DOT <- xtable(as.data.frame(tmp))
caption(DOT) <- "Partial contents of .dot file"
label(DOT) <- "DOT"
print(DOT, include.rownames = FALSE, include.colnames = FALSE, hline.after = c(0, nrow(DOT)))

## ----EI, results = "asis"------------------------------------------------
tab <- read.csv(file = "EdgeInst.csv")
EI <- xtable(tab)
caption(EI) <- "Contents of EdgeInst.csv"
label(EI) <- "EI"
print(EI, include.rownames = FALSE)

## ----E_coli_1a, echo = TRUE, tidy = FALSE--------------------------------
EC1 <- dot2HPD(file = "network_tf_gene.parsed.dot",
	node.inst = NULL,
	edge.inst = "EdgeInst.csv",
	desc = "E coli gene regulatory network (RegulonDB)",
	axis.cols = rep("grey", 3))

## ----E_coli_1b, echo = TRUE, size = "footnotesize"-----------------------
sumHPD(EC1)

## ----E_coli_1c, echo = TRUE, size = "footnotesize"-----------------------
EC2 <- mineHPD(EC1, option = "rad <- tot.edge.count")
sumHPD(EC2)

## ----E_coli_1d, echo = TRUE, size = "footnotesize"-----------------------
EC3 <- mineHPD(EC2, option = "axis <- source.man.sink")
sumHPD(EC3)

## ----E_coli_1e, echo = TRUE, size = "footnotesize"-----------------------
EC4 <- mineHPD(EC3, option = "remove zero edge")
sumHPD(EC4)

## ----E_coli_1f, echo = TRUE----------------------------------------------
edges <- EC4$edges
edgesR <- subset(edges, color == 'red')
edgesG <- subset(edges, color == 'green')
edgesO <- subset(edges, color == 'orange')

edges <- rbind(edgesO, edgesG, edgesR)
EC4$edges <- edges

EC4$edges$weight = 0.5


## ----E_coli_2, fig.cap = "Hive panel of E. coli gene regulatory network", out.width = "0.7\\textwidth", fig.width = 2, fig.height = 6----
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 1)))
#
pushViewport(vplayout(1, 1)) # upper plot

plotHive(EC4, dr.nodes = FALSE, ch = 20,
axLabs = c("source", "sink", "manager"),
axLab.pos = c(40, 75, 35),
axLab.gpar = gpar(fontsize = 6, col = "white", lwd = 2),
arrow = c("degree", 150, 100, 180, 70), np = FALSE)
grid.text("native units", x = 0.5, y = 0.05, default.units = "npc", gp = gpar(fontsize = 8, col = "white"))

popViewport(2)
#
pushViewport(vplayout(2, 1)) # middle plot

plotHive(EC4, dr.nodes = FALSE, method = "rank", ch = 100,
#axLabs = c("source", "sink", "manager"),
#axLab.pos = c(100, 125, 180),
#axLab.gpar = gpar(fontsize = 10, col = "white"),
np = FALSE)
grid.text("ranked units", x = 0.5, y = 0.05, default.units = "npc", gp = gpar(fontsize = 8, col = "white"))

popViewport(2)
#
pushViewport(vplayout(3, 1)) # lower plot

plotHive(EC4, dr.nodes = FALSE, method = "norm", ch = 0.1, axLabs = c("source", "sink", "manager"),
axLab.pos = c(0.1, 0.2, 0.2), axLab.gpar = gpar(fontsize = 6, col = "white"), np = FALSE)
grid.text("normed units", x = 0.5, y = 0.05, default.units = "npc", gp = gpar(fontsize = 8, col = "white"))

