# A. Data preparation
setwd("C:/aErer"); library(diagram)
M <- matrix(nrow = 4, ncol = 4, data = 0)
M[2, 1] <- 1; M[4, 2] <- 2; M[3, 4] <- 3; M[1, 3] <- 4; M

# B. Display the diagram on a screen device
windows(width = 5.3, height = 2.5, family = "serif")
par(mai = c(0, 0, 0, 0))
book <- plotmat(A = M, pos = c(1, 2, 1), curve = 0.3, 
  name = c("Empirical Study", "Proposal", "Manuscript", "R Program"),
  lwd = 1, box.lwd = 1.5, cex.txt = 0.8, arr.type = "triangle",
  box.size = c(0.15, 0.1, 0.1, 0.1), box.cex = 0.75, 
  box.type = c("hexa", "ellipse", "ellipse", "ellipse"),
  box.prop = 0.4, box.col = c("pink", "yellow", "green", "orange"),
  lcol = "purple", arr.col = "purple")
names(book); book[["rect"]]
showDiagram <- recordPlot()
                         
# C. Save the graph on a file device
pdf(file = "fig_showDiagram.pdf", width = 5.3, height=2.5, family="serif")
replayPlot(showDiagram); dev.off()   