## ----2d_data-------------------------------------------------------------
#for the dataset
data('mtcars')
mcor <-cor(mtcars)
# Print mcor and round to 2 digits
round(mcor,digits=2)

## ----2d_xtable, results='asis'-------------------------------------------
library(xtable)
print(xtable(mcor), type='html', comment=F)

## ----2d_plot, dev='png', warning=FALSE-----------------------------------
library(corrplot)
corrplot(mcor)

## ----2d_network_data, dev='png'------------------------------------------
library(igraph)
# Specify edges for a directed graph
gd <-graph(c(1,2, 2,3, 2,4, 1,4, 5,5, 3,6))
plot(gd)
# For an undirected graph
gu <-graph(c(1,2, 2,3, 2,4, 1,4, 5,5, 3,6),directed=FALSE)
# No labels
plot(gu,vertex.label=NA)

