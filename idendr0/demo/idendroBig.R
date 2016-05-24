## Demo of dendrogram + heat map integrated with GGobi scatter plot
## matrix using "big data" (5000 observations).
##

library(idendr0) # hca5000, idendro

# load "big data"
data(hca5000)

# visualize clusters
idendro(hca5000$hx, hca5000$x, ggobi = TRUE)
