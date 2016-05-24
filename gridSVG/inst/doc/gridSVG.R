### R code from vignette source 'gridSVG.Rnw'

###################################################
### code chunk number 1: gridSVG.Rnw:71-74
###################################################
library(grid)
library(gridSVG)



###################################################
### code chunk number 2: gridSVG.Rnw:75-81
###################################################
pushViewport(viewport(gp=gpar(col="black", fill=NA)))
grid.circle(x=0.33, r=unit(2, "inches"), gp=gpar(alpha=0.3, fill="red"))
grid.circle(x=0.67, r=unit(2, "inches"), gp=gpar(alpha=0.3, fill="green"))
popViewport()
gridToSVG()



