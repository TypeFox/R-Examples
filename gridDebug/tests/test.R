
library(grid)
library(gridDebug)

grid.newpage()
grid.rect()
gridTree(grid=TRUE,
         grobNodeAttrs=list(shape="circle",
           fillcolor="black", fontcolor="white"))

grid.newpage()
grid.rect()
gridTreeTips(grid=TRUE,
             grobNodeAttrs=list(shape="circle",
               fillcolor="black", fontcolor="white"))

library(lattice)
xyplot(1:10 ~ 1:10)
grobBrowser()
