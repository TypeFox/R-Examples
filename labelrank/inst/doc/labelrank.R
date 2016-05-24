## ---- echo=FALSE, results='asis'-----------------------------------------
load('../data/lr.nom.rda')
load('../data/lr.num.rda')
load('../data/y.rda')
x <- lr.nom[1:5,]
knitr::kable(cbind(x,y[1:5,]),col.names=c('$x_1$','$x_2$','A','B','C'),row.names=F,caption='Label ranking example')

## ---- echo=T,results='markup'--------------------------------------------
library(labelrank)
x <- lr.nom[1:5,]
ex.y <- y[1:5,]
nb_rank(x,ex.y,x[5,])

## ---- echo=T,results='markup'--------------------------------------------

ex.knn <- lr.num[1:5,]
nn_rank(ex.knn,ex.y,ex.knn[5,])

