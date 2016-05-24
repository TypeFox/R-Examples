## ----echo=FALSE,message=F,comment=F,warning=F----------------------------
require(knitr)

## ---- echo=FALSE---------------------------------------------------------
GE = paste('$\\varepsilon_{',rep(1:4,each=5),
           rep(1:5,length.out=20),'}$',sep='')
GE = matrix(GE,5,4)
colnames(GE) = paste('E',1:4,sep='')
rownames(GE) = paste('G',1:5,sep='')
kable(GE)

## ---- echo=FALSE---------------------------------------------------------
GE = paste('$q_{',rep(1:4,each=5),
           rep(1:5,length.out=20),'}$',sep='')
GE = matrix(GE,5,4)
colnames(GE) = paste('E',1:4,sep='')
rownames(GE) = paste('G',1:5,sep='')
kable(GE)

