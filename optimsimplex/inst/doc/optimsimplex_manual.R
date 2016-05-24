### R code from vignette source 'optimsimplex_manual.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: optimsimplex_manual.Rnw:46-48
###################################################
require(optimsimplex)
options(width=70)


###################################################
### code chunk number 2: optimsimplex_manual.Rnw:125-132
###################################################
coords <- matrix(c(0,1,0,0,0,1),ncol=2)
tmp <- optimsimplex(coords=coords)
s1 <- tmp$newobj
s1
optimsimplex.getallx(s1)
optimsimplex.getn(s1)
optimsimplex.getnbve(s1)


###################################################
### code chunk number 3: optimsimplex_manual.Rnw:143-164
###################################################
rosenbrock <- function(x){
  y <- 100*(x[2]-x[1]^2)^2+(1-x[1])^2
}

mycostf <- function(x, this){
  y <- rosenbrock(x)
  this$nb <- this$nb+1
  return(list(f=y,this=this))
}

mystuff <- list(nb=0)

tmp <- optimsimplex(x0=c(-1.2,1.0), fun=mycostf,method='randbounds',
                    boundsmin=c(-5.0,-5.0), boundsmax=c(5.0,5.0), nbve=5, 
                    data=mystuff)

tmp$newobj

tmp$data

cat(sprintf("Function evaluations: %d\n",tmp$data$nb))


