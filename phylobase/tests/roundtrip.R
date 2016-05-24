library(phylobase)
library(ape)

## set.seed(1)
## t0A <- rcoal(5)
t0 <- read.tree(text="((t4:0.3210275554,(t2:0.2724586465,t3:0.2724586465):0.0485689089):0.1397952619,(t5:0.07551818331,t1:0.07551818331):0.385304634);")
## hack around variability in ape:
##   read.tree() and rcoal() produce sets of
##     elements in different orders
t0 <- unclass(t0)[c("edge","edge.length","tip.label","Nnode")]
class(t0) <- "phylo"

## phylo -> phylo4 -> phylo
t1 <- as(t0,"phylo4")
t5 <- as(t1,"phylo")
stopifnot(identical(t0,t5))

## phylo4 -> phylo4vcov -> phylo4 -> phylo
t2<-as(t1,"phylo4vcov")
t3<-as(t2,"phylo4")
t4<-as(t3,"phylo")
stopifnot(identical(t4$edge,t0$edge) &&
          identical(t4$tip.label,t0$tip.label) &&
          identical(t4$Nnode,t0$Nnode) &&
          max(abs(t4$edge.length-t0$edge.length))<1e-10)

## UNROOTED
t6 <- ape::unroot(t0)
## hack around ape conversion issues:
##  unroot() converts integer to double
storage.mode(t6$edge) <- "integer"
storage.mode(t6$Nnode) <- "integer"
t7 <- as(as(t6,"phylo4"),"phylo")
stopifnot(identical(t6,t7))


## EXPLICIT ROOT EDGE
t8 <- t0
t8$root.edge <- 0.5
t9 <- as(as(t8,"phylo4"),"phylo")
stopifnot(identical(t8,t9))
