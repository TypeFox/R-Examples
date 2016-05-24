library(pcalg)
## setwd("/sfs/u/kalischm/research/packages/unifDAGs/")
## source("aux_general.R")
## source("randDAG.R")

### Check all methods: ----------------------------------------------

## MM hack: extract them from the randDAG() function definition
body. <- body(randDAG)
is.switch <- function(P) !is.symbol(P) && identical(as.symbol("switch"), P[[1]])
switchCall <- body.[vapply(body., is.switch, NA)][[1]]
stopifnot(identical(as.symbol("switch"), switchCall[[1]]))
(rDAGmeths <- names(switchCall)[-c(1:2, length(switchCall))])
rDAGall <- function(n, d, ...)
    sapply(rDAGmeths, function(meth) randDAG(n,d, method=meth, ...),
           simplify=FALSE)
set.seed(37)
rD.10.4 <- rDAGall(10, 4)
## with a low-level warning
rD.10.4 # looks ok

require(graph)
stopifnot(vapply(rD.10.4, isDirected, NA))

stopifnot(identical(
    lapply(rD.10.4, leaves, "out"),
    list(er = "3", regular = c("1", "5", "6"), watts = c("3", "4", "6"),
         bipartite = c("1", "2", "5"), barabasi = c("4", "8"),
         geometric = c("4", "7"), power = c("4", "5", "9"),
         interEr = c("3", "7"))
))

stopifnot(identical(
    lapply(rD.10.4, leaves, "in"),
    list(er = c("1", "4", "7"), regular = c("3", "7", "10"),
         watts = c("1", "8"), bipartite = c("4", "6"),
         barabasi = c("6", "7"), geometric = c("5", "10"),
         power = c("2", "7"), interEr = c("8", "10"))
))

set.seed(47)
rD.12.2 <- rDAGall(12, 2)
stopifnot(vapply(rD.12.2, isDirected, NA),
          vapply(rD.12.2, numNodes, 1) == 12,
          identical(vapply(rD.12.2, numEdges, 1),
                    setNames(c(9, 12, 12, 11, 11, 11, 13, 8), rDAGmeths))
          )

## Use the output here
require(Matrix)
lapply(rD.10.4, function(g) as(as(g, "Matrix"),"nMatrix"))
lapply(rD.12.2, function(g) as(as(g, "Matrix"),"nMatrix"))

##---------------------------------------------------------------------------

## check weights
set.seed(123)
n <- 100
g <- randDAG(n=n,d=3, wFUN=list(runif,min=0,max=1))
g
m <- wgtMatrix(g)
stopifnot(sum(m != 0) == 137)
v <- as.numeric(m)
v <- v[v!=0]
## dput(as.vector(summary(v, digits=7)))
stopifnot(all.equal(as.vector(summary(v, digits=7)),
                    c(0.008103577, 0.2589966, 0.5287397,
                      0.5232445, 0.8159941, 0.9915566)))
ct <- cut(x=v, breaks=seq(0,1,by=0.1))
stopifnot(all.equal(chisq.test(as.numeric(table(ct)), p = rep(0.1,10))$p.value,
                    0.3101796548))

## check generation of negative weights (fixed Bug)
set.seed(123)
tmp1 <- randDAG(3,2,wFUN = list(runif, min = 2, max = 2))
all( unlist(tmp1@edgeData@data) == 2 )
set.seed(123)
tmp2 <- randDAG(3,2,wFUN = list(runif, min = -2, max = -2))
all( unlist(tmp2@edgeData@data) == -2 )

