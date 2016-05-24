
.parminfo <- function(parms){
  UseMethod(".parminfo")
}

.parminfo.pms <- function(parms){
  c(list(cont.names=rownames(parms$mu), disc.names=names(dimnames(parms$p)),disc.levels=dim(parms$p)),
    parms[c("gentype","Ad.idx","Ac.idx")])
}

.parminfo.CGstats <- function(parms){
  list(cont.names=rownames(parms$center), disc.names=names(dimnames(parms$n.obs)),disc.levels=dim(parms$n.obs))
}

.disc.levels <- function(parms){
  dim(parms[[1]])
}
