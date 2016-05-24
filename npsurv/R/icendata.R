######################## #
# Interval-censored data #
######################## #

# Also allows for exact observations included.

icendata = function(x, w=1) {
  if(is.null(x)) return(NULL)
  if(is.icendata(x)) {
    if(all(w == 1)) return(x)
    w = rep(w, length = length(x$t) + nrow(x$o))
    if(length(x$t) > 0) x$wt = x$wt * w[1:length(x$wt)]
    if(nrow(x$o) > 0) x$wo = x$wo * w[length(x$wt)+1:nrow(x$o)]
    return(x)
  }
  z = vector("list", 7)
  names(z) = c("t", "wt", "o", "wo", "i1", "upper", "u")
  if(is.vector(x)) x = cbind(x, x)
  if(!is.matrix(x)) x = as.matrix(x)
  if(ncol(x) == 3) {w = w * x[,3]; x = x[,1:2]}
  if(length(w) != nrow(x)) w = rep(w, len=nrow(x))
  iw = w > 0
  w = w[iw]
  x = x[iw,,drop=FALSE]
  o = order(x[,1], x[,2])
  x = x[o,]
  w = w[o]
  id = c(TRUE, diff(x[,1]) > 0 | diff(x[,2]) > 0)
  id[is.na(id)] = FALSE            # for Inf's
  w = aggregate(w, by=list(group=cumsum(id)), sum)[,2]
  x = x[id,]
  i = x[,1] == x[,2]
  z$t = x[i,1]
  names(z$t) = NULL
  z$wt = w[i]
  z$o = x[!i,1:2,drop=FALSE]
  dimnames(z$o) = list(NULL, c("L","R"))
  z$wo = w[!i]
  z$upper = max(x[,1])
  z$i1 = z$t != z$upper
  z$u = sort(unique(c(0, pmin(c(x[,1], x[,2]), z$upper))))
  class(z) = "icendata"
  z
}

is.icendata = function(x) "icendata" %in% class(x)

# is.rightcensored.icendata = function(x) all(x$o[,2] == Inf)

expand.icendata = function(x) {
  if(!is.icendata(x)) x = icendata(x)
  z = vector("list", 7)
  names(z) = c("t", "wt", "o", "wo", "i1", "upper", "u")
  z$upper = x$upper
  if(length(x$t) > 0) {
    z$t = rep(x$t, x$wt)
    z$wt = rep(1, length(z$t))
    z$i1 = z$t != z$upper
  }
  else z$t = z$wt = numeric(0)
  if(nrow(x$o) > 0) {
    z$o = cbind(rep(x$o[,1], x$wo), rep(x$o[,2], x$wo))
    z$wo = rep(1, nrow(z$o))
    colnames(z$o) = c("L","R")
  }
  else {z$o = matrix(nrow=0, ncol=2); z$wo = numeric(0)}
  z$u = x$u
  class(z) = "icendata"
  z
}

