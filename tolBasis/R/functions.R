
# extract elements of a vector or a list repeating elements when necessary
.extract.cycle <- function(vector, index) {
  if(is.null(vector)) return(rep(NA, length(index)))
  len <- length(vector)
  if(len==0) return(rep(NA, length(index)))
  vector[[((index-1)%%len)+1]]
}

# merge two lists by name
.merge.list.2 <- function(x, y) {
  if(length(x)==0) return(y)
  if(length(y)==0) return(x)
  i = match(names(y), names(x))
  i = is.na(i)
  if(any(i)) x[names(y)[which(i)]] = y[which(i)]
  x
}
# merge two or more lists by name
.merge.list <- function(x, ...) {
  narg <- nargs()
  if(narg==1) return(x)
  arg <- c(list(x), list(...))
  aux <- .merge.list.2(arg[[narg-1]], arg[[narg]])
  if(narg==2) return(aux)
  do.call(".merge.list", c(arg[1:(narg-2)], list(aux)))
}

# extract a sublist of arguments using a prefix
.extract.prefix <- function(arg.list, prefix) {
  if(length(arg.list)==0) return(arg.list)
  sel <- sapply(names(arg.list), function(x) length(grep(paste0(prefix,".*"), x))>0)
  if(length(sel)==0) return(list())
  arg.sel <- arg.list[sel]
  names(arg.sel) <-  sapply(names(arg.sel), function(x) { substring(x, nchar(prefix)+2) })
  arg.sel
}

# exclude arguments from a list using a prefix
.exclude.prefixes <- function(arg.list, prefixes) {
  sel <- sapply(names(arg.list), function(x) 
    !any(sapply(prefixes, function(prefix) length(grep(paste0(prefix,".*"), x))>0)))
  if(length(sel)==0) list() else arg.list[sel]
}

# filter arguments from a list using names
.filter.names <- function(arg.list, nmes) {
  arg.names <- names(arg.list)
  if(is.null(arg.names)) arg.names <- rep("", length(arg.list))
  sel <- sapply(arg.names, function(x) 
    any(sapply(nmes, function(name) x==name)))
  if(length(sel)==0) list() else arg.list[sel]
}

# exclude arguments from a list using names
.exclude.names <- function(arg.list, nmes) {
  arg.names <- names(arg.list)
  if(is.null(arg.names)) arg.names <- rep("", length(arg.list))
  sel <- sapply(arg.names, function(x) 
    !any(sapply(nmes, function(name) x==name)))
  if(length(sel)==0) list() else arg.list[sel]
}

.flatten <- function(lst) {
  new <- list()
  if(length(lst)==0) return(list())
  for(i in 1:length(lst)) {
    if(inherits(lst[[i]], "list")) new <- c(new, .flatten(lst[[i]]))
    else new <- c(new, lst[i])
  }
  new
}
 
.exclude.null <- function(l) l[!sapply(l, is.null)]
.exclude.na <- function(l) l[!sapply(l, is.na)]
