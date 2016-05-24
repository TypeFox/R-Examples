
.median_survival_time <- function(x) {
    minmin <- function(y, xx) {
        if (any(!is.na(y) & y==.5)) {
            if (any(!is.na(y) & y <.5))
                .5*(min(xx[!is.na(y) & y==.5]) + min(xx[!is.na(y) & y<.5]))
            else
                .5*(min(xx[!is.na(y) & y==.5]) + max(xx[!is.na(y) & y==.5]))
        } else   min(xx[!is.na(y) & y<=.5])
    }
    med <- suppressWarnings(minmin(x$surv, x$time))
    return(med)
}

### get the recursive index
### obj is of class "partynode"
.get_path <- function(obj, i) {

    idx <- c()
    recFun <- function(node, i) {
        if (id_node(node) == i) return(NULL)
        idx <<- c(idx, which(names(unclass(node)) == "kids"))
        kid <- sapply(kids_node(node), id_node)
        nextid <- max(which(kid <= i))
        idx <<- c(idx, nextid)
        return(recFun(node[[nextid]], i))
    }
    out <- recFun(obj, i)
    return(idx)
}

### <TH> shall we export this functionality?
"nodeids<-" <- function(obj, value) UseMethod("nodeids<-")

"nodeids<-.party" <- function(obj, value) {

  id0 <- nodeids(obj)
  id1 <- as.integer(value)
  stopifnot(identical(id1, 1:length(id0)))

  idxs <- lapply(id0, .get_path, obj = node_party(obj))
  x <- unclass(obj)
  ni <- which(names(x) == "node")
  nm <- x$names
  for (i in 1:length(idxs))
      x[[c(ni, idxs[[i]])]]$id <- id1[i]
  class(x) <- class(obj)
  if (!is.null(nm)) 
      names(x) <- nm[id0]

  return(x)
}

"nodeids<-.constparty" <- function(obj, value) {

  id0 <- nodeids(obj)
  cls <- class(obj)
  class(obj) <- "party"
  nodeids(obj) <- value
  id1 <- nodeids(obj)
  obj$fitted[["(fitted)"]] <- 
      id1[match(fitted(obj)[["(fitted)"]], id0)]
  class(obj) <- cls
  obj
}
### </TH>
