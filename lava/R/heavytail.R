##' @export
`heavytail` <- function(x,...) UseMethod("heavytail")
##' @export
"heavytail<-" <- function(x,...,value) UseMethod("heavytail<-")

##' @export
"heavytail<-.lvm" <- function(x,...,value) {
  if (inherits(value,"formula")) {
    return(heavytail(x,all.vars(value),...))
  }
  heavytail(x, value, ...)
}

##' @export
`heavytail.lvm` <-
function(x,var=NULL,df=1,...) {
  if (is.null(var)) {
    htidx <- x$attributes$heavytail
    if (length(htidx)>0 && any(htidx!=0)) {
      res <- htidx[htidx>0]
      attributes(res)$couple <- unlist(x$attributes$heavytail.couple)[htidx>0]
      return(res)
    }
    return(NULL)
  }
  couples <- attributes(heavytail(x))$couple
  newval <- 1
  if (length(couples)>0) newval <- max(couples)+1
  x$attributes$heavytail.couple[var] <- newval
  x$attributes$heavytail[var] <- df
  return(x)
}

heavytail.init.hook <- function(x,...) {
  x$attributes$heavytail <- list()
  x$attributes$heavytail.couple <- list()
  return(x)
}

heavytail.sim.hook <- function(x,data,...) {
  n <- nrow(data)
  hvar <- heavytail(x)
  if (length(hvar)==0) return(data)
  couples <- unique(attributes(hvar)$couple)
  h.type <- list()
  for (j in couples)
    h.type <- c(h.type, list( hvar[(which(attributes(hvar)$couple==j))]))
  for (i in seq_along(couples)) {
    df <- hvar[[i]][1]
    Z <- rchisq(n,df=df)/df
    for (v in names(h.type[[i]])) {
      data[,v] <- data[,v]/sqrt(Z)
    }
  }
  return(data)
}
