#segments.panel.rcomp <- function(X,Y,...,steps=30) {
#  function(what,...) {
#    if( !is.null(colnames(what)) )
#      what <- colnames(what)
#    X <- gsi.margin(X,what)
#    Y <- gsi.margin(Y,what)
#    segments.rcomp(X,Y,...,steps=steps)
#  }
#}

#segments.panel.aplus <- function(X,Y,...,steps=30) {
#  function(what,...) {
#    if( !is.null(colnames(what)) )
#      what <- colnames(what)
#    X <- X[,what]
#    Y <- Y[,what]
#    segments.aplus(X,Y,...,steps=steps)
#  }
#}

#segments.panel.rplus <- function(X,Y,...,steps=30) {
#  function(what,...) {
#    if( !is.null(colnames(what)) )
#      what <- colnames(what)
#    X <- X[,what]
#    Y <- Y[,what]
#    segments.rplus(X,Y,...,steps=steps)
#  }
#}

#segments.panel.rmult <- function(X,Y,...,steps=30) {
#  function(what,...) {
#    if( !is.null(colnames(what)) )
#      what <- colnames(what)
#    X <- X[,what]
#    Y <- Y[,what]
#    segments.rplus(X,Y,...,steps=steps)
#  }
#}

gsi.margin <- function(X,...) UseMethod("gsi.margin",X)

gsi.margin.acomp <- function(X,what,...,margin="acomp") {
  if( margin == "sub" )
    acomp(X,what)
  else if( margin=="rcomp" )
      rcompmargin(X,what)
  else if( margin=="acomp")
    acompmargin(X,what)
  else {
    if( !is.numeric(what) )
      what <- match(what,colnames(X))
    if( !is.numeric(margin))
      margin <- match(margin,colnames(X))
    acomp(X,c(what,margin))
  }
}

gsi.margin.rcomp <- function(X,what,...,margin="rcomp") {
  if( margin == "sub" )
    acomp(X,what)
  else if( margin=="rcomp" )
    rcompmargin(X,what)
  else if( margin=="acomp")
    acompmargin(X,what)
  else {
    if( !is.numeric(what) )
      what <- match(what,colnames(X))
    if( !is.numeric(margin))
      margin <- match(margin,colnames(X))
    rcomp(X,c(what,margin))
  }
}

gsi.margin.aplus <- function(X,what,...) {
  aplus(X,what)
}

gsi.margin.rplus <- function(X,what,...) {
  rplus(X,what)
}
