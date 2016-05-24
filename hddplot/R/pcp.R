"pcp" <-
function(x=USArrests,
           varscores=TRUE,
           cases="rows",
           center="vars",
           standardize=FALSE,
           scale.cases=1,
           log=FALSE,
           sc=1,
           reflect=c(1,1))
{
  x <- as.matrix(x)
  avv <- 0
  sdv <- 1
  casedim <- 2-as.logical(cases=="rows")
  vardim <- 3-casedim
  ## casedim=1 if rows are cases; otherwise casedim=2
  ## scale.cases=0 gives a pure rotation of the variables
  ## scale.cases=1 weights a/c the singular values
  ncases <- dim(x)[casedim]
  nvar <- dim(x)[vardim]
  if(is.null(sc))sc <- dim(x)[casedim]-1
  if(log)x <- log(x, base=2)
  if(standardize){
    avv <- apply(x, vardim, mean)
    sdv <- apply(x, vardim, sd)        
    x <- sweep(x, vardim, avv,"-")
    x <- sweep(x, vardim, sdv,"/")
  }
  else if(as.logical(match("vars", center, nomatch=0))){
    avv <- apply(x,vardim, mean)
    x <- sweep(x, vardim, avv,"-")}

  svdx <- La.svd(x)
  h <- NULL
  if(cases=="rows"){
    g <- sweep(svdx$u, 2, svdx$d^scale.cases, "*")*sqrt(sc)
    if(varscores)
      h <- t((svdx$d^(1-scale.cases)* svdx$vt ))/sqrt(sc)
  }
  else if(cases=="columns"){
    g <- sweep(t(svdx$vt), 2, svdx$d^scale.cases, "*")*sqrt(sc)
    if(varscores)
      h <- sweep(svdx$u, 2, svdx$d^(1-scale.cases),"*")/sqrt(sc)
  }
  invisible(list(g=g, rotation=h, av=avv, sdev=svdx$d/sqrt(ncases-1)))
}

