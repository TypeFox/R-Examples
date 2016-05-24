# updated in february 2013 by Raimon
pairwisePlot <- function(X,Y,...) UseMethod("pairwisePlot",X)

# updated in february 2013 by Raimon                                 
pairwisePlot.default <- function(X,Y=X,...,xlab=deparse(substitute(X)),ylab=deparse(substitute(Y)),nm=c(length(Y),length(X)),panel=plot,
                                 add.line=FALSE, line.col=2,add.robust=FALSE,rob.col=4) {
  ylab
  xlab
  if( is.data.frame(X) ) {
    X<-X # data.frames are lists already
  } else if( is.list(X) ) {
    if( !missing(xlab) || is.null(names(X)) )
      names(X) <- if( length(xlab) < length(X) ) paste(xlab,1:length(X),sep="") else xlab
  } else {
    if( length(dim(X))== 0 ) {
      X <- t(t(X))
      colnames(X) <- xlab
    }
    if( is.null(colnames(X)) )
      colnames(X) <- if( length(xlab) < ncol(X) ) paste(xlab,1:ncol(X),sep="") else xlab
    cn <- 1:ncol(X)
    names(cn)<-colnames(X)
    X <- lapply(cn,function(i) X[,i])
  }
                                        #else {
                                        #stop("Unhandled type of Y in pairwisePlot:",class(Y))
                                        #}
  if( is.data.frame(Y) ){
    Y<-Y # data.frames are lists already
  } else if( is.list(Y) ) {
    if( !missing(ylab) || is.null(names(Y)) )
      names(Y) <- if( length(ylab) < length(Y) ) paste(ylab,1:length(Y),sep="") else ylab
  } else {
    if( length(dim(Y))== 0 ) {
      Y <- t(t(Y))
      colnames(Y) <- ylab
    }
    if( is.null(colnames(Y)) )
      colnames(Y) <- if( length(ylab) < ncol(Y) ) paste(ylab,1:ncol(Y),sep="") else ylab
    cn <- 1:ncol(Y)
    names(cn)<-colnames(Y)
    Y <- lapply(cn,function(i) Y[,i])
  }
                                        #else {
                                        # stop("Unhandled type of Y in pairwisePlot:",class(Y))
                                        #}
  if( !is.null(nm) ) {
    opar <- par(mfrow=nm)
    on.exit(par(opar))
  }
  withformula = length( grep("formula", methods(panel) ) )>0
  for(j in 1:length(Y) )
    for(i in 1:length(X)) {
      # if there is a formula interface for the desired panel, use it
      if(withformula){
        panel(Y[[j]]~X[[i]],
            xlab=names(X)[i],
            ylab=names(Y)[j],...)
      }else{
        panel(X[[i]], Y[[j]],
              xlab=names(X)[i],
              ylab=names(Y)[j],...)
      }  
      if( !any(is.factor(Y[[j]]), is.factor(X[[i]]) ) ){
       if(add.line) abline( lm(Y[[j]]~X[[i]]), col=line.col )
       if(add.robust) abline(lmrob(Y[[j]]~X[[i]])$coefficients, col=rob.col, lwd=2)
      } 
    }
}

# updated in february 2013 by Raimon
pwlrPlot = function(x, y, ..., add.line=FALSE, line.col=2, add.robust=FALSE, rob.col=4){
  if(is.null(y)){
    stop("argument y cannot be empty")
  }
  if(is.null(x)){
    stop("argument x cannot be empty")
  }
  tk = "acomp" %in% c(class(x), class(y))
  if(!tk){
    stop("one of x or y must be an acomp composition!")
  }
  X = list(x,y)
  # check which of the two is not compositional, and it is a vector
  funcond = function(z) any(is.data.frame(z) & ncol(z)==1, 
                            is.matrix(z) & ncol(z)==1, 
                            (is.factor(z) || is.numeric(z)) && class(z)!="acomp" )
  tk = sapply(X, funcond)
  covarisfactor = is.factor(X[tk][[1]])
  if(!any(tk)){
    stop("one of x or y must be a unique covariable, either a one-column data.frame or matrix, a factor or a vector")
  }
  if( any(add.line, add.robust) & covarisfactor){
    add.line = FALSE
    add.robust = FALSE
    warning("regression line meaningless with factor covariable; add.line and add.robust set to FALSE")
  }
  # apply a pwlr to the compositional part
  D = ncol(X[!tk][[1]])
  partnames = colnames(X[!tk][[1]])
  covariablename = c(deparse(substitute(x)), deparse(substitute(y)))[tk]
  pwlr = function(i,j) log( X[!tk][[1]][,i]/X[!tk][[1]][,j] )
  X[!tk][[1]] <- mapply(i=rep(1:D,times=D),j=rep(1:D, each=D), pwlr)
  # ensure that the covariable is represented as a vector, then repeat as many times as D*D
  aux =  unlist(X[tk][[1]], recursive=FALSE)
  dim(aux) = NULL
  X[tk][[1]] = rep(aux, times=D*D)
  dim(X[tk][[1]]) = dim(X[!tk][[1]])
  # set graphical parameters now, store prior back as default
  opar <- par()
  on.exit(par(opar))
  par(mfrow=c(D,D), mar=c(1,1,0,0),oma=c(3,3,2,2))
  # plot!
  for(i in 1:D){
    for(j in 1:D){
      k = (i-1)*D+j
      if(i!=j){
        plot( X[[2]][,k]~X[[1]][,k], xaxt=ifelse(i==D,"s","n"), yaxt=ifelse(j==1,"s","n"),... )
        if(add.line) abline(lm( X[[2]][,k]~X[[1]][,k] ), col=line.col, lwd=2)
        if(add.robust) abline(lmrob(X[[2]][,k]~X[[1]][,k] )$coefficients, col=rob.col, lwd=2)
        if(!is.factor(X[[1]][,k]) ) axis(side=3, labels=(i==1))
        if(i==1 & is.factor(X[[1]][,k]) ) axis(side=3, at=1:length(levels( X[[1]][,k])), labels=abbreviate(levels( X[[1]][,k] )))
        axis(side=4, labels=(j==D))
      }else{
        plot( 0, 0, xaxt="n", yaxt="n", bty="n", type="n", xlab="", ylab=""  )
        text(0,0, labels=partnames[i],cex=1.5)
      }
    }
  }
  titols = rep("pairwise logratio", 2)
  titols[tk] = covariablename
  mtext(side=1:2, titols, cex=1.25, outer=TRUE, line=1.5)
}