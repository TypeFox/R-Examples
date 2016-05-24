plot.lp.evmOpt <- function(x, main=NULL,
         pch= 1, ptcol =2 , cex=.75, linecol = 4 ,
         cicol = 1, polycol = 15, ...){
  family <- x$family
  x <- x$link
  
  if(dim(x)[1] == 1){
    stop("Need range of covariate values to plot linear predictors")
  }
  if(!any(colnames(x) == "phi.lo") ){
    stop("Please use ci.fit=TRUE in call to predict, to calculate confidence intervals")
  }

  makelp <- function(x, family){
      p <- family$param
      res <- vector("list", length=length(p))
      names(res) <- p
      for (i in 1:length(p)){
          res[[i]] <- as.matrix(x[, paste0(p[i], c("", ".lo", ".hi"))])
      }
      res
  }
  Ests <- makelp(data.frame(x), family)
  Names <- family$param
  cn <- colnames(x)
  which <- cn != "mu" & cn != "phi"    & cn != "xi" &
           cn != "mu.lo" & cn != "mu.hi" &
           cn != "phi.lo" & cn != "phi.hi" &
           cn != "xi.lo"  & cn != "xi.hi" &
           cn != "mu.se" & cn != "phi.se" & cn != "xi.se"

  X <- x[,which]
  if(is.null(dim(X))){
     X <- matrix(X)
     dimnames(X) <- list(dimnames(x)[[1]],dimnames(x)[[2]][which])
  }

  for(i in 1:length(Names)){
    for(j in 1:dim(X)[2]){
      if(length(unique(Ests[[i]][,1])) > 1){
        if(length(unique(X[,j])) > 1){
          ord <- order(X[,j])
          x <- X[ord,j]
          y <- Ests[[i]][ord,]
          plot(x, y[,1],type="n",ylab=Names[i],xlab=colnames(X)[j],main=main,ylim=range(y))

          if (polycol != 0){
            polygon(c( x,        rev(x)),
                    c(y[,2],rev(y[,3])),
                    col=polycol, border = FALSE) # Close polygon
          } else {
            lines(x, y[,2], col = cicol)
            lines(x, y[,3], col = cicol)
          }

          lines(x, y[,1], col = linecol[ 1 ] )
        }
      }
    }
  }
  invisible()
}

plot.lp.evmSim <- function(x, type="median", ...){
  if(dim(x$link)[1] == 1){
    stop("Need range of covariate values to plot linear predictors")
  }
  p <- x$family$param
  np <- length(p)
# re-format to same column structure as lp.evmOpt x
  ColIndexMeans <- 1+4*(0:(np-1))
  if(casefold(type) == "median"){
    offset <- 1
  } else if(casefold(type) == "mean") {
    offset <- 0
  } else {
    stop("type must be \"mean\" or \"median\" ")
  }
  which <- c(ColIndexMeans + offset,rep(1:2,np) + rep(ColIndexMeans+1,each=2), (4*np+1): dim(x$link)[2])
  x$link <- x$link[,which]
  colnames(x$link)[1:(3*np)] <-  c(p,paste(rep(p,each=2),rep(c(".lo",".hi"),np),sep=""))

  plot.lp.evmOpt(x,...)
}

plot.lp.evmBoot <- plot.lp.evmSim

