rda.plotmat <- function(mat, main='', len=5, ncolor=100, se=TRUE, 
                        pos=NULL, minpos=NULL, nice=FALSE){
  if(!(is.numeric(mat) & is.matrix(mat))){
    stop("You must provide an numeric matrix to plot.")
  }

  if(len <= 1) {
    stop("The length of labels must be > 1.")
  }
  
  if(nice & ((ncol(mat) <= 3) || (nrow(mat) <= 3))){
    alph <- as.numeric(dimnames(mat)[[1]])
    delt <- as.numeric(dimnames(mat)[[2]])    
    if((ncol(mat) <= 3) & (nrow(mat) > 3)) {
      par(mfrow=c(ncol(mat), 1))
      for(i in 1:ncol(mat)){
        plot(alph, mat[, i], xlab=expression(paste(alpha)), ylab='', 
             main=substitute(paste(list(x), ": ", Delta) == list(y), 
                             list(x=main, y=delt[i])))
      }
      par(mfrow=c(1, 1))
    }
    if((ncol(mat) > 3) & (nrow(mat) <= 3)) {
      par(mfrow=c(nrow(mat), 1))
      for(i in 1:nrow(mat)){
        plot(delt, mat[i, ], xlab=expression(paste(Delta)), ylab='', 
             main=substitute(paste(list(x), ": ", alpha) == list(y), 
                             list(x=main, y=alph[i])))
      }
      par(mfrow=c(1, 1))
    }
    if((ncol(mat) <= 3) & (nrow(mat) <= 3)) {
      cat("The dimensions of your data are too small to plot.\n")
      print(mat)
    }      
  }
  else{
    old.par <- par(no.readonly = TRUE)

    a <- nrow(mat)
    alpha <- as.numeric(dimnames(mat)[[1]])
    delta <- as.numeric(dimnames(mat)[[2]])
    jump.alpha <- (max(alpha)-min(alpha))/(len-1)
    jump.delta <- (max(delta)-min(delta))/(len-1)
  
    layout(matrix(c(1, 2, 1, 3), 2, 2), heights=c(8, 2), widths=c(3, 7))
    image(delta, alpha, t(mat[rev(seq(a)), ]), yaxt='n',
          ylim=c(min(alpha)-jump.alpha*0.1, max(alpha)+jump.alpha*0.1),
          xlab=expression(paste(Delta)),
          ylab=expression(paste(alpha)),
          main=main, col=terrain.colors(ncolor))
    alpha.lab <- round(seq(min(alpha), max(alpha), len=len), 3)
    axis(2, at=alpha.lab, labels=rev(alpha.lab), line=1)

    if(se){
      if(is.null(nrow(pos))) {
        warning("The one-standard error boundary points are not plotted.")
      }
      else {
        for(i in 1:nrow(pos)){
          points(delta[pos[i, 2]], max(alpha)-alpha[pos[i, 1]], pch=16, 
                 cex=1, col=4)
        }
      }
      if(is.null(nrow(minpos))) {
        warning("The mimimal CV error points are not plotted.")
      }
      else {
        for(i in 1:nrow(minpos)){
          points(delta[minpos[i, 2]], max(alpha)-alpha[minpos[i, 1]], pch=1, 
          cex=1.5, col=2)
        }
      }
    }
    
    par(mar=c(5, 2, 0, 0))
    plot(1:10, rep(0, 10), type='n', axes=FALSE, xlab="",
         ylab="")
    legend(3, 1, pch=c(1, 16), col=c(2, 4), legend=c("Min CV Err.", 
           "1-SE CV Err."), cex=0.9)
 
    par(mar=c(5, 0, 0, 2))
    plot(seq(min(mat), max(mat), len=ncolor),
         rep(0, ncolor), pch=15, axes=FALSE, xlab="Color Code",
         ylab="", col=terrain.colors(ncolor))
    axis(1, at=seq(min(mat), max(mat), len=len),
         labels=round(seq(min(mat), max(mat), len=len), 3))
    on.exit(par(old.par))
    invisible()
  }
}
  
