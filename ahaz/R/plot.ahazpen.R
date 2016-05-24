"plot.ahazpen"<-function(x, xvar = c("norm","lambda"), labels = FALSE, df = TRUE,
                          ylab = "Regression coefficients",xlab = xname,...)
  {
    ## Purpose: plot regularization path from 'ahazpen' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   x    : 'ahazpen' object
    ##   type : 'lambda' - first axis is lambda (log scale)
    ##           'coefficients' - first axis is L1 norm
    ##   label: show beta indices in right margin
    ##   df   : show degrees of freedom in top margin
    ##   ylab : label for y axis
    ##   ylab : label for x axis ("L1 norm" for xvar="norm", "lambda" otherwise)
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen

    xvar <- match.arg(xvar)

    if (length(x$lambda) == 1)
      stop("No coefficient paths to plot!")

    include <- ahaz.nzcoef(x$beta)
    
    if(x$nvars>1)
      beta<-t(as.matrix(x$beta[include,]))
    else
      beta<-matrix(x$beta[include,],ncol=1)
    

    if(!sum(include))
      stop("All coefficient paths identically zero!")
    
    switch(xvar,
           "norm"={xname<-"L1 norm";z <- apply(beta, 1, function(x){sum(abs(x))});log=""},
           "lambda"={xname<-expression(lambda);  z <- x$lambda; log="x"}
           )


    if(is.null(list(...)$type))
      matplot(z, beta,type = "l", xlab = xlab,ylab = ylab,log = log,lty=1,...)
    else
      matplot(z, beta,xlab = xlab,ylab = ylab,log = log,lty=1,...)
 
    if(labels)
      axis(4, at = beta[nrow(beta),], labels = as.character(include))
    if(df){
        mtext("# nonzero coefficients", side = 3, line = 2, adj = .5, cex = 0)
        vv <- (c(1,diff(x$df))!=0)
        axis(side = 3, at = z[vv], labels = paste(x$df[vv]), tick = TRUE, line = 0)
      }
  }



         
