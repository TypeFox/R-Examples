## effect.sizes=x,y
plot.meta <- function(x, effect.sizes, add.margin=0.1, interval=0.95,
                      main="Effect Sizes and their Confidence Ellipses",
                      axis.labels=paste("Effect size ", effect.sizes, sep=""),                      
                      study.col="black", study.pch=19, study.min.cex=0.8, study.weight.plot=FALSE,
                      study.ellipse.plot=TRUE, study.ellipse.col="black", study.ellipse.lty=2,
                      study.ellipse.lwd=0.5,
                      estimate.col="blue", estimate.pch=18, estimate.cex=2,
                      estimate.ellipse.plot=TRUE, estimate.ellipse.col="red",
                      estimate.ellipse.lty=1, estimate.ellipse.lwd=2,
                      randeff.ellipse.plot=TRUE, randeff.ellipse.col="green",
                      randeff.ellipse.lty=1, randeff.ellipse.lwd=2,
                      univariate.plot=TRUE, univariate.lines.col="gray",
                      univariate.lines.lty=3, univariate.lines.lwd=1,
                      univariate.polygon.width=0.02, univariate.polygon.col="red",
                      univariate.arrows.col="green", univariate.arrows.lwd=2, diag.panel=FALSE,
                      xlim=NULL, ylim=NULL, ...) {

  if (!is.element("meta", class(x)))
    stop("\"x\" must be a class of \"meta\".")

  no.y <- x$no.y
  if (no.y==1) stop("There must be at least TWO effect sizes.\n")
  if (x$no.x!=0) warning("There are predictors in the model.\nThe plot is based on the intercepts.\n")
  if (missing(effect.sizes)) effect.sizes <- seq(1, no.y)

  ## Expand main to match the no. of plots
  no.of.plots <- length(effect.sizes)*(length(effect.sizes)-1)/2
  if (length(main)==1) main <- rep(main, no.of.plots)
   
  ## Procedure for more than 2 effect sizes
  if (length(effect.sizes)>2) {
   
    if (diag.panel==FALSE) {
      mat <- vec2symMat(seq(1, no.y*(no.y-1)/2))
      mat[upper.tri(mat)] <- 0
     } else {
      mat <- vec2symMat(seq(1, no.y*(no.y-1)/2), diag=FALSE)
      mat[upper.tri(mat)] <- 0
      Diag(mat) <- seq(no.y*(no.y-1)/2+1, no.y*(no.y+1)/2) 
    }

    ## Save the default par
    par.defaults <- par(no.readonly = TRUE) 
    layout(mat, respect=TRUE)
    
    n.plots <- 1
    for (i in 1:(length(effect.sizes)-1))
      for (j in (i+1):length(effect.sizes)) {
        my.effects <- effect.sizes[c(i,j)]
        ## Pass everything except effect.sizes=my.effects and main=main[n.plots]
        ## Added xlim and ylim arguments in v0.8-5. 
        plot.meta(x, effect.sizes=my.effects,
                  add.margin=add.margin, interval=interval, main=main[n.plots], axis.labels=axis.labels[c(i,j)],
                  study.col=study.col, study.pch=study.pch, study.min.cex=study.min.cex,
                  study.weight.plot=study.weight.plot, study.ellipse.plot=study.ellipse.plot,
                  study.ellipse.col=study.ellipse.col, study.ellipse.lty=study.ellipse.lty,
                  study.ellipse.lwd=study.ellipse.lwd,
                  estimate.col=estimate.col, estimate.pch=estimate.pch, estimate.cex=estimate.cex,
                  estimate.ellipse.plot=estimate.ellipse.plot, estimate.ellipse.col=estimate.ellipse.col,
                  estimate.ellipse.lty=estimate.ellipse.lty, estimate.ellipse.lwd=estimate.ellipse.lwd,
                  randeff.ellipse.plot=randeff.ellipse.plot, randeff.ellipse.col=randeff.ellipse.col,
                  randeff.ellipse.lty=randeff.ellipse.lty, randeff.ellipse.lwd=randeff.ellipse.lwd,
                  univariate.plot=univariate.plot, univariate.lines.col=univariate.lines.col,
                  univariate.lines.lty=univariate.lines.lty, univariate.lines.lwd=univariate.lines.lwd,
                  univariate.polygon.width=univariate.polygon.width,
                  univariate.polygon.col=univariate.polygon.col,
                  univariate.arrows.col=univariate.arrows.col, univariate.arrows.lwd=univariate.arrows.lwd,
                  diag.panel=FALSE, xlim=xlim, ylim=ylim, ...)
        n.plots <- n.plots+1
      }
    ## Reset the default par
    par(par.defaults)
    
  ## End of procedure for length(effect.sizes)>2   
  } else {
    ## Only two effect sizes
    if ( (length(effect.sizes)==2)&(diag.panel==TRUE) ) {
      layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(1,1), c(1,1), respect=TRUE)
    }
    ## Label of "Intercept" is used in meta object  
    my.effects <- paste("Intercept", effect.sizes, sep="")
    ES <- coef(x)[my.effects]
    ACov <- vcov(x)[my.effects, my.effects]
    # Fixme: diff nos of rand.effects or fixed-effects only (done?)
    my.rand <- vech(outer(effect.sizes, effect.sizes,
                          function(x,y) { paste("Tau2_",x,"_",y,sep="")}))
    # Names of all parameter estimates
    para.names <- summary(x$mx.fit)$parameters$name
    RE <- rep(NA, 3)
    # Extract variance components in the parameter estimates
    RE[my.rand %in% para.names] <- coef(x)[my.rand[my.rand %in% para.names]]
    RE <- vec2symMat(RE)
    #### Fixme
    if (all(!is.na(Diag(RE))) & is.na(RE[2,1]) ) RE[2,1] <- RE[1,2] <- 0.0
    dimnames(RE) <- list(my.effects, my.effects)

    ## y, x, var(y), cov(x,y), var(x)
    ## Index for positions of the conditional sampling variances
    rand.ind <- vec2symMat(seq(1, no.y*(no.y+1)/2))+ no.y
    ellipse.pt <- lapply(split(x$data, 1:nrow(x$data)),
                         function(x) {rand.val <- c(x[rand.ind[effect.sizes[1],effect.sizes[1]]],
                                                    x[rand.ind[effect.sizes[1],effect.sizes[2]]],
                                                    x[rand.ind[effect.sizes[2],effect.sizes[2]]])
                                      ellipse( vec2symMat(rand.val),
                                                centre=c(x[effect.sizes[1]],x[effect.sizes[2]]))})
    
    ## study.cex is proportional to sqrt(1/det(ACOV)).
    ## Thus, the plot area (cex) is proportional to 1/det(ACOV).
    if (study.weight.plot==TRUE) {
      weigh <- sapply(split(x$data, 1:nrow(x$data)),                   
                      function(x) {rand.val <- c(x[rand.ind[effect.sizes[1],effect.sizes[1]]],
                                                 x[rand.ind[effect.sizes[1],effect.sizes[2]]],
                                                 x[rand.ind[effect.sizes[2],effect.sizes[2]]])
                                   sqrt(1/det( vec2symMat(rand.val) ))})
      study.cex <- weigh*study.min.cex/min(weigh)
    } else {
      study.cex <- study.min.cex
    }

    ## Added xlim and ylim arguments in v0.8-5
    if (is.null(xlim)) {
        xlim <- sapply(ellipse.pt, function (x) { x <- data.frame(x); x$x })
        xlim <- c(min(xlim, na.rm=TRUE), max(xlim, na.rm=TRUE)) + c(-add.margin, 0)
    }
    if (is.null(ylim)) {
        ylim <- sapply(ellipse.pt, function (x) { x <- data.frame(x); x$y })
        ylim <- c(min(ylim, na.rm=TRUE), max(ylim, na.rm=TRUE)) + c(-add.margin, 0)
    }
    
    # remove incomplete x or y
    complete <- with(x, complete.cases(data[, effect.sizes[1]], data[, effect.sizes[2]]))
    my.x <- x$data[complete, effect.sizes[1]]
    my.y <- x$data[complete, effect.sizes[2]]
    plot(my.y~my.x, xlim=xlim, ylim=ylim, col=study.col, pch=study.pch, cex=study.cex,
         xlab=axis.labels[1], ylab=axis.labels[2], main=main, ...)
    if (study.ellipse.plot==TRUE) {
      for (i in 1:length(ellipse.pt)) {
        points(ellipse.pt[[i]], type="l", col=study.ellipse.col, lty=study.ellipse.lty,
               lwd=study.ellipse.lwd)
        }
    }
  
    points(x=ES[1], y=ES[2], col=estimate.col, pch=estimate.pch, cex=estimate.cex)
    if (estimate.ellipse.plot==TRUE) {
      points(ellipse(ACov, centre=ES), type="l", col=estimate.ellipse.col,
             lty=estimate.ellipse.lty, lwd=estimate.ellipse.lwd)
    }  
    ## Plot ellipse for random effects only if no missing in RE
    if ( randeff.ellipse.plot==TRUE && all(!is.na(Diag(RE))) ) {
      points(ellipse(RE, centre=ES), type="l", col=randeff.ellipse.col,
             lty=randeff.ellipse.lty, lwd=randeff.ellipse.lwd)
    }
    ## Plot univariate if at least one RE is present
    if ( univariate.plot==TRUE ) {
      intervals <- c((1-interval)/2, (1+interval)/2)      
      x.m <- ES[1]
      x.se.limit <- qnorm(intervals, mean=x.m, sd=sqrt(ACov[1,1]))
      abline(v=x.se.limit, col=univariate.lines.col, lty=univariate.lines.lty, lwd=univariate.lines.lwd)
      polygon(c(x.se.limit[1],x.m,x.se.limit[2],x.m),
              c(ylim[1], ylim[1]-univariate.polygon.width, ylim[1], ylim[1]+univariate.polygon.width),
              col=univariate.polygon.col)
      # fixme
      if ( !is.na(RE[1,1]) ) {
        x.tau.limit <- qnorm(intervals, mean=x.m, sd=sqrt(RE[1,1]))
        abline(v=x.tau.limit, col=univariate.lines.col, lty=univariate.lines.lty, lwd=univariate.lines.lwd)
        arrows(x.tau.limit[1], ylim[1], x.tau.limit[2], ylim[1], code=3, length=0.1,
               col=univariate.arrows.col, lwd=univariate.arrows.lwd)
      }
      y.m <- ES[2]
      y.se.limit <- qnorm(intervals, mean=y.m, sd=sqrt(ACov[2,2]))
      abline(h=y.se.limit, col=univariate.lines.col, lty=univariate.lines.lty, lwd=univariate.lines.lwd)
      polygon(c(xlim[1], xlim[1]-univariate.polygon.width, xlim[1], xlim[1]+univariate.polygon.width),
              c(y.se.limit[1],y.m,y.se.limit[2],y.m),  col=univariate.polygon.col)      
      # fixme
      if ( !is.na(RE[2,2]) ) {
        y.tau.limit <- qnorm(intervals, mean=y.m, sd=sqrt(RE[2,2]))      
        abline(h=y.tau.limit, col=univariate.lines.col, lty=univariate.lines.lty, lwd=univariate.lines.lwd)
        arrows(xlim[1], y.tau.limit[1], xlim[1], y.tau.limit[2], code=3, length=0.1,
               col=univariate.arrows.col, lwd=univariate.arrows.lwd)
      }
    }
  }

}
