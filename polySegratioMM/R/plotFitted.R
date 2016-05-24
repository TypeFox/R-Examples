plotFitted <-
  function(seg.ratios, summary.mixture, add.random.effect=TRUE,
           theoretical=FALSE, model=NULL, theory.col="red",
           xaxis=c("logit","raw"), ylim=NULL, NCLASS=NULL, n.seq=100,
           xlab="logit(Segregation Ratio)", ylab="Density", density.plot=FALSE,
           fitted.lwd=2, fitted.col="blue", bar.col="lightgreen", cex=1,
           warnings = FALSE, main=NULL, ...)
{
  ## Purpose: plot out histogram of observed segregation ratios on
  ##          logit scale along with scaled density of fitted
  ##          components corresponding to dosage classes using trellis

  ## Arguments:
  ## seg.ratios:   segregation ratios as class 'segRatio' (optional)
  ## summary.mixture:     mcmc summary data produce by 'readJagsMix'
  ## proportions: proportions of markers in each dosage class
  ## n.individuals: no. individuals per marker which is determined from
  ##                seg.ratios or set to 200 if not supplied
  ## xlim:        c(lower,upper) limits for segregation ratios
  ## NCLASS:      number of classes for histogram
  ## ylab:        y-axis label
  ## density.plot: add smoothed  density to plot
  ## xlab:        x-axis label
  ## eg ... ylab, main etc

  if (class(seg.ratios) != "segRatio")
    stop("'seg.ratios' must be of class 'segRatio'")
  
  if (class(summary.mixture) != "summarySegratioMCMC")
    stop("'summary.mixture' must be of class 'summarySegratioMCMC'")

  var.names <- rownames(summary.mixture$statistics)

  ## extract model parameters
  
  eta <- summary.mixture$statistics[grep("P\\[",var.names), "Mean"]
  if (warnings) {
    cat("P:\n")
    print(eta)
  }
  mu <-  summary.mixture$statistics[grep("mu\\[",var.names), "Mean"]
  names.sigma <- grep("sigma\\[",var.names)
  if (length(names.sigma)==0) {
    sigma <-  summary.mixture$statistics["sigma", "Mean"]
  } else {
    sigma <-  summary.mixture$statistics[grep("sigma\\[",var.names), "Mean"]
  }
  
  if (length(sigma)==1) { 
    if (warnings) cat("Common variances assumed\n")
    sigma <- rep(sigma,length(eta))
  } else {
    if (warnings) cat("Common variances not assumed\n")
  }

  ## if random.effect then add sigmab to sigma's

  if (add.random.effect){
    random.effect <- grep("sigmab",var.names)
    if (length(random.effect)==1){
      sigmab <- summary.mixture$statistics[grep("sigmab",var.names), "Mean"]
      sigma <- sigma + sigmab
    }
  }

  ## what to do if xaxis=="raw" ?? not yet!!
  pl.type <- match.arg(xaxis)
  if (pl.type=="raw")
    warning("Option \"raw\" not yet available for 'xaxis'")

  ## else when you wotk out raw!!    
  ## seg ratios on logit scale 
    
  y <- gtools::logit(seg.ratios$seg.ratio)
  if (length(ylim)!=2)
    ylim <- c(min(y),max(y))

  if (length(NCLASS)==0)  # my dodgey way of getting more classes 
    NCLASS <- min(max(nclass.Sturges(y),round(length(y)/6)),25)

  if(theoretical==TRUE) {
    if (class(model)=="modelSegratioMM") {
      ttt <- plotTheoretical(ploidy.level=model$E.segRatio$ploidy.level,
                             n.components=model$n.components,
                             n.individuals=seg.ratios$n.individuals,
                             proportions=eta, xaxis=pl.type)

##      if (pl.type=="logit") {
##        tmp.hist <- hist(seg.ratios$seg.ratio, nclass=NCLASS, plot = FALSE)
##        area <-  sum(tmp.hist$mids * tmp.hist$density)
##        tmp.hist <- hist(y, nclass=NCLASS, plot = FALSE)
##        area.logit <-  sum(tmp.hist$mids * tmp.hist$density)
##        warning("Adding theoretical density on logit - inaccurate in right tail")
##      }
      y.theory <- ttt$panel.args[[1]]$y # *area.logit/area
      x.theory <- ttt$panel.args[[1]]$x
    } else {
      stop('To plot theoretical distribution, supply model as class "modelSegratioMM"')
    }
  }
  
  ## function for overall fitted distribution
  
  totalCurve <- function(x, eta, mu, sigma, ...) {

    y <- eta[1]*dnorm(x, mean=mu[1],sd=sigma[1])
    for (i in 2:length(eta)) {
      y <-  y + eta[i]*dnorm(x, mean=mu[i],sd=sigma[i])
    }
    return(y)
  }
 
  specific.set <- list(axis.text=list(col=fitted.col),
                       axis.line=list(col=fitted.col),
                       axis.text=list(cex=0.8*cex),
                       par.main.text=list(cex=1.2*cex,col=fitted.col),
                       par.xlab.text=list(cex=cex,col=fitted.col),
                       par.ylab.text=list(cex=cex,col=fitted.col))

  test <-  histogram(y, nint=NCLASS)
  width <- test$panel.args.common$breaks[2]-test$panel.args.common$breaks[1]
  MULT <- width*test$packet.sizes

  if (length(main)==0) {
    main=paste("Fitted mixture density:",
      main=deparse(substitute(seg.ratios)))
  }
  
  xhist <- lattice::histogram(y, nint=NCLASS, main=main, ylab=ylab, xlab=xlab,
                     type="density",col=bar.col,
                     par.settings=specific.set,
                     panel = function(x,nint=nint,main=main,ylab=ylab,
                       xlab=xlab,type=type,col=col,
                       ...) {
                       panel.histogram(x,nint=nint,type=type,col=col, ...)
                       panel.curve( totalCurve(x, eta, mu, sigma) ,
                                   col = fitted.col, lty = 1,
                                   lwd = fitted.lwd, n=n.seq)
                       for (kk in 1:length(eta)){
                         panel.curve(eta[kk]*dnorm(x, mean=mu[kk],
                                                   sd=sigma[kk]),
                                     col = fitted.col, lty = kk+1,
                                     lwd = fitted.lwd, n=n.seq)
                       }
                       if (theoretical==TRUE) {
                         panel.lines(x.theory, y.theory, col=theory.col,
                                     lty=1, lwd = fitted.lwd)
                       }
                       if (density.plot)
                         panel.densityplot(x,...)
                     })
  return(xhist)
}

