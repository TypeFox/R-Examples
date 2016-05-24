plotTheoretical <-
  function(ploidy.level=8, seg.ratios=NULL,
           n.components=NULL, expected.segratio=NULL,
           proportions=c(0.65,0.2,0.1,0.03,0.01,0.01, 0, 0),
           n.individuals=200, xaxis=c("raw","logit"), 
           type.parents=c("heterogeneous","homozygous"),
           xlim=c(0,1), NCLASS=NULL,
           xlab="Segregation Ratio", ylab="Density", density.plot=FALSE,
           fitted.lwd=2, fitted.col="blue", cex=1, warnings = TRUE,
           main=NULL, ...)
{
  ## Purpose: plot theoretical density and histogram of observed
  ##          segregation ratios (if provided) on observed p-value scale
  ##          with components corresponding to dosage classes

  ## Arguments:
  ## ploidy.level: the number of homologous chromosomes, either as numeric
  ##                or as a character string
  ## seg.ratios:   segregation ratios as class 'segRatio' (optional)
  ## n.components: no. of components in model (default: ploidy.level/2)
  ## expected.segratio: expected segregation ratios. determined by theory
  ##                    (Haldane, 1930) given ploidy level if not supplied
  ## proportions: proportions of markers in each dosage class
  ## n.individuals: no. individuals per marker which is determined from
  ##                seg.ratios or set to 200 if not supplied
  ## xlim:        c(lower,upper) limits for segregation ratios
  ## NCLASS:      number of classes for histogram
  ## ylab:        y-axis label
  ## density.plot: add smoothed  density to plot
  ## xlab:        x-axis label
  ## eg ... ylab, main etc

  if (length(seg.ratios)==0) {
    plot.sr <- FALSE
  } else {
    plot.sr <- TRUE
    if (class(seg.ratios) != "segRatio") {
      stop("'seg.ratios' must be of class 'segRatio'")
    }
    n.individuals <- seg.ratios$n.individuals
    y <- seg.ratios$seg.ratio
    if (length(NCLASS)==0)
      NCLASS <- min(max(nclass.Sturges(y),round(length(y)/6)),25)
  }

  type <- match.arg(type.parents)
  pl.type <- match.arg(xaxis)

  if (pl.type == "raw") {
    x.seq <- c(0:n.individuals)/n.individuals  # plotting binomial densities
  } else {
    x.seq <- c(1:(n.individuals-1))/n.individuals  # plotting on logit scale
  }
  
  if (length(expected.segratio) == 0) {
    E.segRatio <- expected.segRatio(ploidy.level, type.parents=type)
  } else {
    E.segRatio$ratio <- expected.segratio
  }
  Ep <- E.segRatio$ratio

  if(length(n.components)==0)
    n.components <- length(E.segRatio$ratio)
    
  if (length(proportions) != n.components) {
    if (length(proportions)< n.components) {
      eta <- proportions
      n.components <- length(eta)
      if (warnings)
        warning(paste("'proportions' truncated to length:",n.components))
    } else {
      eta <- proportions[1:n.components]
      if (warnings)
        warning(paste("'n.components' set to:",n.components))
    }
  } else {
    eta <- proportions
  }
  
  if (sum(eta) != 1) {
    eta <- eta/sum(eta)            # normalise to sum to 1
    cat("Warning: component proportions normalised, now:\n")
    print(eta)
  }

  specific.set <- list(axis.text=list(col=fitted.col),
                       axis.line=list(col=fitted.col),
                       axis.text=list(cex=0.8*cex),
                       par.main.text=list(cex=1.2*cex,col=fitted.col),
                       par.xlab.text=list(cex=cex,col=fitted.col),
                       par.ylab.text=list(cex=cex,col=fitted.col))


  if (length(main)==0) {
    main=paste("Theoretical Binomial mixture density for", n.individuals,
      "individuals with ploidy:",ploidy.level)
  }
  
  ## generate line data for components and total assuming Binomial distribution
  
  comp.lines <- vector("list", n.components+1)
  names(comp.lines) <- paste(c("Total",1:n.components))
  comp.lines[[1]] <- 0*x.seq
  for (kk in 1:n.components){
    comp.lines[[kk+1]] <- eta[kk]*n.individuals*
      dbinom(x.seq*n.individuals,size=n.individuals,prob=Ep[kk])
    comp.lines[[1]] <-  comp.lines[[1]] +  comp.lines[[kk+1]]
  }

  ## modify x axis if set to logit - also lines so as area right
  ## if plotting histogram of segeregation ratios
  ## NB: This is a dodge and so may not work for small n or truncated x-axis
  
  if (pl.type == "logit") {
    width <- x.seq[2]-x.seq[1]
    ##    old.total <- width*sum(comp.lines[[1]][-1])  # should be ~1
    orig.seq <- x.seq
    x.seq <- gtools::logit(x.seq)
    xlab=paste("logit(",xlab,")",sep="")
    
    new.area <- old.area <- rep(0,length(comp.lines))
    new.width <- diff(x.seq)
    for (iii in 1:length(comp.lines)){
      old.area[iii] <- width*sum(comp.lines[[iii]][-1])
      new.area[iii] <- sum(new.width*
                           comp.lines[[iii]][-1])
    }
    ##old.total <- old.area[1]
    ##new.total <- new.area[1]
    ##comp.lines[[1]] <- old.total/new.total * comp.lines[[1]]
    ## recalculate 'Total' line as components will change
    comp.lines[[1]] <- 0*comp.lines[[1]]
    for (iii in 1:n.components) {
      comp.lines[[iii+1]] <-  (old.area[iii+1]/new.area[iii+1]) *
        comp.lines[[iii+1]]
      comp.lines[[1]] <- comp.lines[[1]] + comp.lines[[iii+1]]
    }
    warning("Binomial mixture density is only approximate on logit scale")
    
    if (plot.sr) {
      ## first attempt - too global??
      ##tmp.hist <- hist(y, nclass=NCLASS, plot = FALSE)
      ##area <-  sum(tmp.hist$mids * tmp.hist$density)
      y <- gtools::logit(y)
      ##tmp.hist <- hist(y, nclass=NCLASS, plot = FALSE)
      ##area.logit <-  sum(tmp.hist$mids * tmp.hist$density)

      if (identical(xlim , c(0,1)))
        xlim <- c(min(y),max(y))
    } else {
      if (identical(xlim , c(0,1)))
        xlim <- c(min(x.seq),max(x.seq))
    }
  }

  if (plot.sr) {
    
    xhist <- lattice::histogram(y, nint=NCLASS, xlab=xlab, ylab=ylab,
                       type="density", col="lightgreen", main=main,
                       par.settings=specific.set,
                       panel = function(x, nint=nint, # main=main,
                         type=type,col=col, dp=density.plot,
                         seq=x.seq, cl=comp.lines,
                         ...) {
                         panel.histogram(x, nint=nint, type=type,col=col, ...)
                         for (kk in 1:length(cl)){
                           panel.lines(seq, cl[[kk]], col = fitted.col, lty=kk,
                                       lwd = fitted.lwd)
                         }                        
                         if (dp)
                           panel.densityplot(x,...)
                       } )
    return(xhist)
  } else {
    
    dplot <- lattice::xyplot( comp.lines[[1]] ~ x.seq, xlim=xlim, type = "l",
                    col = fitted.col, lty=1, lwd = fitted.lwd,
                    ylab=ylab, xlab=xlab, main=main,
                    seq=x.seq, cl=comp.lines, # par.settings=specific.set,
##                    key= legend(0.8,4, c("Total","1","2","3"),lty=c(1:4)),
                    panel = function(x, xlab=xlab, 
                      type=type, col=col, dp=density.plot,
                      seq=x.seq, cl=comp.lines,
                      ...) {
##              panel.xyplot(x, type=type, col=col, lty=lty, lwd=lwd, ...)
                      panel.xyplot(x, type=type, ...)
                      for (kk in 1:(length(cl)-1)){
                        panel.lines(seq, cl[[kk]], col = fitted.col, lty=kk,
                                    type="l", lwd = fitted.lwd)
                      }                        
                    } )
    return(dplot)
  }
}
