######################################################################
#                                                                    #
# RELIABILITY DIAGRAM FOR A COLLECTION OF PROBABILITY FORECASTS      #
#                                                                    #
######################################################################
ReliabilityDiagram <- 
function(probs, obs, bins=10, nboot=500, 
         plot=FALSE, plot.refin=TRUE, 
         cons.probs=c(0.025, 0.975),attributes=FALSE) 
{
  #
  # Plot reliability diagram for a probability forecast
  #
  # Usage: ReliabilityDiagram(probs, obs, nbins, nboot)
  #
  # Arguments:
  #
  #    probs ... vector of length N, probs[k] has the predicted probability for
  #              the event obs[k] 
  #    obs ... obs[k] = 1 if the event happened at instance k, obs[k] = 0
  #            otherwise
  #    bins ... either scalar: number of equidistant bins to discretize the
  #                            forecast probabilities,
  #             or a vector: user-defined breakpoints of the bins; the `hist`
  #                          function will produce errors if these are not valid
  #    nboot ... number of bootstrap resamples for estimating consistency bars
  #              if nboot==0, no resampling is done and NAs are returned as 
  #              consistency bars
  #    plot ... boolean; whether to plot the reliability diagram
  #    plot.refin ... boolean; whether to plot the small refinement histogram
  #                   in lower right corner
  #    cons.probs ... a 2-vector, lower and upper confidence limit
  #    attributes ... logical, whether to plot attributes (polygon, 
  #                   no-resolution and no-skill line) 
  #                   
  # Return value:
  #
  #    a data frame of K+1 rows with the following columns:
  #
  #      * p.avgs     ... in-bin averages of the forecast probabilities
  #      * cond.probs ... observed conditional frequency of event, given i
  #      * cbar.lo    ... lower limit consistency of consistency bar[i], as specified by user
  #      * cbar.hi    ... upper limit consistency of consistency bar[i], as specified by user
  #
  # Author: 
  #
  #    Stefan Siegert 
  #    s.siegert@exeter.ac.uk 
  #    December 2013
  #
  # Example:
  #
  #    N <- 1000
  #    p <- rbeta(N, 1, 3)
  #    y <- rbinom(N, 1, p)
  #    rd <- ReliabilityDiagram(p, y, plot=TRUE)
  #    print(rd)
  #
  #
  # change log:
  #
  #  2015/06/15
  #  * added option `attributes` 
  #
  #  2013/12/02
  #  * manual definition of bin-breaks
  #  * manual definition of consistency intervals
  #  * sanity checks
  #  * parallel option for resampling
  #
  #  2013/10/31:
  #  * return summary data as data frame
  #  * added options `plot` and `plot.refin`
  #
  #  2013/08/20:
  #  * points are plotted at in-bin-averages, not at bin centres
  #  * legend has been removed
  #  * consistency bars have been added, calculated by a resampling technique
  #  * see Broecker (2007) http://dx.doi.org/10.1175/WAF993.1 for details
  #  * the bars are pointwise 2.5% ... 97.5% intervals around the hypothesis of reliability
  #  * dependency on package "verification" was removed
  #
  # Author: Stefan Siegert <s.siegert@exeter.ac.uk>
  #
  # based on previous version by Caio Coelho and the routine 
  # reliability.plot.default of the R-package `verification`
  #


  # sanity checks
  if (class(probs) == "data.frame") {
    probs <- c(as.matrix(probs))
  }
  if (class(obs) == "data.frame") {
    obs <- c(as.matrix(obs))
  }
  stopifnot(length(probs) == length(obs))
  stopifnot(nboot >= 0)
  stopifnot(all(probs >= 0), all(probs <= 1), all(obs %in% c(0,1)))
  stopifnot(length(cons.probs) == 2, all(cons.probs >= 0), all(cons.probs <= 1))

  # some definitions and corrections
  n <- length(obs)
  nboot <- floor(nboot)
  cons.probs <- sort(cons.probs)


  #############################################
  # reliability analysis
  #############################################
  # estimate refinement function
  if (length(bins) == 1) {
    nbins <- floor(bins)
    brx <- seq(0, 1, length.out=nbins+1) + 
         c(-.1, rep(0, nbins-1), .1)
  } else {
    nbins <- length(bins) - 1
    bins <- sort(bins)
    stopifnot(min(bins)<= 0 & max(bins) >= 1)
    brx <- bins
  }
  h <- hist(probs, breaks=brx, plot=FALSE)$counts        

  # estimate calibration function
  g <- hist(probs[obs==1], breaks=brx, plot=FALSE)$counts
  obar.i <- g / h 
  obar.i[ is.nan(obar.i) ] <- NA
  
  # calculate in-bin averages
  p.bins <- as.numeric(cut(probs, breaks=brx, include.lowest=TRUE))
  p.avgs <- sapply(seq(nbins), 
                   function(ii) mean(probs[p.bins == ii], na.rm=TRUE))
  p.avgs[ is.nan(p.avgs) ] <- NA



  #############################################
  # consistency resampling (broecker and smith 2007)
  #############################################
  if (nboot) {
    resamp.mat <- matrix(nrow=0, ncol=nbins)
    # the resampling function
    sample.rel.diag <- function(dummy=0) {
      p.hat <- sample(x=probs, size=n, replace=TRUE)
      x.hat <- rbinom(n=n, size=1, prob=p.hat)
      hh <- hist(p.hat, breaks=brx, plot=FALSE)$counts        
      gg <- hist(p.hat[x.hat==1], breaks=brx, plot=FALSE)$counts
      return(gg / hh)
    }
    l <- replicate(nboot, sample.rel.diag())
    resamp.mat <- t(l)
    cons.bars <- apply(resamp.mat, 2, 
                       function(z) quantile(z, cons.probs, na.rm=TRUE))
  } else {
    cons.bars <- matrix(NA, ncol=nbins, nrow=2)
  }


  #############################################
  # plot the reliability diagram
  #############################################
  if (plot) {
    # reliability plot
    old.par <- par(no.readonly = TRUE) 
    on.exit(par(old.par))
    plot(NULL, xlim = c(0,1), ylim = c(0,1),
       xlab= "Forecast probability",
       ylab="Observed relative frequency")
   if (attributes) {
       obs.clim<- sum(g)/sum(h)
       a   <- (1-obs.clim)/2 + obs.clim
       b   <- obs.clim / 2
       x.p <- c(obs.clim, obs.clim, 1, 1, 0, 0)
       y.p <- c(0, 1, 1, a, b, 0)
       polygon(x.p, y.p, col = "#e6e6e6")
       abline(h=obs.clim,lty=2)
       text( 0.9, obs.clim, "No resolution", pos = 3)
       text(0.9, obs.clim + (a-b)*(0.9 - obs.clim), "No skill", pos = 1,
     	   srt = atan( a - b )/(2*pi)*360 )
    }
    # consistency bars
    for (i in 1:length(p.avgs)) {
        lines(rep(p.avgs[i], 2), cons.bars[, i], col="#CCCCCC", lwd=6)
    }
    # reliability points and diagonal
    points(p.avgs, obar.i, col = "black", pch = 1, lwd=2, type="b")
    lines(c(0,1), c(0,1), lty=1)
    if (plot.refin) {
      # refinement histogram in lower corner
      pp<- par("plt")
      par("plt" = c(pp[2] - 0.2 , pp[2],  pp[3], pp[3]+ 0.2) )
      par(new = TRUE)
      barplot(h, axes = FALSE, axisnames = FALSE)
      axis(4)
      box() 
    }
  }

  #############################################
  # return data
  #############################################
  ret.df <- data.frame(p.avgs=p.avgs, cond.probs=obar.i, 
                       cbar.lo=cons.bars[1,], cbar.hi=cons.bars[2,])
  return(ret.df)
}


