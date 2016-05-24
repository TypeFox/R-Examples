"pfactor.bernstein" <-
function(para, x=NULL, n=NULL,
            bern.control=NULL,
            poly.type=c("Bernstein", "Kantorovich"),
            stat.type =c("Mean", "Median"),
            fix.lower=NULL, fix.upper=NULL,
            lmr.dist=NULL, lmr.n=c("3","4", "5"),
            nsim=500, plot.em=TRUE, pfactors=NULL,
            p.lo=.Machine$double.eps, p.hi=1) {

    if(! is.null(bern.control)) {
         poly.type <- bern.control$poly.type
         stat.type <- bern.control$stat.type
         fix.lower <- bern.control$fix.lower
         fix.upper <- bern.control$fix.upper
    }

   if(is.null(para)) {
      warning("para is NULL, returning NA")
      return(NA)
   }
   if(p.hi > 1) {
      warning("The upper threshold of the pfactors is too high regardless of the real sequence given or default")
      return(NA)
   }
   if(p.lo < .Machine$double.eps) {
      warning("The lower threshold of the pfactors is too small regardless of the real sequence given or default")
      return(NA)
   }

   poly.type  <- match.arg(poly.type)
   stat.type  <- match.arg(stat.type )
   lmr.n      <- as.integer(match.arg(lmr.n))
   if(stat.type == "Median") {
     stat.func <- function(...) { return(median(...)) }
   } else {
     stat.func <- function(...) { return(mean(...))   }
   }
   lmr <- NULL
   if(! is.null(x)) {
      lmr <- lmoms(x, nmom=5)
      n <- length(x)
   } else {
      lmr <- par2lmom(para)
      if(length(lmr$L1) == 1) lmr <- lmorph(lmr)
      if(is.null(n)) {
          warning("Data not provided and a sample length was not provided in lieu, returning NA")
          return(NA)
      }
   }
   if(is.null(lmr.dist)) {
      lmr.dist <- lmr$ratios[lmr.n]
   }
   mom.type <- "<NA>"
   mom.text <- "<NA>"
   if(lmr.n == 3) {
      mom.text <- "(Tau3 - Tau3smooth)"
      mom.type <- "Tau3"
   } else if(lmr.n == 4) {
      mom.text  <- "(Tau4 - Tau4smooth)"
      mom.type <- "Tau4"
   } else {
      mom.text  <- "(Tau5 - Tau5smooth)"
      mom.type <- "Tau5"
   }
   lmr.lab <- paste(c(stat.type," error statistic ",mom.text),sep="",collapse="");

   if(is.null(pfactors)) {
      ps1  <- 10^(seq(log10(0.01), log10(0.20), by=0.07)) # p-factors to explore
      ps2  <- 10^(seq(log10(0.20), log10(0.80), by=0.04)) # p-factors to explore
      ps   <- c(ps1, ps2); # merge the two vectors
   } else {
      ps <- pfactors
   }
   ps <- ps[ps > p.lo & ps <= p.hi]

   err      <- rep(NA, nsim)
   err.stat <- rep(NA, length(ps))
   if(length(ps) > 1) message("This function can be extremely CPU intensive depending on the nsim for the ",
                              length(ps)," p-factors to test:");
   for(j in 1:length(ps)) {
      if(length(ps) > 1) message("-",j, appendLF=FALSE)
      p <- ps[j] # set the p-factor
      for(i in 1:nsim) {
         X <- rlmomco(n, para) # draw a random sample of the distribution
         # Compute the L-moments by Bernstein or Kantorovich polynomials
         lmr <- lmoms.bernstein(X, poly.type=poly.type,
                                   bound.type="Carv",
                                   fix.lower=fix.lower, fix.upper=fix.upper, p=p)
         err[i] <- (lmr.dist - lmr$ratios[lmr.n]) # The error
      } # END the simulation loop
      err.tmp <- stat.func(err, na.rm=TRUE)
      if(length(ps) == 1) return(err.tmp)
      err.stat[j] <- err.tmp
   } # END the p-factor loop
   if(length(ps) > 1) message("-DONE")

   lowess   <- lowess(log10(ps), err.stat) # smooth the computed errors
   crossing <- 10^(approx(lowess$y, lowess$x, xout=0)$y) # zero crossing

   if(plot.em) {
      outer <- max(abs(err.stat), abs(lowess$y))
      xlim <- c(min(ps), max(ps)); ylim <- c(-outer,outer)
      plot(ps, err.stat, type="n", log="x", xlab="p-factor", lmr.lab,
           xlim=xlim, ylim=ylim)
      lines(xlim, c(0,0), lty=2) # draw the origin
      lines(ps, err.stat)        # draw the actual simulation run
      lines(ps, lowess$y, col=2) # draw the smooth to the simulation run
      lines(c(crossing,crossing), ylim, lty=2)
   }

   zz <- list(pfactor=crossing,
              poly.type=poly.type, stat.type=stat.type, lmom.type=mom.type,
              fix.lower=fix.lower, fix.upper=fix.upper,
              source="pfactor.bernstein",
              ps=ps, err.stat=err.stat, err.smooth=lowess$y)
   return(zz)
}


