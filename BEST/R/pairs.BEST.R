pairs.BEST <-
function(x, nPtToPlot = 1000, col="skyblue", ...) {
    # Plot the parameters pairwise, to see correlations

  # Sanity checks:
  if(!inherits(x, "data.frame"))
    stop("x is not a valid BEST object")
  if(ncol(x) == 3 && all(colnames(x) == c("mu","nu","sigma"))) {
    oneGrp <- TRUE
  } else if (ncol(x) == 5 && all(colnames(x) == c("mu1", "mu2","nu","sigma1","sigma2"))) {
    oneGrp <- FALSE
  } else {
    stop("x is not a valid BEST object")
  }

    nuCol <- which(colnames(x) == "nu")
    mcmcChain <- cbind(x[, -nuCol], log10(x$nu))
    #plotIdx = floor(seq(1, nrow(mcmcChain),by=nrow(mcmcChain)/nPtToPlot)) #TODO Use length.out
    plotIdx = floor(seq(1, nrow(mcmcChain), length.out=nPtToPlot)) #TODO Use length.out

    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt) # cex.cor is now cruft?
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    if(oneGrp) {
          labels <- c( expression(mu), 
                     expression(sigma), 
                     expression(log10(nu)) )
    } else {
          labels <- c( expression(mu[1]) , expression(mu[2]) , 
                     expression(sigma[1]) , expression(sigma[2]) , 
                     expression(log10(nu)) )
    }
    pairs( mcmcChain[plotIdx, ] , labels=labels,
            lower.panel=panel.cor, col=col,  ... )
  }
