funnel.limitmeta <- function(x,
                             ##
                             pch=21,
                             cex=1,
                             col="black",
                             bg="darkgray",
                             ##
                             lwd=1,
                             ##
                             pch.adjust=18,
                             cex.adjust=1.5,
                             col.adjust="gray",
                             bg.adjust="gray",
                             ##
                             line=TRUE,
                             xmin.line,
                             xmax.line,
                             lty.line=1,
                             lwd.line=lwd,
                             col.line="gray",
                             ##
                             shrunken=FALSE,
                             pch.shrunken=22,
                             cex.shrunken=1,
                             col.shrunken="black",
                             bg.shrunken="white",
                             ##
                             lty.connect=1,
                             lwd.connect=0.8,
                             col.connect="black",
                             ##
                             backtransf=x$backtransf,
                             ...){
  
  
  meta:::chkclass(x, "limitmeta")
  
  
  TE <- x$TE
  seTE <- x$seTE
  ##
  TE.limit <- x$TE.limit
  seTE.limit <- x$seTE.limit
  ##
  minTE <- min(TE, na.rm=TRUE)
  maxTE <- max(TE, na.rm=TRUE)
  x.incr <- (maxTE-minTE)/1000
  ##
  TE.adjust <- x$TE.adjust
  ##
  tau <- x$tau
  alpha.r <- x$alpha.r
  beta.r <- x$beta.r
  ##
  sm <- x$sm
  
  
  if (alpha.r < 0){
    if (missing(xmin.line))
      xmin.line <- minTE
    if (missing(xmax.line))
      xmax.line <- TE.adjust - x.incr
  }
  if (alpha.r > 0){
    if (missing(xmin.line))
      xmin.line <- TE.adjust + x.incr
    if (missing(xmax.line))
    xmax.line <- maxTE
  }
  
  
  if (backtransf & meta:::is.relative.effect(sm)){
    TE <- exp(TE)
    TE.limit <- exp(TE.limit)
    TE.adjust <- exp(TE.adjust)
  }
  
  
  ##
  ## Generate funnel plot
  ##
  funnel(x$x, pch=pch, cex=cex, col=col, bg=bg, lwd=lwd,
         backtransf=backtransf, ...)
  
  
  ##
  ## Add line for adjustment method beta0
  ##
  if (line){
    if (x$method.adjust=="beta0"){
      if (backtransf & meta:::is.relative.effect(sm)){
        curve(sqrt((log(x)-beta.r)^2 / alpha.r^2 - tau^2),
              from=exp(xmin.line), to=exp(xmax.line),
              lty=lty.line, col=col.line, lwd=lwd.line, add=TRUE)
      }
      else{
        curve(sqrt((x-beta.r)^2 / alpha.r^2 - tau^2),
              from=xmin.line, to=xmax.line,
              lty=lty.line, col=col.line, lwd=lwd.line, add=TRUE)
      }
    }
  }
  
  
  ##
  ## Add adjusted treatment effect
  ##
  points(TE.adjust, 0, pch=pch.adjust, cex=cex.adjust, col=col.adjust, bg=bg.adjust)
  
  
  ##
  ## Add lines
  ##
  if (shrunken)
    segments(TE, seTE, TE.limit, seTE.limit,
             lty=lty.connect, lwd=lwd.connect, col=col.connect)
  
  
  ##
  ## Plot studies again
  ##
  points(TE, seTE, pch=pch, cex=cex, col=col, bg=bg)
  
  
  ##
  ## Add shrunken estimates
  ##
  if (shrunken)
    points(TE.limit, seTE.limit,
           pch=pch.shrunken, cex=cex.shrunken, col=col.shrunken, bg=bg.shrunken)
  
  
  invisible(NULL)
}
