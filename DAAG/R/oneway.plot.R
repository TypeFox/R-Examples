"oneway.plot" <-
function (obj, axisht = 6, xlim=NULL, xlab=NULL,
              lsdht = 1.5, hsdht = 0.5, textht=axisht-2.5,
              oma=rep(1,4), angle=80, alpha = 0.05)
{
    if(prod(par()$mfrow)==1){
          opar <- par(mar = rep(0, 4), oma = oma, xpd=TRUE)
          on.exit(par(mar=opar$mar,oma=opar$oma,xpd=FALSE))
      } else
      {par(xpd=TRUE)
       on.exit(par(xpd=FALSE))
   }
    b <- coef(obj)
    est <- b[1] + c(0, b[-1])
    sed <- summary.lm(obj)$coef[-1, 2]
    sed.min <- min(sed)
    sed.max <- max(sed)
    sed.rms <- sqrt(mean(sed.min^2 + sed.max^2))
    if (sed.max - sed.min > 0.1 * sed.rms) {
        show.sed <- FALSE
        cat("\nDesign is unbalanced.  SEDs depend on the treatments compared.\n")
    }
    if(is.null(xlim)){xlim <- range(est)
                      xlim <- xlim + c(-0.05, 0.05) * diff(xlim)
                  }
    plot(xlim[1], 0, xlim = xlim, ylim = c(0, 1), type="n", axes=FALSE,
         xlab="", ylab="", mgp=c(2,0.5,0))
    axisht <- axisht-round(par()$mar[1])
    textht <- textht-round(par()$mar[1])
    lsdht <- lsdht-round(par()$mar[1])
    hsdht <- hsdht-round(par()$mar[1])
    chh <- par()$cxy[2]*par()$cex
    chw <- par()$cxy[1]*par()$cex
    lines(xlim, rep(axisht*chh, 2))
    axis(1, tck = 0.02, pos = axisht*chh, at = est, labels = FALSE)
    axis(1, pos = axisht*chh)
    if(!is.null(xlab)) text(mean(par()$usr[1:2]), textht*chh,
                            labels=xlab)
    trtnam <- all.names(obj$call$formula)[3]
    trtlev <- obj$xlevels[[trtnam]]
    xpos <- bounce(est, d = chw)
    text(xpos, rep(axisht*chh, length(xpos)) + 0.85 * chh, trtlev,
         srt =  angle, adj = 0)
    df <- obj$df
    talpha <- qt(1 - alpha/2, df)
    lsd <- talpha * sed.rms
    tukey <- qtukey(1 - alpha, nmeans = length(est), df)/sqrt(2)
    hsd <- tukey * sed.rms
    est.min <- min(est)
    est.max <- max(est)
    adjtxt <- 0
    if (est[1] + hsd <= est.max) {
        hsdlim <- c(est[1], est[1] + hsd)
        lsdlim <- c(est[1], est[1] + lsd)
    }
    else if (est[1] - hsd >= est.min) {
        hsdlim <- c(est[1] - hsd, est[1])
        lsdlim <- c(est[1] - lsd, est[1])
        adjtxt <- 1
    }
    else {
        hsdlim <- c(est[1], est[1] + hsd)
        lsdlim <- c(est[1], est[1] + lsd)
    }
    if(!is.null(lsdht)){
        lines(lsdlim, rep(lsdht*chh, 2))
        text(lsdlim[2 - adjtxt] + (0.5 - adjtxt) * chw, lsdht*chh, "LSD",
             adj = adjtxt)}
    if(!is.null(hsdht)){
        lines(hsdlim, rep(hsdht*chh, 2))
        text(hsdlim[2 - adjtxt] + (0.5 - adjtxt) * chw, hsdht*chh, "Tukey HSD",
             adj = adjtxt)}
    print(par()$mfg)
    invisible()
}
