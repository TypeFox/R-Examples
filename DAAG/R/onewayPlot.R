"onewayPlot" <-
function (obj, trtnam="trt", axisht = 6, xlim = NULL, xlab = NULL, 
            lsdht = 1.5, hsdht = 0.5, textht = axisht - 2.5, oma = rep(1, 
                                                               4),
            angle = 80, alpha = 0.05) 
{
  if (prod(par()$mfrow) == 1) {
    opar <- par(mar = rep(0, 4), oma = oma, xpd = TRUE)
    on.exit(par(mar = opar$mar, oma = opar$oma, xpd = FALSE))
  }
  else {
    par(xpd = TRUE)
    on.exit(par(xpd = FALSE))
  }
  term.labels <- attr(summary(obj), "term.labels")
  if(is.null(trtnam))trtnam <- term.labels[1]
  blist <- summary.lm(obj)$coef
  nam <- rownames(blist)
  introw <- match("(Intercept)",nam, nomatch=0)
  if(introw>0)const <- blist[introw,1] else const <- 0
  nch <- nchar(trtnam)
  facrows <- substring(nam, 1, nch)%in%trtnam
  est <- c(0, blist[facrows, 1]) + const
  sedB <- blist[facrows, 2]
  sed.min <- min(sedB)
  sed.max <- max(sedB)
  sed.rms <- sqrt(mean(c(sed.min, sed.max)^2))
  if (sed.max - sed.min > 0.1 * sed.rms) {
    show.sed <- FALSE
    cat("\nDesign is unbalanced.  SEDs depend on the treatments compared.\n")
  }
  df <- obj$df
  talpha <- qt(1 - alpha/2, df)
  lsd <- talpha * sed.rms
  tukey <- qtukey(1 - alpha, nmeans = length(est), df)/sqrt(2)
  hsd <- tukey * sed.rms
  est.min <- min(est)
  est.max <- max(est)
  xtratex <- 1+1.3*strwidth("Tukey HSD", "figure")
  if (is.null(xlim)) {
    xlim <- range(est)
    xlim <- xlim + c(-0.05, 0.05) * diff(xlim)
  } 
  if (est[1] + hsd <= est.max) {
    hsdlim <- c(est[1], est[1] + hsd)
    lsdlim <- c(est[1], est[1] + lsd)
  }
  else if (est[1] - hsd >= est.min) {
    hsdlim <- c(est[1] - hsd, est[1])
    lsdlim <- c(est[1] - lsd, est[1])
    adjtxt <- 1
  }
  else if(xtratex*hsd<diff(xlim)){
    hsdlim <- c(est.min, est.min + hsd)
    lsdlim <- c(est.min, est.min + lsd)
    adjtxt <- 1
  }
  else {
    xlim[2] <- xlim[1]+xtratex*hsd
    hsdlim <- c(est.min, est.min + hsd)
    lsdlim <- c(est.min, est.min + lsd)
    adjtxt <- 1
  }
  plot(xlim[1], 0, xlim = xlim, ylim = c(0, 1), type = "n", 
       axes = FALSE, xlab = "", ylab = "", mgp = c(2, 0.5, 0))
  axisht <- axisht - round(par()$mar[1])
  textht <- textht - round(par()$mar[1])
  lsdht <- lsdht - round(par()$mar[1])
  hsdht <- hsdht - round(par()$mar[1])
  chh <- par()$cxy[2] * par()$cex
  chw <- par()$cxy[1] * par()$cex
  lines(xlim, rep(axisht * chh, 2))
  axis(1, tck = 0.02, pos = axisht * chh, at = est, labels = FALSE)
  axis(1, pos = axisht * chh)
  if (!is.null(xlab)) 
    text(mean(par()$usr[1:2]), textht * chh, labels = xlab)
  trtlev <- obj$xlevels[[trtnam]]
  xpos <- bounce(est, d = chw)
  text(xpos, rep(axisht * chh, length(xpos)) + 0.85 * chh, 
       trtlev, srt = angle, adj = 0)
  adjtxt <- 0
  if (!is.null(lsdht)) {
    lines(lsdlim, rep(lsdht * chh, 2))
    text(lsdlim[2 - adjtxt] + (0.5 - adjtxt) * chw, lsdht * 
         chh, "LSD", adj = adjtxt)
  }
  if (!is.null(hsdht)) {
    lines(hsdlim, rep(hsdht * chh, 2))
    text(hsdlim[2 - adjtxt] + (0.5 - adjtxt) * chw, hsdht * 
         chh, "Tukey HSD", adj = adjtxt)
  }
  print(par()$mfg)
  invisible()
}

