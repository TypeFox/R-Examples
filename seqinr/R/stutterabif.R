stutterabif <- function(abifdata, chanel, poswild, 
  datapointbefore = 70, datapointafter = 20,
  datapointsigma = 3.5,
  chanel.names = c(1:4,105),
  DATA = paste("DATA", chanel.names[chanel], sep = "."),
  maxrfu = 1000,
  method = "monoH.FC",
  pms = 6,
  fig = FALSE){
  	
  xseq <- floor((poswild-datapointbefore):(poswild+datapointafter))
  yseq <- abifdata$Data[[DATA]][xseq]
  #
  # The baseline for the signal is taken as the most common value:
  #
  baseline <- baselineabif(yseq, maxrfu = maxrfu)
  #
  # Values below baseline are forced to baseline:
  #
  yseq[yseq < baseline] <- baseline
  #
  # Set baseline to zero:
  #
  yseq <- yseq - baseline
  #
  # First pass to fit a gaussian on wild allele:
  #
  ytheo.one <- function(p, x) p[1]*dnorm(x, p[2], p[3])
  obj.one <- function(p) sum(( yseq - ytheo.one(p, xseq))^2)
  guess.one <- numeric(3)
  guess.one[1] <- datapointsigma*sqrt(2*pi)*max(yseq)
  guess.one[2] <- poswild
  guess.one[3] <- datapointsigma
  suppressWarnings(nlmres.one <- nlm(obj.one, p = guess.one))
  #
  # Second pass to estimate stutter allele parameters. We do
  # not take into account wild peak allele data.
  #
  est <- nlmres.one$estimate
  censored <- floor(est[2] - 6*est[3]) # cut at mu - 6 sigma
  ytheo.two <- function(p, x) p[1]*dnorm(x, p[2], p[3])
  used <- ( xseq < censored )
  obj.two <- function(p) sum(( yseq[used] - ytheo.two(p, xseq[used]))^2)
  guess.two <- numeric(3)
  guess.two[1] <- datapointsigma*sqrt(2*pi)*max(yseq[used])
  guess.two[2] <- xseq[which.max(yseq[used])]
  guess.two[3] <- datapointsigma
  suppressWarnings(nlmres.two <- nlm(obj.two, p = guess.two))
  est2 <- nlmres.two$estimate
  #
  # Fit a spline:
  # 
  spfun <- splinefun(xseq, yseq, method = method)
  xx <- seq(min(xseq), max(xseq), le = 500)
  yy <- spfun(xx)  
  #
  # Compute output result:
  #
  
  x1inf <- est2[2] - pms*est2[3]
  x1sup <- est2[2] + pms*est2[3]
  l1 <- optimize(spfun, interval = c(x1inf, x1sup), maximum = TRUE)$maximum
  h1 <- spfun(l1)
  s1 <- integrate(spfun, x1inf, x1sup)$value
  
  x2inf <- est[2] - pms*est[3]
  x2sup <- est[2] + pms*est[3]
  l2 <- optimize(spfun, interval = c(x2inf, x2sup), maximum = TRUE)$maximum
  h2 <- spfun(l2)
  s2 <- integrate(spfun, x2inf, x2sup)$value

  rh <- h1/h2
  rs <- s1/s2
  
  p <- list(p1 = est2[1], mu1 = est2[2], sd1 = est2[3],
            p2 = est[1], mu2 = est[2], sd2 = est[3],
            xseq = xseq, yseq = yseq)
  if(fig){
    plot(xseq, yseq,
      main = paste("rh = ", round(rh, 5), "rs =", round(rs, 5)),
      ylab = "RFU", las = 1,
      xlab = "Time in datapoint units")
    abline(h = 0, lty = 3)
    abline(v = c(l1, l2), lty=2, col = c("red", "blue"))
    abline(h = c(h1, h2), lty=2, col = c("red", "blue"))

    #lines(xx, yy, col = "red")
    #
    # Stutter:
    #
    xx <- seq(x1inf, x1sup, le = 500)
    yy <- spfun(xx)
    lines(xx, yy, col = "red", lwd = 2)
    #
    # Wild Allele:
    #
    xx <- seq(x2inf, x2sup, le = 500)
    yy <- spfun(xx)
    lines(xx, yy, col = "blue", lwd = 2)
    #
    # legend:
    #
    legend("topleft", inset = 0.01, legend = c("Stutter product",
      "Corresponding allele"), lty = 1, lwd = 2, col = c("red", "blue"),
      bg = "white", title = "Datapoints used for:")
  }
  return(list(rh = rh, rs = rs, h1 = h1, h2 = h2, s1 = s1, s2 = s2, p = p))
}

