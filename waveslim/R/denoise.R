manual.thresh <- function(wc, max.level=4, value, hard=TRUE)
{
  wc.fine <- wc[["d1"]]
  factor <- median(abs(wc.fine)) / .6745

  wc.shrink <- wc

  if(hard) {
    # Hard thresholding
    for(i in names(wc)[1:max.level]) {
      wci <- wc[[i]]
      unithresh <- factor * value
      wc.shrink[[i]] <- wci * (abs(wci) > unithresh)
    }
  }
  else {
    # Soft thresholding
    for(i in names(wc)[1:max.level]) {
      wci <- wc[[i]]
      unithresh <- factor * value
      wc.shrink[[i]] <- sign(wci) * (abs(wci) - unithresh) * 
        (abs(wci) > unithresh)
    }
  }
  wc.shrink
}

universal.thresh <- function(wc, max.level=4, hard=TRUE)
{
  n <- length(idwt(wc))

  wc.fine <- wc[["d1"]]
  factor <- median(abs(wc.fine)) / .6745

  wc.shrink <- wc

  if(hard) {
    # Hard thresholding
    for(i in names(wc)[1:max.level]) {
      wci <- wc[[i]]
      unithresh <- factor * sqrt(2 * log(n))
      wc.shrink[[i]] <- wci * (abs(wci) > unithresh)
    }
  }
  else {
    # Soft thresholding
    for(i in names(wc)[1:max.level]) {
      wci <- wc[[i]]
      unithresh <- factor * sqrt(2 * log(n))
      wc.shrink[[i]] <- sign(wci) * (abs(wci) - unithresh) *
        (abs(wci) > unithresh)
    }
  }
  wc.shrink
}

universal.thresh.modwt <- function(wc, max.level=4, hard=TRUE)
{
  n <- length(wc[[1]])
  wc.fine <- wc[["d1"]]
  factor <- sqrt(2) * median(abs(wc.fine)) / .6745
  wc.shrink <- wc
  j <- 1
  if(hard) {
    ## Hard thresholding
    for(i in names(wc)[1:max.level]) {
      wci <- wc[[i]]
      unithresh <- factor * sqrt(2 * log(n)) / 2^(j/2)
      wc.shrink[[i]] <- wci * (abs(wci) > unithresh)
      j <- j+1
    }
  }
  else {
    ## Soft thresholding
    for(i in names(wc)[1:max.level]) {
      wci <- wc[[i]]
      unithresh <- factor * sqrt(2 * log(n)) / 2^(j/2)
      wc.shrink[[i]] <- sign(wci) * (abs(wci) - unithresh) *
        (abs(wci) > unithresh)
      j <- j+1
    }
  }
  wc.shrink
}

sure.thresh <- function(wc, max.level=4, hard=TRUE)
{
  wc.shrink <- wc

  sure <- function(t, x) {
    ax <- sort(abs(x))    
    num <- match(FALSE, ax <= t, nomatch = length(ax) + 1) - 1
    length(ax) - 2 * num + sum(pmin(ax, t)^2)
  }

  for(i in names(wc)[1:max.level]) {
    wci <- wc[[i]]
    ni <- length(wci)
    factor <- median(abs(wci)) / .6745
    xi <- wci / factor
    sxi <- sort(abs(xi))^2
    s <- cumsum(sxi) + ((ni - 1):0) * sxi
    risk <- (ni - (2 * (1:ni)) + s) / ni
    surethresh <- sqrt(sxi[order(risk)[1]])
    if(hard) {
      ## Hard thresholding
      wc.shrink[[i]] <- wci * (abs(xi) > surethresh)
    }
    else {
      ## Soft thresholding
      wc.shrink[[i]] <- sign(wci) * (abs(wci) - factor*surethresh) *
        (abs(xi) > surethresh)
    }
  }
  return(wc.shrink)
}

hybrid.thresh <- function(wc, max.level = 4, verbose = FALSE, seed = 0)
{
  shrinkit <- function(coeffs, thresh)
    sign(coeffs) * pmax(abs(coeffs) - thresh, 0)

  sure <- function(t, x) {
    ax <- sort(abs(x))    
    num <- match(FALSE, ax <= t, nomatch = length(ax) + 1) - 1
    length(ax) - 2 * num + sum(pmin(ax, t)^2)
  }

  wc.shrink <- wc
  n <- length(unlist(wc))
  nlev <- log(n + 1, 2) - 1
  i <- 1
  iloc <- 1
  while(i <= max.level) {
    ## Extract current level coefficients from all wavelet coefficients
    raw <- wc[[names(wc)[i]]]
    d <- length(raw)
    ## Test:  if the variance is small enough, just use threshold sqrt(2logd)
    if((sum(raw^2) - d)/d <= sqrt(i^3/2^i)) {
      if(verbose)
        cat(paste("At level ", i, " the threshhold is sqrt(2log(d)): ", 
                  sqrt(2 * log(d)), "\n", sep = ""))
      wc.shrink[[names(wc)[i]]] <- shrinkit(wc[[names(wc)[i]]], sqrt(2*log(d)))
    }
      else {
        ## Generate random subset
        if(length(seed) != 1)
          .Random.seed <- seed
        Iset <- sort(sample(d, d/2))
        rawI <- raw[Iset] / (median(abs(raw[Iset])) / .6745)
        rawIp <- raw[ - Iset] / (median(abs(raw[ - Iset])) / .6745)
        ggI <- sort(abs(rawI))
        ggIp <- sort(abs(rawIp))        
        ## Calculate SURE for all possible thresholds
        surevecI <- sapply(c(ggI[ggI < sqrt(2 * log(d))], 0,
                             sqrt(2 * log(d))), sure, ggI)
        surevecIp <- sapply(c(ggIp[ggI < sqrt(2 * log(d))], 0, 
                              sqrt(2 * log(d))), sure, ggIp)    
        ## Threshold that minimizes risk
        llI <- length(surevecI)
        llIp <- length(surevecIp)       
        ## The minimum occurs either at sqrt(2logd), 
        if(min(surevecI) == surevecI[llI])
          threshI <- sqrt(2 * log(d))
        else if(min(surevecI) == surevecI[llI - 1])
          threshI <- 0
        else threshI <- ggI[match(min(surevecI), surevecI)]     
        ## or at 0, 
        if(min(surevecIp) == surevecIp[llIp])
          threshIp <- sqrt(2 * log(d))
        else if(min(surevecIp) == surevecI[llIp - 1])
          threshIp <- 0
        else
          threshIp <- ggIp[match(min(surevecIp), surevecIp)]
        ## or at 0, 
        if(verbose) {
          cat(paste("At level ", i, ", threshold1 is ", threshI, "\n",
                    sep = ""))
          cat(paste("At level ", i, ", threshold2 is ", threshIp,
                    "\n", sep = ""))
        }
        ## Perform shrinking
        newI <- shrinkit(rawI, threshIp)
        newIp <- shrinkit(rawIp, threshI)
        new <- rep(0, d)
        new[Iset] <- newI
        new[ - Iset] <- newIp
        wc.shrink[[names(wc)[i]]] <- new
      }
    ## Otherwise, go through all this stuff
    iloc <- iloc + 2^i
    i <- i + 1
  }
  wc.shrink
}

da.thresh <- function(wc, alpha=.05, max.level=4, verbose=FALSE,
                      return.thresh=FALSE) {
  onebyone2 <- function(dat, alpha) {
    kolsmi.chi2 <- function(dat) {
      n <- length(dat)
      return(max(abs(cumsum(dat)-(1:n)*sum(dat)/n))/sqrt(2*n))
    }
    crit <- c(seq(0.28,1.49,by=.01), seq(1.50,2.48,by=.02))
    alph <- c(.999999,.999996,.999991,.999979,.999954,.999909,.999829,
              .999697,.999489,.999174,.998715,.998071,.997192,.996028,
              .994524,.992623,.990270,.987410,.983995,.979978,.975318,
              .969983,.963945,.957186,.949694,.941466,.932503,.922817,
              .912423,.901344,.889605,.877240,.864282,.850771,.836775,
              .822247,.807323,.792013,.776363,.760418,.744220,.727811,
              .711235,.694529,.677735,.660887,.644019,.627167,.610360,
              .593628,.576998,.560495,.544143,.527959,.511970,.496192,
              .480634,.465318,.450256,.435454,.420930,.406684,.392730,
              .379072,.365714,.352662,.339918,.327484,.315364,.303556,
              .292060,.280874,.270000,.259434,.249174,.239220,.229566,
              .220206,.211140,.202364,.193872,.185658,.177718,.170050,
              .162644,.155498,.148606,.141962,.135558,.129388,.123452,
              .117742,.112250,.106970,.101896,.097028,.092352,.087868,
              .083568,.079444,.075495,.071712,.068092,.064630,.061318,
              .058152,.055128,.052244,.049488,.046858,.044350,.041960,
              .039682,.037514,.035448,.033484,.031618,.029842,.028154,
              .026552,.025030,.023588,.022218,.019690,.017422,.015390,
              .013574,.011952,.010508,.009223,.008083,.007072,.006177,
              .005388,.004691,.004078,.003540,.003068,.002654,.002293,
              .001977,.001703,.001464,.001256,.001076,.000921,.000787,
              .000671,.000572,.000484,.000412,.000350,.000295,.000250,
              .000210,.000178,.000148,.000126,.000104,.000088,.000074,
              .000060,.000051,.000042,.000035,.000030,.000024,.000020,
              .000016,.000013,.000011,.000009)
    if(alpha < min(alph) || alpha > max(alph))
      stop("alpha =",alpha,"is out of range")
    ind <- match(TRUE, alpha > alph)
    critval <- crit[ind-1]+(alph[ind-1]-alpha)*(crit[ind]-crit[ind-1]) / 
      (alph[ind-1]-alph[ind])
    i <- length(dat)
    cc <- kolsmi.chi2(dat)
    while(cc[length(cc)] > critval && i > 1) {
      i <- i-1
      cc <- c(cc,kolsmi.chi2(dat[sort(order(dat)[1:i])]))
    }
    return(cc)
  }
  getthrda2 <- function(dat, alpha) {
    a <- onebyone2(dat, alpha)
    if(length(a) == length(dat))
      if(1 - pchisq(min(dat),1) < alpha)
        ggg <- 0
      else
        ggg <- sqrt(min(dat))
    else
      ggg <- sqrt(max(dat[sort(order(dat)[1:(length(dat)-length(a)+1)])]))
    return(ggg)
  }
  shrinkit <- function(coeffs, thresh)
    sign(coeffs) * pmax(abs(coeffs) - thresh, 0)

  if(alpha <= .000009 || alpha >= .999999)
    stop("alpha out of range")
  ans <- wc
  n <- length(unlist(wc))
  nlev <- log(n+1, 2)-1
  i <- 1
  iloc <- 1
  while(i <= max.level) {
    gg <- wc[[names(wc)[i]]]
    thresh <- getthrda2(gg^2,alpha)
    if(verbose)
      cat(paste("At level ",i,", the threshold is ",thresh, "\n",sep=""))
    if(return.thresh)
      if(i == nlev)
        rt <- thresh
      else
        rt <- c(thresh, rt)
    else
      ans[[names(wc)[i]]] <- shrinkit(wc[[names(wc)[i]]], thresh)
    iloc <- iloc + 2^i
    i <- i+1
  }
  if(return.thresh)
    return(rt)
  else
    return(ans)
}
