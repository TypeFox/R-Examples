unstrip <- function(x)
{
   dd <- dim(x)
   y <- x
   if(length(dd)==2){
      dd2 <- dd[2]
      if(dd2==1) y<- c(x[,1])
      if(dd2==2) y<- cbind(c(x[,1]),c(x[,2]))
      if(dd2>2) y<- cbind(c(x[,1]),c(x[,2]),c(x[,3]))
      if(dd2>3)for(i in 4:dd2) y <- cbind(y,c(x[,i]))
      y
   }
   if(length(dd)==1 || length(dd)==0){
      y <- c(unlist(c(unlist(x))))
      names(y) <- NULL
   }
   y
}
   
hare <- function(data, delta, cov, penalty, maxdim, exclude,
   include, prophaz = FALSE, additive = FALSE, linear, fit, silent = TRUE)
{
# get the parameters from the C-program
   call <- match.call()
        if(!missing(data)) data <- unstrip(data)
        if(!missing(delta)) delta <- unstrip(delta)
        if(!missing(cov)) cov <- unstrip(cov)
        if(!missing(exclude)) exclude <- unstrip(exclude)
        if(!missing(include)) include <- unstrip(include)
   MAXKNOTS <- -3
   MAXSPACE <- -3
   z <- .C("sharex",
      mk = as.integer(MAXKNOTS),
      ms = as.integer(MAXSPACE),
      PACKAGE = "polspline")
   MAXKNOTS <- z$mk
   MAXSPACE <- z$ms   # 
# a few elementary data checks
   if(missing(data))
      stop("there has to be data")
   if(length(data) < 25)
      stop("not enough data")
   if(min(data) < 0)
      stop("negative data")
   if(missing(delta))
      delta <- data - data + 1
   if(length(data) != length(delta))
      stop("data and delta have different length")
   dd <- abs(delta - 0.5)
   if(min(dd) < 0.5 || max(dd) > 0.5)
      stop("delta not all 0 or 1")
   ndata <- length(data)   #
# dealing with the covariates, sorting the cases 
# first if there are no covariates
   if(missing(cov)) {
      ncov <- 0
      cov <- 0
      iia <- order(data)
      delta <- delta[iia]
      data <- data[iia]
   }
   else {
      if(length(cov) == ndata)
         cov <- matrix(cov, ncol = 1, nrow = ndata)
      if(length(cov[, 1]) != ndata)
         stop("covariates not ndata * ncov matrix")
      ncov <- length(cov[1,  ])
      cov <- cbind(cov, 1)
      y <- cbind(data, cov)
      keys <- 1:ndata
      for(i in (ncov + 2):1)
         keys <- keys[sort.list(y[keys, i])]
      data <- data[keys]
      delta <- delta[keys]
      cov <- cov[keys, 1:ncov]
   }
   if(additive) {
      if(!missing(exclude))
         stop("cannot have exclude and additive")
      if(!missing(include))
         stop("cannot have include and additive")
      prophaz <- FALSE
      include <- c(0, 0)
   }
   if(missing(exclude) + missing(include) == 0)
      stop("only 1 from exclude and include allowed")
   vexclude <- 0   #
# using exclude
   if(missing(exclude) == FALSE) {
      if(length(exclude) == 2)
         exclude <- matrix(exclude, ncol = 2, nrow = 1)
      if(length(exclude[1,  ]) != 2)
         stop("exclude has wrong shape")
      if(min(exclude) < 0 || max(exclude) > ncov)
         stop("exclude has wrong values")
      vexclude <- as.vector(t(exclude))
      vexclude <- c(length(vexclude)/2, vexclude)   #
# proportional hazards model
      if(prophaz && ncov > 0) {
         vexclude <- c(vexclude, as.vector(rbind(1:ncov, 0)))
         vexclude[1] <- vexclude[1] + ncov
      }
   }
#
# using include
   if(missing(include) == FALSE || additive) {
      if(length(include) == 2)
         include <- matrix(include, ncol = 2, nrow = 1)
      if(length(include[1,  ]) != 2)
         stop("include has wrong shape")
      if(min(include) < 0 || max(include) > ncov)
         stop("include has wrong values")
      include <- t(apply(include, 1, sort))
      if(length(include) == 2)
         include <- matrix(include, ncol = 2, nrow = 1)
      if(prophaz)
         include <- include[include[, 1] > 0,  ]
      vexclude <- as.vector(t(include))
      vexclude <- c( - length(vexclude)/2, vexclude)
   }
#
# using proprtional hazards
   if(missing(include) && missing(exclude) && prophaz && ncov > 0) 
         vexclude <- c(ncov, as.vector(rbind(1:ncov, 0)))   
   # set parameters
   mindist <- 5
   if(missing(penalty))
      penalty <- log(ndata)
   if(missing(maxdim)) {
      maxdim <- floor(6 * (ndata)^0.2)
      if(maxdim > MAXSPACE - 1)
         maxdim <- MAXSPACE - 1
      maxdim <-  - maxdim
   }
   if(maxdim > MAXSPACE - 1) {
      maxdim <- MAXSPACE - 1
      print(paste("maximum dimension reduced to", maxdim))
   }
   lins <- rep(0, MAXSPACE)
   if(!missing(linear)) {
      linear[linear <= 0] <- ncov + 1
      linear[linear > ncov + 1] <- ncov + 1
      lins[linear] <- 1
   }
   if(additive)
      vexclude <- c(-1, 0, 0)   # do it
   fitter <- 0
   bbtt <- matrix(0, ncol = 6, nrow = abs(maxdim))
   cckk <- matrix(0, ncol = (MAXKNOTS + 1), nrow = (ncov + 1))
   if(!missing(fit)) {
                if(class(fit)!="hare")
         stop("fit is not a hare object")
      fitter <- fit$ndim
      if(fit$ncov != ncov)
         stop("ncov and fit's ncov are different")
      bbtt[1:fit$ndim,  ] <- fit$fcts
      bbtt <- as.vector(t(bbtt))
      bbtt[is.na(bbtt)] <- -1
      a1 <- length(fit$knots[1,  ])
      cckk[, 1:a1] <- fit$knots
      cckk <- as.vector(cckk)
      cckk[is.na(cckk)] <- -1
   }
   z <- .C("share",
      as.integer(ncov),
      ndim = as.integer(ndata),
      as.double(data),
      as.integer(delta),
      as.double(cov),
      as.double(penalty),
      as.integer(mindist),
      as.integer(maxdim),
      bbtt = as.double(bbtt),
      cckk = as.double(cckk),
      as.integer(vexclude),
      as.integer(lins),
      as.integer(silent),
      logl = as.double(rep(0, MAXSPACE)),
      as.integer(fitter),
      ad = as.integer(rep(0, MAXSPACE)),
      as.integer(0),   # 
      PACKAGE = "polspline")
# organize bbtt and cckk
   maxdim <- abs(maxdim)
   z$bbtt <- matrix(z$bbtt, nrow = maxdim, ncol = 6, byrow = TRUE)[1:z$ndim,  
      ]
   z$cckk <- matrix(z$cckk, nrow = ncov + 1, ncol = MAXKNOTS + 1, byrow = 
      TRUE)
   z$cckk <- z$cckk[, 1:(1 + max(z$cckk[, 1]))]
   z$cckk <- matrix(z$cckk, nrow = ncov + 1)
   l1 <- max(z$cckk[, 1])
   for(i in 1:(ncov + 1))
      if(z$cckk[i, 1] != l1) z$cckk[i, (z$cckk[i, 1] + 2):(l1 + 1)] <- 
            NA
   if(l1 > 0 && ncov > 0)
      dimnames(z$cckk) <- list(c("T", 1:ncov), c("K", 1:l1))
   if(l1 > 0 && ncov == 0)
      dimnames(z$cckk) <- list(c("T"), c("K", 1:l1))
   if(l1 == 0 && ncov > 0)
      dimnames(z$cckk) <- list(c("T", 1:ncov), "K")
   if(l1 == 0 && ncov == 0)
      dimnames(z$cckk) <- list(c("T"), "K")
   l1 <- max((1:MAXSPACE)[z$logl > -1e+100])
   z$bbtt <- matrix(z$bbtt, ncol = 6)
   dimnames(z$bbtt) <- list(1:(z$ndim), c("dim1", "knot1", "dim2", "knot2",
      "beta", "SE"))
   z$bbtt[z$bbtt[, 3] == -1, 3:4] <- NA
   if(is.na(l1)) {
      z$logl <- z$logl[1]
      z$ad <- z$ad[1]
   }
   else {
      z$logl <- z$logl[1:l1]
      z$ad <- z$ad[1:l1]
   }
   z$ad[z$logl < -1e+100] <- NA
   z$logl[z$logl < -1e+100] <- NA
   z$logl <- cbind(z$logl, z$ad)
   dimnames(z$logl) <- list(NULL, c("log-lik", "A/D"))
   ranges <- NA
   if(ncov == 1)
      ranges <- matrix(range(cov), ncol = 1, nrow = 2)
   if(ncov > 1)
      ranges <- apply(cov, 2, range)   # done
   fit <- list(call = call, ncov = ncov, ndim = z$ndim, fcts = z$bbtt, knots 
       = z$cckk, penalty = penalty, max = max(data), ranges = ranges,
      logl = z$logl, sample = ndata)
        class(fit) <- "hare"
        fit
}
plot.hare <- function(x, cov, n = 100, which = 0, what = "d", time, add = FALSE,
   xlim, xlab, ylab, type, ...)
{
    if(class(x)!="hare")
       stop("x is not a hare object")
    if(!missing(cov))cov <- unstrip(cov)
    if(!missing(time))time <- unstrip(time)
      fit <- x
   nocov <- 0
   if(fit$ncov == 0)
      nocov <- 1
   else {
           if(length(cov) != fit$ncov)
      stop("covariates are wrong")
        }
   if(which == 0) {
      if(missing(xlim)) {
         if(nocov == 0) {
            u1 <- qhare(0.01, cov, fit)
            u2 <- qhare(0.99, cov, fit)
         }
         else {
            u1 <- qhare(0.01, fit = fit)
            u2 <- qhare(0.99, fit = fit)
         }
         u3 <- 1.1 * u1 - 0.1 * u2
         u2 <- min(u2, fit$max)
         u4 <- 1.1 * u2 - 0.1 * u1
         if(u3 < 0)
            u3 <- 0
         else if(u4/u3 > 5)
            u3 <- 0
      }
      else {
         u3 <- xlim[1]
         u4 <- xlim[2]
      }
      xx <- (0:(n - 1))/(n - 1) * (u4 - u3) + u3
      if(fit$ncov > 0)
         yy <- cov
   }
   else {
      if(which < 0 || which > fit$ncov)
         stop("which is wrong")
      if(missing(time))
         stop("time is missing")
      if(missing(xlim)) {
         u3 <- fit$ranges[1, which]
         u4 <- fit$ranges[2, which]
      }
      else {
         u3 <- xlim[1]
         u4 <- xlim[2]
      }
      xx <- (0:(n - 1))/(n - 1) * (u4 - u3) + u3
      yy <- matrix(cov, ncol = fit$ncov, nrow = n, byrow = TRUE)
      yy[, which] <- xx
      xx <- time
   }
   iwhat <- 0
   if(what == "d" || what == "D")
      iwhat <- 3
   if(what == "h" || what == "H")
      iwhat <- 2
   if(nocov == 0)
      yy <- xhare(iwhat, xx, yy, fit)
   else yy <- xhare(iwhat, xx, arg4 = fit)
   if(what == "s" || what == "S")
      yy <- 1 - yy
   if(missing(xlab))
      xlab <- ""
   if(missing(ylab))
      ylab <- ""
   if(missing(type))
      type <- "l"
   xx <- (0:(n - 1))/(n - 1) * (u4 - u3) + u3
   if(!add)
      plot(xx, yy, xlab = xlab, ylab = ylab, type = type, ...)
   else lines(xx, yy, type = type, ...)
}
print.hare <- function(x,...)
{
      summary.hare(x)
}
summary.hare <- function(object,...)
{
    if(class(object)!="hare")
       stop("object is not a hare object")
             fit <- object
   s3 <- as.vector(t(fit$logl))
   s3[is.na(s3)] <- 0
   s1 <- as.vector(t(fit$fcts))
   s2 <- as.vector(fit$knots)
   s1[is.na(s1)] <- -1
   s2[is.na(s2)] <- -1
   .C("ssumm",
      as.double(fit$penalty),
      as.integer(fit$sample),
      as.double(s3),
      as.integer(length(s3)/2),
      as.double(s2),
      as.double(s1),
      as.integer(fit$ndim),
      as.integer(fit$ncov),
      PACKAGE = "polspline")
   invisible()
}
dhare <- function(q,cov,fit)
{
    if(class(fit)!="hare")
       stop("fit is not a hare object")
   xhare(3, q,cov,fit)
}
hhare <- function(q,cov,fit)
{
    if(class(fit)!="hare")
       stop("fit is not a hare object")
   xhare(2, q,cov,fit)
}
phare <- function(q,cov,fit)
{
    if(class(fit)!="hare")
       stop("fit is not a hare object")
   xhare(0, q,cov,fit)
}
qhare <- function(p,cov,fit)
{
    if(class(fit)!="hare")
       stop("fit is not a hare object")
   xhare(1, p,cov,fit)
}
rhare <- function(n, cov,fit)
{
    if(class(fit)!="hare")
       stop("fit is not a hare object")
   xhare(1, runif(n), cov,fit)
}
xhare <- function(arg1,arg2,arg3,arg4)
{
# mainly messing with the covariates 
    iwhat <- arg1
    if(!missing(arg2))arg2 <- unstrip(arg2)
    if(!missing(arg3))arg3 <- unstrip(arg3)
    q <- arg2
    cov <- arg3
    fit <- arg4
    if(class(fit)!="hare")
       stop("fit is not a hare object")
        zz <- 0
   if(missing(arg4)) {
                zz <- 7
      fit <- cov
      if(is.null(fit$ncov))
         stop("fit missing")
   }
   if(fit$ncov == 0) {
      if(!missing(arg3) && zz==0)
         stop("there should be no covariates")
      else cov <- 0
   }
   else {
      if(is.matrix(cov) == FALSE)
         cov <- matrix(cov, ncol = fit$ncov)
      nd <- length(cov[, 1])
      nc <- length(cov[1,  ])
      nq <- length(q)
      if(nc != fit$ncov)
         stop("not the right number of covariates")
      if(nd != 1 && nq != 1 && nd != nq)
         stop("no matching number of cases")
      if(nq == 1)
         q <- rep(q, nd)
      if(nd == 1 && nq != 1)
         cov <- matrix(cov, nrow = nq, ncol = nc, byrow = TRUE)
   }
   fit$fcts <- as.vector(t(fit$fcts))
   fit$fcts[is.na(fit$fcts)] <- -1
   fit$knots <- as.vector(fit$knots)
   fit$knots[is.na(fit$knots)] <- 0
   z <- .C("sphare",
      as.integer(fit$ncov),
      as.integer(fit$ndim),
      as.integer(length(q)),
      as.double(cov),
      as.integer(iwhat),
      q = as.double(q),
      as.double(fit$knots),
      as.double(fit$fcts),
      PACKAGE = "polspline")
   z$q
}
heft <- function(data, delta, penalty, knots, leftlin, shift,
   leftlog, rightlog, maxknots, mindist, silent = TRUE)
{
   call <- match.call()
        if(!missing(data))data <- unstrip(data)
        if(!missing(delta))delta <- unstrip(delta)
        if(!missing(knots))knots<- unstrip(knots)
   if(missing(leftlin))leftlin<-2
   leftlin<-leftlin*1
   nx <- -1
   z <- .C("sheftx",
      z = as.integer(nx),
      PACKAGE = "polspline")
   lgth <- z$z
   lgth <- 40
   if(missing(mindist))
      mindist <- 5
   if(mindist < 2) {
      warning("mindist reset to 2")
      mindist <- 2
   }
   if(missing(delta))
      delta <- data - data + 1
   if(length(data) != length(delta))
      stop("data and delta have different length")
   if(min(data) < 0)
      stop("negative data")
   if(min(data) == 0) {
      if(!missing(leftlog)){
      if(leftlog != 0) 
         stop("** hard-zeros, leftlog has to be 0 **")
      }
      else{
         leftlog <- 0
         warning("*** hard zeros: leftlog set to 0 ***")
      }
      if(leftlin==2){
         warning("*** hard zeros: leftlin set to TRUE ***")
         leftlin <- 1
       }
   }
   leftlin <- (leftlin==1)
   dd <- abs(delta - 0.5)
   if(min(dd) < 0.5 || max(dd) > 0.5)
      stop("delta not all 0 or 1")
   delta <- delta[order(data)]
   data <- sort(data)
   nx <- length(data)
   if(!missing(maxknots) && !missing(knots) && maxknots < length(knots))
      stop("maxknots is smaller than length(knots)")
   if(missing(maxknots))
      maxknots <- 0
   if(missing(penalty))
      penalty <- log(nx)
   if(maxknots > lgth - 5) {
      maxknots <- lgth - 5
      warning(paste("maxknots reduced to", maxknots))
   }
   if(!missing(shift))
      if(shift <=  - min(data))
         stop("shift too small")
   if(missing(shift))
      shift <- quantile(data[delta==1], 0.75)
   nknots <- 0
   iauto <- 0
   if(!missing(knots)) {
      nknots <- length(knots)
      if(nknots > lgth - 5)
         stop(paste("nknots can be at most", lgth - 5))
      iauto <- 2
      uu <- knots[2:nknots] - knots[1:(nknots - 1)]
      if(min(uu) < 0)
         stop("knots not in sequence")
      if(knots[1] < 0)
         stop("knot 1 is negative")
      knots <- c(knots, rep(0, lgth - nknots))
   }
   if(iauto < 2)
      knots <- rep(0, lgth)
   error <- c(1, rep(0, 20))
   if(silent != TRUE)
      error[7] <- 37
   tails <- c(0, 0, 0, 0, 1)
   if(!missing(leftlog) || min(data) == 0) {
      if(leftlog <= -1)
         stop("leftlog should be smaller than -1")
      tails[1] <- 1
      tails[2] <- leftlog
   }
   if(!missing(rightlog)) {
      if(rightlog < -1)
         stop("rightlog should be at least -1")
      tails[3] <- 1
      tails[4] <- rightlog
   }
   if(leftlin)
      tails[5] <- 0
   z <- .C("sheft",
      as.integer(nx),
      as.double(data),
      as.integer(delta),
      nk = as.integer(nknots),
      knots = as.double(knots),
      as.double(penalty),
      tails = as.double(tails),
      as.integer(iauto),
      logl = as.double(rep(0, lgth)),
      theta = as.double(rep(0, lgth)),
      iknots = as.integer(rep(0, lgth)),
      error = as.integer(error),
      as.double(shift),
      as.integer(maxknots),
      ad = as.integer(rep(0, lgth)),
      as.integer(mindist),
      PACKAGE = "polspline")
   error <- z$error
   if(z$nk < -100) error[2] <- 1
   z$logl[abs(z$logl) < 1e-100] <- 0
   z$logl[z$ad == 2] <- 0
   if(error[2] == 0){
      fit <-list(call = call, knots = z$knots[1:z$nk], logl = z$logl[2:(z$nk
         + 1)], thetak = z$theta[1:z$nk], thetap = z$theta[z$nk +
         (1:2)], thetal = z$theta[z$nk + (3:4)], penalty = 
         penalty, shift = shift, sample = length(data), logse = 
         z$tails[c(2, 4)], max = max(data), adddel = z$ad[2:(z$nk
          + 1)])
           class(fit) <- "heft"
           fit
        }
   else {
      print("sorry......")
      invisible()
   }
}
plot.heft <- function(x, n = 100, what = "d", add = FALSE, xlim, xlab, ylab, type,
    ...)
{
    if(class(x)!="heft")
       stop("x is not a heft object")
        fit <- x
   if(missing(xlim)) {
      u2 <- min(qheft(0.99, fit), fit$max)
      u3 <- 0
      u4 <- 1.1 * u2
      xlim <- c(u3, u4)
   }
   u3 <- xlim[1]
   u4 <- xlim[2]
   xx <- (0:(n - 1))/(n - 1) * (u4 - u3) + u3
   if(u3 == 0)
      xx <- (1:n)/n * u4
   yy <- c(-10, -10)
   if(what == "d" || what == "D")
      yy <- dheft(xx, fit)
   if(what == "h" || what == "H")
      yy <- hheft(xx, fit)
   if(what == "f" || what == "F" || what == "p" || what == "P")
      yy <- pheft(xx, fit)
   if(what == "s" || what == "S")
      yy <- 1-pheft(xx, fit)
   if(yy[1] < -8)
      stop("What is wrong? Well: what is wrong.")
   if(missing(xlab))
      xlab <- ""
   if(missing(ylab))
      ylab <- ""
   if(missing(type))
      type <- "l"
   if(!add)
      plot(xx, yy, xlim = xlim, xlab = xlab, ylab = ylab, type = type,
         ...)
   else lines(xx, yy, type = type, ...)
}
print.heft <- function(x,...)
{
      summary.heft(x)
}
summary.heft <- function(object,...)
{
    if(class(object)!="heft")
       stop("object is not a heft object")
        fit <- object 
   ul <- fit$penalty
   um <- fit$sample
   ll <- fit$logl
   kk <- (1:length(ll))
   kk <- kk[fit$ad != 2]
   ll <- ll[fit$ad != 2]
   ad <- fit$ad[fit$ad != 2]
   bb <- -2 * ll + ul * (kk-2)
   if(fit$thetal[1]!=0) bb <- bb+ul
   if(fit$thetal[2]!=0) bb <- bb+ul
   if(fit$thetap[2]!=0) bb <- bb+ul
   if(fit$thetap[2]==0 && min(kk)==2) bb <- bb+ul
   cc1 <- bb
   cc2 <- bb
   cc2[1] <- Inf
   cc1[length(bb)] <- 0
   if(length(bb) > 1) {
      for(i in 1:(length(bb) - 1)) {
         cc1[i] <- max((ll[(i + 1):(length(bb))] - ll[i])/(kk[(i +
            1):(length(bb))] - kk[i]))
         cc2[i + 1] <- min((ll[1:i] - ll[i + 1])/(kk[1:i] - kk[i +
            1]))
      }
   }
   c3 <- cc2 - cc1
   cc1[c3 < 0] <- NA
   cc2[c3 < 0] <- NA
   uu <- cbind(kk, ad, ll, bb, 2 * cc1, 2 * cc2)
   ww <- rep("", length(bb))
   dimnames(uu) <- list(ww, c("knots", "A(0)/D(1)", "loglik", "AIC", 
      "minimum penalty", "maximum penalty"))
   print(round(uu, 2))
   cat(paste("the present optimal number of knots is ", kk[bb == min(bb)], 
      "\n"))
   if(ul == log(um))
      cat(paste("penalty(AIC) was the default: BIC=log(samplesize): log(",
         um, ")=", round(ul, 2), "\n"))
   else cat(paste("penalty(AIC) was ", round(ul, 2), 
         ", the default (BIC) ", "would have been", round(log(um
         ), 2), "\n"))
   if(min(kk) == 3 && fit$thetap[2] != 0) {
      cat(paste("models with fewer than", kk[1], "knots", 
         "can be fitted, but they are not optimal for the\n"))
      cat(paste("present choice of penalty - choose penalty in", 
         "heft larger to see these fits\n"))
   }
   if(min(kk) > 3) {
      cat(paste("models with fewer than", kk[1], "knots", 
         "can be fitted, but they are not optimal for the\n"))
      cat(paste("present choice of penalty - choose penalty in", 
         "heft larger to see these fits\n"))
   }
   uuu <- matrix(NA, ncol = 3, nrow = 2, dimnames = list(c("left tail", 
      "right tail"), c("theta", "SE", "t")))
   uuu[, 1] <- fit$thetal
   if(fit$logse[1] > 0) {
      uuu[1, 2] <- fit$logse[1]
      uuu[1, 3] <- abs(fit$thetal[1]/fit$logse[1])
   }
   if(fit$logse[2] > 0) {
      uuu[2, 2] <- fit$logse[2]
      uuu[2, 3] <- abs(fit$thetal[2]/fit$logse[2])
   }
   print(round(uuu, 2))
   invisible()
}
dheft <- function(q, fit)
{
    if(class(fit)!="heft")
       stop("fit is not a heft object")
   y <- hheft(q, fit)
   z <- 1 - pheft(q, fit)
   y * z
}
hheft <- function(q, fit)
{
    if(class(fit)!="heft")
       stop("fit is not a heft object")
        q <- unstrip(q)
   y <- fit$thetap[1] + q * fit$thetap[2] + fit$thetal[1] * log(q/(q + fit$
      shift)) + fit$thetal[2] * log(q + fit$shift)
   for(i in 1:length(fit$knots)) {
      if(fit$thetak[i] != 0)
         y <- y + fit$thetak[i] * ((abs(q - fit$knots[i]) + q - 
            fit$knots[i])/2)^3
   }
   exp(y)
}
pheft <- function(q, fit)
{
    if(class(fit)!="heft")
       stop("fit is not a heft object")
        q <- unstrip(q)
   sq <- rank(q)
   q <- sort(q)
   z <- .C("heftpq",
      as.double(c(fit$knots,rep(0,100))),
      as.double(c(fit$shift,rep(0,100))),
      as.double(c(fit$thetak,rep(0,100))),
      as.double(c(fit$thetal,rep(0,100))),
      as.double(c(fit$thetap,rep(0,100))),
      as.integer(1),
      pp = as.double(q),
      as.double(q),
      as.integer(length(fit$knots)),
      as.integer(length(q)),
      PACKAGE = "polspline")
   zz <- z$pp[sq]
   zz[q < 0] <- 0
   zz
}
qheft <- function(p, fit)
{
    if(class(fit)!="heft")
       stop("fit is not a heft object")
        p <- unstrip(p)
   sp <- rank(p)
   p <- sort(p)
   z <- .C("heftpq",
      as.double(c(fit$knots,rep(0,100))),
      as.double(c(fit$shift,rep(0,100))),
      as.double(c(fit$thetak,rep(0,100))),
      as.double(c(fit$thetal,rep(0,100))),
      as.double(c(fit$thetap,rep(0,100))),
      as.integer(0),
      as.double(p),
      qq = as.double(p),
      as.integer(length(fit$knots)),
      as.integer(length(p)),
      PACKAGE = "polspline")
   zz <- z$qq[sp]
   zz[p < 0] <- NA
   zz[p == 0] <- 0
   zz[p == 1] <- Inf
   zz[p > 1] <- NA
   zz
}
rheft <- function(n, fit)
{
    if(class(fit)!="heft")
       stop("fit is not a heft object")
   pp <- runif(n)
   qheft(pp, fit)
}
oldlogspline.to.logspline <- function(obj,data)
{
   nobj <- list()
   nobj$call <- obj$call
   if(is.null(obj$call))nobj$call <- "translated from oldlogspline"
   nobj$knots <- sum(obj$coef[-(1:2)]!=0)
   nobj$coef.pol <- obj$coef[1:2]
   nobj$coef.kts <- obj$coef[-(1:2)]
   nobj$coef.kts <- nobj$coef.kts[nobj$coef.kts!=0]
   nobj$knots <- obj$knots[obj$coef[-(1:2)]!=0]
   nobj$maxknots <- length(obj$coef)-2
   nobj$penalty <- obj$penalty
   nobj$bound <- obj$bound
   nobj$samples <- obj$sample
   nobj$logl <- obj$logl[obj$logl!=0]
   lx <- length(nobj$logl)
   nobj$logl <- cbind(nobj$maxknots+1-(lx:1),c(rep(2,lx-1),1),nobj$logl)
   if(!missing(data))nobj$range <- obj$range
   else {
      lx <- 1/(nobj$samples+1)
      nobj$range <- qoldlogspline(c(lx,1-lx),obj)
   }
   nobj$mind
   class(nobj) <- "logspline"
   nobj
}
poldlogspline <- function(q, fit)
{
    if(class(fit)!="oldlogspline")
       stop("fit is not an oldlogspline object")
        q <- unstrip(q)
    sq <- rank(q)
    q <- sort(q)
    z <- .C("pqlsd",
        as.double(fit$coef),
        as.double(fit$knots),
        as.double(fit$bound),
        as.integer(1),
        pp = as.double(q),
        as.double(q),
        as.integer(length(fit$knots)),
        as.integer(length(q)),
      PACKAGE = "polspline")
    zz <- z$pp[sq]
    if(fit$bound[1] > 0)
        zz[q<fit$bound[2]] <- 0
    if(fit$bound[3] > 0)
        zz[q>fit$bound[4]] <- 1
    zz
}

qoldlogspline <- function(p, fit)
{
    if(class(fit)!="oldlogspline")
       stop("fit is not an oldlogspline object")
        p <- unstrip(p)
    sp <- rank(p)
    p <- sort(p)
    z <- .C("pqlsd",
        as.double(fit$coef),
        as.double(fit$knots),
        as.double(fit$bound),
        as.integer(0),
        as.double(p),
        qq = as.double(p),
        as.integer(length(fit$knots)),
        as.integer(length(p)),
      PACKAGE = "polspline")
    zz <- z$qq[sp]
    zz[p<0] <- NA
    zz[p>1] <- NA
    zz
}

roldlogspline <- function(n, fit)
{
    if(class(fit)!="oldlogspline")
       stop("fit is not an oldlogspline object")
    pp <- runif(n)
    qoldlogspline(pp, fit)
}

doldlogspline <- function(q, fit)
{
    x <- q
    if(class(fit)!="oldlogspline")
       stop("fit is not an oldlogspline object")
        q <- unstrip(q)
    y <- fit$coef[1] + x * fit$coef[2]
    for(i in 1:length(fit$knots)) {
        if(fit$coef[i+2] != 0)
            y <- y + fit$coef[i+2] * ((abs(x - fit$knots[i]) +
                                x - fit$knots[i])/2)^3
    }
    y <- exp(y)
    if(fit$bound[1] > 0)
            y[x < fit$bound[2]] <- 0
    if(fit$bound[3] > 0)
            y[x > fit$bound[4]] <- 0
    y
}

plot.oldlogspline <- function(x, n = 100, what = "d", xlim, xlab = "", ylab = "", type = "l", add = FALSE, ...)
{
    fit <- x
    if(class(fit)!="oldlogspline")
       stop("fit is not an oldlogspline object")
        if(missing(xlim)) {
                u1 <- qoldlogspline(0.01, fit)
                u2 <- qoldlogspline(0.99, fit)
                u3 <- 1.1 * u1 - 0.1 * u2
                u4 <- 1.1 * u2 - 0.1 * u1
        }
        else {
                u3 <- xlim[1]
                u4 <- xlim[2]
        }
        xx <- (0:(n - 1))/(n - 1) * (u4 - u3) + u3
        if(what == "d" || what == "D")
                yy <- doldlogspline(xx, fit)
        if(what == "f" || what == "F" || what == "p" || what == "P")
                yy <- poldlogspline(xx, fit)
        if(what == "s" || what == "S")
                yy <- 1 - poldlogspline(xx, fit)
        if(what == "h" || what == "H")
                yy <- doldlogspline(xx, fit)/(1 - poldlogspline(xx, fit))
        if(missing(xlab))
                xlab <- ""
        if(missing(ylab))
                ylab <- ""
        if(missing(type))
                type <- "l"
        if(add==FALSE)plot(xx, yy, xlab = xlab, ylab = ylab, type = type, ...)
        else lines(xx,yy, type = type, ...)
}
print.oldlogspline <- function(x,...)
{
   summary.oldlogspline(x)
}
summary.oldlogspline <- function(object,...)
{
    if(class(object)!="oldlogspline")
       stop("fit is not an oldlogspline object")
    fit <- object
    if(fit$delete==FALSE)stop(paste("summary.oldlogspline can only provide",
       "information if delete in oldlogspline is TRUE"))
    ul <- fit$penalty
    um <- fit$sample
    ll <- fit$logl
    kk <- (1:length(ll))
    kk <- kk[ll != 0] + 2
    ll <- ll[ll != 0]
    error<-FALSE
    rr <- ll[1:(length(ll)-1)]-ll[2:length(ll)]
    if(length(ll)>1 && max(rr)>0)error<-TRUE
    bb <- -2 * ll + ul * kk
    cc1 <- bb
    cc2 <- bb
    cc2[1] <- Inf
    cc1[length(bb)] <- 0
    if(length(bb) > 1) {
        for(i in 1:(length(bb) - 1)) {
            cc1[i] <- max((ll[(i + 1):(length(bb))] - ll[i])/(
                    kk[(i + 1):(length(bb))] - kk[i]))
            cc2[i + 1] <- min((ll[1:i] - ll[i + 1])/(kk[1:i] - kk[i + 1]))
        }
    }
    c3 <- cc2 - cc1
    cc1[c3 < 0] <- NA
    cc2[c3 < 0] <- NA
    uu <- cbind(kk, ll, bb, 2 * cc1, 2 * cc2)
    ww <- rep("", length(bb))
    if(error){
    cat("Warning - imprecision in loglikelihood (possibly due to heavy tails)\n")
    cat("the output of summary.oldlogspline might not be correct\n")
    }
    dimnames(uu) <- list(ww, c("knots", "loglik", "AIC", "minimum penalty",
        "maximum penalty"))
    print(round(uu, 2))
    cat(paste("the present optimal number of knots is ", kk[bb== min(bb)],"\n"))
    if(ul == log(um))
        cat(paste("penalty(AIC) was the default: BIC=log(samplesize): log(",
                um, ")=", round(ul, 2),"\n"))
    else
        cat(paste("penalty(AIC) was ", round(ul, 2),", the default (BIC) ",
                "would have been", round(log(um), 2),"\n"))
    if(min(kk) > 3 && fit$delete==TRUE){
        cat(paste( "models with fewer than", kk[1],"knots ", 
                  "can be fitted, but they are not optimal for\n"))
        cat(paste("the present choice of penalty - choose penalty in",
                  "oldlogspline larger\nto see these fits\n"))
    }
    if(min(kk) > 3 && fit$delete==3)
        cat(paste("models with fewer than", kk[1],"knots ",
                    "were not fitted because of convergence problems\n"))
      
    invisible()
}

oldlogspline <- function(uncensored, right, left, interval, lbound, ubound,
        nknots, knots, penalty, delete = TRUE)
{
    nsample <- rep(0, 6)
    # interval is the nterval censored data - a matrix with two columns
    if(!missing(uncensored))uncensored <- unstrip(uncensored)
    if(!missing(right))right <- unstrip(right)
    if(!missing(left))left <- unstrip(left)
    if(!missing(interval))interval <- unstrip(interval)
    if(!missing(knots))knots <- unstrip(knots)
    if(!missing(interval)) {
        if(length(interval[1,  ]) != 2)
            stop("interval must have two columns")
        if(min(abs(interval[, 1] - interval[, 2])) < 0) stop(
                   "not all lower bounds smaller than upper bounds")
        nsample[3] <- length(interval)/2
        nsample[1] <- length(interval)/2
        # grouping boundaries can not be beyond the boundaries of the density
        if(!missing(lbound))
            interval[interval[, 1] < lbound, 1] <- lbound
        if(!missing(ubound))
            interval[interval[, 2] > ubound, 2] <- ubound
        sample <- as.vector(t(interval))
        ror <- order(interval[,1],interval[,2])
        if(nsample[3]>1){
      ro1 <- interval[ror[(1:(nsample[3]-1))],1]==interval[ror[2:nsample[3]],1]
      ro2 <- interval[ror[(1:(nsample[3]-1))],2]==interval[ror[2:nsample[3]],2]
            nsample[6] <- nsample[3]-sum(ro1+ro2==2)
        }
        else nsample[6] <- 1
    }
# uncensored is the uncensored data
    if(!missing(uncensored)) {
        uncensored2 <- uncensored[!is.na(uncensored)]
        u2 <- length(uncensored) - length(uncensored2)
        if(u2 > 0)
            print(paste("***", u2, " NAs ignored in uncensored"))
        uncensored <- uncensored2
        if(nsample[1] > 0)
            sample <- c(uncensored, sample)
        if(nsample[1] == 0)
            sample <- uncensored
        nsample[1] <- length(uncensored) + nsample[1]
        nsample[2] <- length(uncensored)
        uncensored <- sort(uncensored)
        if(nsample[2]>1)
            nsample[6] <- sum(uncensored[2:nsample[2]] !=
                uncensored[1:(nsample[2]-1)]) + 1 + nsample[6]
        else 
            nsample[6] <- nsample[6]+1
    }
# we can not run on only right or left censored data
        if(nsample[1] == 0) stop("you either need uncensored or interval censored data")
        # right is the right censored data
        if(!missing(right)) {
                if(nsample[1] > 0)
                        sample <- c(sample, right)
                if(nsample[1] == 0)
                        sample <- right
                nsample[1] <- length(right) + nsample[1]
                nsample[4] <- length(right)
                right <- sort(right)
                if(nsample[4]>1){
                nsample[6] <- sum(right[2:nsample[4]]!=right[1:(nsample[4]-1)])+
                          1 + nsample[6]
                }
                else nsample[6] <- nsample[6]+1
        }
# left is the left censored data
        if(!missing(left)) {
                if(nsample[1] > 0)
                        sample <- c(sample, left)
                if(nsample[1] == 0)
                        sample <- left
                nsample[1] <- length(left) + nsample[1]
                nsample[5] <- length(left)
                left <- sort(left)
                if(nsample[5]>1){
                nsample[6] <- sum(left[2:nsample[5]]!=left[1:(nsample[5]-1)])+
                          1 + nsample[6]
                }
                else nsample[6] <- nsample[6]+1
        }
# the default for penalty is bic: log(length(sample))
        if(missing(penalty)) penalty <- log(nsample[1])
        n1 <- 4 * nsample[1]^0.2 + 1
        if(!missing(nknots))
                n1 <- nknots + 1
        if(!missing(knots)) n1 <- length(knots) + 1      # user provides knots
        if(!missing(knots)) {
                nknots <- length(knots)
                knots <- sort(knots)
                iautoknot <- 0
                if(knots[1] > min(sample))
                        stop("first knot must be smaller than smallest sample")

                if(knots[nknots] < max(sample))
                        stop("last knot should be larger than largest sample")

        }
        else {
                if(missing(nknots))
                        nknots <- 0
                knots <- vector(mode = "double", length = max(nknots, 50))
                iautoknot <- 1
        }
        xbound <- c(1, 0, 0, 0, 0)
        if(!missing(lbound)) {
                xbound[2] <- 1
                xbound[3] <- lbound
                if(lbound > min(sample))
                        stop("lbound should be smaller than smallest sample")
        }
        if(!missing(ubound)) {
                xbound[4] <- 1
                xbound[5] <- ubound
                if(ubound < max(sample))
                        stop("ubound should be larger than largest sample")
        }
# SorC will carry the error messages - in code form
        SorC <- vector(mode = "integer", length = 35)
        SorC[1] <- 1    # the actual function call
        nsample[6] <- nsample[6]-1
        z <- .C("logcensor",
                as.integer(delete),
                as.integer(iautoknot),
                as.double(sample),
                as.integer(nsample),
                bd = as.double(xbound),
                SorC = as.integer(SorC),
                nk = as.integer(nknots),
                kt = as.double(knots),
                cf = as.double(c(knots, 0, 0)),
                as.double(penalty),
                as.double(sample),
                as.double(sample),
                logl = as.double(rep(0, n1 + 1)),
      PACKAGE = "polspline")
        bound <- c(z$bd[2], z$bd[3], z$bd[4], z$bd[5])
        SorC <- z$SorC  # error messages
        if(abs(SorC[1]) > 2) {
                for(i in 3:abs(SorC[1]))
                        cat(paste("===> warning: knot ", SorC[i - 1],
                                " removed - double knot\n"))
                if(SorC[1] < 0)
                        SorC[1] <- -1
                if(SorC[1] == 23)
                        SorC[1] <- -3
        }
        if(abs(SorC[1]) > 3) {
                cat("* several double knots suggests that your data is *\n")
                cat("* strongly rounded, attention might be required   *\n")
                SorC[1] <- 1
        }
        if(SorC[1] == -3)
                stop("* too many double knots")
        if(SorC[1] == -1 && SorC[28] == 0)
                stop("* no convergence")
        if(SorC[28] > 0)
                cat(paste("* convergence problems, smallest number of knots",
                        " tried is ", SorC[28] + 1," *\n"))
        if(SorC[1] == 2)
                stop("* sample is too small")
        if(SorC[1] == -2)
                stop(paste("* too many knots, at most ", SorC[2],
                        "knots possible"))
        if(SorC[22] == 1) {
                cat("possible discontinuity at lower end\n")
                cat(paste("consider rerunning with lbound=", z$kt[1],
         "\n"))

        }
        if(SorC[22] == 3) {
                cat("possible infinite density at lower end\n")
                cat("running program with fewer knots\n")
        }
        if(SorC[21] == 1)
                cat("running with maximum degrees of freedom\n")
        if(SorC[25] >0)
               cat("* problems are possibly due to a very heavy right tail *\n")
        if(SorC[24] >0)
                cat("* problems are possibly due to a very heavy left tail *\n")
        if(SorC[23] == 3) {
                cat("possible infinite density at upper end\n")
                cat("running program with fewer knots\n")
        }
        if(SorC[23] == 1) {
                cat("possible discontinuity at upper end\n")
                cat(paste("consider rerunning with ubound=", z$kt[z$nk],
         "\n"))

        }
        if(delete && SorC[28]>0)delete<-3
        coef <- z$cf[1:(z$nk + 2)]
        uu <- 3:z$nk
        if(delete == FALSE)uu <- 1
        fit <- list(coef = coef, knots = z$kt[1:z$nk], bound = bound, logl = z$logl[
                uu], penalty = penalty, sample = nsample[1], delete = delete)
        class(fit) <- "oldlogspline"
        fit
}

lspec <- function(data, period, penalty, minmass, knots, maxknots, atoms, 
        maxatoms, maxdim, odd = FALSE, updown = 3,silent=TRUE)
{
   call <- match.call()
        if(!missing(data))data <- unstrip(data)
        if(!missing(period))period <- unstrip(period)
        if(!missing(knots))knots <- unstrip(knots)
        if(!missing(atoms))atoms <- unstrip(atoms)
   if(missing(period) && missing(data))
      stop(" either data or period should be specified ")
   if(!missing(period) && !missing(data))
      stop(" only one of data or period should be specified ")
   if(!missing(period))
      ny <- 2 * length(period)
   if(missing(period)) {
      ny <- length(data)
      period <- Mod(fft(data))^2/(ny * 2 * pi)
      period <- period[1:floor((length(period) + 2)/2)]
      odd <- TRUE
      if(floor(ny/2) == ny/2)
         odd <- FALSE
   }
        else{
         if(odd) ny <- ny + 1
           period <- c(1,period)
        }
   if(min(period) <= 0)
      stop(" all period elements should be larger than 0 ")
   if(length(period) < 10)
      stop("too few observations")
   z <- .C("tspspsx",
      z = as.integer(rep(-1, 12)),
      PACKAGE = "polspline")
   lgth <- z$z[1]
   nx <- length(period)
   if(missing(penalty))
      penalty <- log(nx - 1)
   dimatt <- 0
   ktsatt <- 1
   spkatt <- 1
   nknots <- 0
   natoms <- 0
   if(!missing(maxknots)){
      maxknots <- max(1, maxknots)
        }
   else {
      maxknots <- -1
      ktsatt <- 0
   }
   if(!missing(maxatoms)){
      maxatoms <- max(0, maxatoms)
        }
   else {
      maxatoms <- -1
      spkatt <- 0
   }
   if(missing(minmass)){
      if(!missing(data))
         minmass <- var(data)*(-log(1-0.95^(1/nx))-1)/ny
      else{
         minmass <- mean(period[2:length(period)])*2*pi
         minmass <- minmass*(-log(1-0.95^(1/nx))-1)/ny
                }
        }
   minmass <- minmass * ny /(2*pi)
   if(!missing(knots)) {
      nknots <- length(knots)
      if(nknots>1){
      uu <- knots[2:nknots] - knots[1:(nknots - 1)]
      if(min(uu) <= 0)
         stop("knots not in sequence")
      }
      if(knots[1] < 0)
         stop("knot 1 too small")
      if(knots[nknots] > pi)
         stop("last knot too large")
      knots <- c(knots, rep(0, lgth - nknots))
      if(ktsatt * maxknots < ktsatt * nknots)
         stop("more knots than maxknots")
   }
        else{
      knots <- rep(0, lgth)
        }
   if(!missing(atoms)) {
      natoms <- length(atoms)
      atoms <- round((atoms * ny)/(2 * pi))
      if(natoms>1){
      uu <- atoms[2:natoms] - knots[1:(natoms - 1)]
      if(min(uu) <= 0)
         stop("atoms not in sequence or too close")
      }
      if(atoms[1] < 1)
         stop("atom 1 too small")
      if(atoms[natoms] > ny/2)
         stop("last atom too large")
      atoms <- c(atoms, rep(0, lgth - natoms))
      if(spkatt * maxatoms < spkatt * natoms)
         stop("more atoms than maxatoms")
   }
        else{
      atoms <- rep(0, lgth)
       }
   u1 <- max(nknots, 1, maxknots) + max(natoms, maxatoms)
   if(u1 > lgth - 5)
      stop("too many dimensions")
   if(!missing(maxdim)) {
      dimatt <- 1
      if(u1 > maxdim)
         stop("maxdim too small for other specifications")
      if(maxdim > lgth - 5)
         stop(paste("maxdim can be at most", lgth - 5))
   }
        else{
      maxdim <- max(4 * nx^0.2, 15, u1)
        }
   dims <- c(nx, maxdim, dimatt, maxknots, ktsatt, nknots, maxatoms, 
      spkatt, natoms, odd, updown, 1*silent, 0)
   z <- .C("tspsps",
      dims = as.integer(dims),
      as.double(period),
      knots = as.double(c(knots,rep(0,nx))),
      atoms = as.integer(c(atoms,rep(0,nx))),
      as.double(penalty),
      logl = as.double(rep(0, lgth)),
      theta = as.double(rep(0, lgth)),
      ad = as.integer(rep(0, lgth)),
      minmass = as.double(minmass),
      PACKAGE = "polspline")
   dims <- z$dims
   minmass <- minmass /( ny /(2*pi))
   if(dims[12] == 1)
      stop(paste("numerical problems -\n", 
         "probably too many knots or knots too close together", 
         " or a very sharp atom"))
   if(dims[12] == 2)
      stop("no convergence")
   z$logl[abs(z$logl) < 1e-100] <- 0
   z$logl[z$ad == 2] <- 0
   mass <- ((z$theta[(dims[6] + 4) + (1:dims[9])]) * 2 * pi)/ny
   atoms <- (z$atoms[1:dims[9]] * 2 * pi)/ny
   if(dims[9] == 0) {
      mass <- 0
      atoms <- 0
   }
   thetak <- z$theta[5:(dims[6] + 4)]
   knots <- z$knots[1:dims[6]]
   if(dims[6] == 0) {
      thetak <- 0
      knots <- 0
   }
   logl <- z$logl[dims[6]+dims[9]]
   fit <- list(call = call, thetap = z$theta[1:4], nknots = dims[6], knots = 
      knots, thetak = thetak, natoms = dims[9], atoms = atoms, 
      mass = mass, penalty = penalty, minmass = minmass,
      sample = ny, logl = logl, updown = dims[11])
        class(fit) <- "lspec"
        fit
}
clspec <- function(lag, fit, cov = TRUE, mm)
{
        if(class(fit)!="lspec")
          stop("fit is not an lspec object")
        if(!missing(lag))lag <- unstrip(lag)
   llag <- abs(lag)
   if(max(abs(round(llag) - llag)) > 0.01)
      stop("some lags are not integer")
   if(missing(mm)) {
      mm <- max(c(1024, fit$sample, max(llag + 1)))
      mm <- 2^(1 + floor(log(mm - 0.1)/log(2)))
   }
   if(mm < max(llag + 1))
      stop("mm too small")
   rr <- dlspec(((0:mm) * pi)/mm, fit)$d
   rr <- c(rr, rr[mm:2])
   rr <- (Re(fft(rr)) * pi)/mm
   rr <- rr[llag + 1]
        if(fit$natoms>0){
   for(i in 1:fit$natoms) {
      rr <- rr + 2 * cos(lag * fit$atoms[i]) * fit$mass[i]
   }
   }
   if(cov == FALSE)
      rr <- rr/rr[1]
   rr
}
dlspec <- function(freq, fit)
{
        if(class(fit)!="lspec")
          stop("fit is not an lspec object")
        if(!missing(freq))freq <- unstrip(freq)
   freq <- freq - floor(freq/(2 * pi)) * 2 * pi
   freq[freq > pi] <- 2 * pi - freq[freq > pi]
   y <- rep(fit$thetap[1], length(freq)) + freq * fit$thetap[2]
   y <- y + freq^2 * fit$thetap[3] + freq^3 * fit$thetap[4]
   if(fit$nknots > 0) {
      for(i in 1:fit$nknots) {
         z <- freq - fit$knots[i]
         y[z > 0] <- y[z > 0] + z[z > 0]^3 * fit$thetak[i]
      }
   }
   d1 <- exp(y)
   modfreq <- round((freq * fit$sample)/(2 * pi))
   modmatch <- round((fit$atoms * fit$sample)/(2 * pi))
   uu <- rep(0, round(fit$sample/2) + 2)
   uu[modmatch] <- fit$mass
   uu <- c(NA, uu)
   l1 <- uu[modfreq+1]
   modfreq <- ((2 * pi)/fit$sample) * modfreq
   list(d = d1, modfreq = modfreq, m = l1)
}
plspec <- function(freq, fit, mm)
{
        if(class(fit)!="lspec")
          stop("fit is not an lspec object")
        if(!missing(freq))freq <- unstrip(freq)
        if(missing(mm)) {
                mm <- max(c(4096, fit$sample))
                mm <- 2^floor(log(mm - 0.1)/log(2))
        }
        ff <- freq[freq >=  - pi]
        ff <- ff[ff <= pi]
        gg <- c(abs(ff), pi)
        uu <- (c((1:mm) - 0.5) * pi)/mm
        tt <- dlspec(uu, fit)$d
        ss <- cumsum(tt)/mm
        ss <- (c(0, 0, ss, ss[mm]) * pi)
        tt <- (gg * mm)/pi
        vv <- floor(tt)
        tt <- tt - vv
        tt <- (1 - tt) * ss[vv + 2] + tt * ss[vv + 3]
        if(fit$natoms > 0) {
                for(i in 1:fit$natoms)
                        tt[gg >= fit$atoms[i]] <- tt[gg >= fit$atoms[i]] + fit$
                                mass[i]
        }
        if(length(gg) < length(freq) + 1 || is.na(gg[1]))
           warning("plspec is only valid for frequencies between -pi and pi")
        ss <- rep(NA, length(freq))
        ss[abs(freq) <= pi] <- tt[ - length(tt)] + tt[length(tt)]
        ss[freq < 0] <- 2 * tt[length(tt)] - ss[freq < 0]
        if(fit$natoms > 0) {
                for(i in 1:fit$natoms)
                        ss[freq ==  - fit$atoms[i]] <- ss[freq ==  - fit$atoms[
                                i]] + fit$mass[i]
        }
        ss
}
rlspec <- function(n, fit, mean = 0, cosmodel = FALSE, mm)
{
        if(class(fit)!="lspec")
          stop("fit is not an lspec object")
        if(missing(mm)) {
                mm <- max(c(1024, fit$sample, n))
                mm <- 2^(1 + floor(log(mm - 0.1)/log(2)))
        }
        if(mm < max(n/2 + 1))
                stop("mm too small")
        rr <- (dlspec(((0:mm) * pi)/mm, fit)$d * pi)/(2 * mm)
        rr[1] <- rr[1] * 2
        rr[mm + 1] <- rr[mm + 1] * 2
        rr <- sqrt(rr)
        uu <- rnorm(rr, 0, rr)
        uu <- c(uu, uu[mm:2])
        vv <- rnorm(rr, 0, rr)
        vv <- c(vv,  - vv[mm:2])
        vv[c(1, (mm + 1))] <- 0
        uu <- uu + vv * (1i)
        uu <- Re(fft(uu))
        uu <- uu[1:n] + mean
        if(fit$natoms > 0) {
      cc <- runif(1)*2*pi-pi
      if(cosmodel) aa <- 2*sqrt(fit$mass)
      else
                aa <- 2 * rnorm(fit$natoms, 0, sqrt(fit$mass))
                aa[fit$atoms == pi] <- 2 * aa[fit$atoms == pi]
                for(i in 1:fit$natoms)
                        uu <- uu + aa[i] * cos((1:n) * fit$atoms[i] + pi * cc)
        }
        uu
}
plot.lspec <- function(x, what = "b", n, add = FALSE, xlim, ylim, xlab, ylab, type, ...)
{
        fit <- x
        if(class(fit)!="lspec")
          stop("fit is not an lspec object")
   if(add) {
      plim <- (par()$usr)[1:2]
      if(!missing(xlim)) {
         plim[1] <- max(xlim[1], plim[1])
         plim[2] <- min(xlim[2], plim[2])
      }
   }
   else {
      plim <- c(0, pi)
      if(what =="p"||what=="P"||what=="f"||what=="F")plim[1]<- -pi
      if(!missing(xlim)) {
         plim[1] <- xlim[1]
         plim[2] <- xlim[2]
      }
   }
   if(missing(xlab))
      xlab <- ""
   if(missing(ylab))
      ylab <- ""
   if(what == "l" || what == "L") {
      if(missing(type))
         type <- "h"
      if(fit$natoms>0){
      x5 <- c(-fit$atoms,fit$atoms)
      tt <- round(plim[2]/(2*pi))+1
      vv <- round(plim[1]/(2*pi))-1
      x1 <- x5
      for(i in vv:tt)
         if(i!=0)x1 <- c(x1, x5+i*2*pi)
      y1 <- dlspec(x1,fit)$m
      y1 <- y1[x1 <= plim[2]]
      x1 <- x1[x1 <= plim[2]]
      y1 <- y1[x1 > plim[1]]
      x1 <- x1[x1 > plim[1]]
      x1 <- c(x1[1], x1)
      y1 <- c(0, y1)
      if(!add)
         plot(x1, y1, xlim = plim, xlab = xlab, ylab = ylab, 
            type = type, ...)
      else lines(x1, y1, type = type, ...)
      abline(h = 0)}
      else{
      if(add) abline(h=0)
      else plot(plim,c(0,0),xlab = xlab, ylab = ylab,type="l",...)
                }
   }
   if(what == "d" || what == "D" || what == "b" || what == "B") {
      if(missing(type))
         type <- "l"
      if(missing(n))
         n <- max(100, fit$sample + 1)
      xx <- (0:(n - 1))/(n - 1) * (plim[2] - plim[1]) + plim[1]
      yy <- dlspec(xx, fit)$d
      if(fit$natoms == 0)
         what <- "d"
      if(missing(ylim))ylim<-range(yy)
   }
   if(what == "b" || what == "B") {
      type <- "l"
                
      x5 <- c(-fit$atoms,fit$atoms)
      tt <- round(plim[2]/(2*pi))+1
      vv <- round(plim[1]/(2*pi))-1
      x3 <- x5
      for(i in vv:tt)
         if(i!=0)x3 <- c(x3, x5+i*2*pi)
      y3 <- dlspec(x3, fit)
      y3 <- max(yy)*1.1
      if(fit$nknots==1)y3 <- 2*y3
      if(!missing(ylim))y3 <- ylim[2]
      x2 <- x3 
      y2 <- dlspec(x2, fit)$d
      x4 <- x3 
      y4 <- y2
      for(i in 1:length(x3)) {
         yy <- c(yy[xx < x2[i]], y2[i], y3, y4[i], yy[xx > x4[
            i]])
         xx <- c(xx[xx < x2[i]], x2[i], x3[i], x4[i], xx[xx > x4[
            i]])
      }
      if(missing(ylim))ylim<-range(yy)
      yy <- yy[xx >= plim[1]]
      xx <- xx[xx >= plim[1]]
      yy <- yy[xx <= plim[2]]
      xx <- xx[xx <= plim[2]]
      y2 <- y2[x2 >= plim[1]]
      x2 <- x2[x2 >= plim[1]]
      y2 <- y2[x2 <= plim[2]]
      x2 <- x2[x2 <= plim[2]]
   }
   if(what == "f" || what == "F" || what == "p" || what == "P") {
      if(!missing(xlim)){
         if(xlim[1]< -pi || xlim[2]>pi)
      stop("for this plot the range cannot strecth beyond (-pi,pi)")
      }
      if(missing(xlim)){
         plim[1] <- max(plim[1],-pi)
         plim[2] <- min(plim[2],pi)
      }
      if(missing(type))
         type <- "l"
      if(missing(n))
         n <- max(100, fit$sample + 1)
      xx <- (0:(n - 1))/(n - 1) * (plim[2] - plim[1]) + plim[1]
      yy <- plspec(xx, fit)
      if(missing(ylim))ylim<-range(yy)
      if(fit$natoms > 0) {
         x2 <- fit$atoms
         y3 <- plspec(x2, fit)
         y2 <- y3 - fit$mass
         for(i in 1:fit$natoms) {
            yy <- c(yy[xx < x2[i]], y2[i], y3[i], yy[xx > 
              x2[i]])
            xx <- c(xx[xx < x2[i]], x2[i], x2[i], xx[xx > 
              x2[i]])
         }
         x2 <- -fit$atoms
         y3 <- plspec(x2, fit)
         y2 <- y3 - fit$mass
         for(i in 1:fit$natoms) {
            yy <- c(yy[xx < x2[i]], y2[i], y3[i], yy[xx > 
              x2[i]])
            xx <- c(xx[xx < x2[i]], x2[i], x2[i], xx[xx > 
              x2[i]])
         }
         yy <- yy[xx >= plim[1]]
         xx <- xx[xx >= plim[1]]
         yy <- yy[xx <= plim[2]]
         xx <- xx[xx <= plim[2]]
      }
   }
   if(what != "l" && what != "L") {
      if(!add)
         plot(xx, yy, xlim = plim, xlab = xlab, ylab = ylab, 
            type = type, ylim = ylim, ...)
      else lines(xx, yy, type = type, ...)
      if(what =="b" || what=="B")points(x2,y2)
   }
   invisible()
}
print.lspec <- function(x,...)
{
        summary.lspec(x)
}
summary.lspec <- function(object,...)
{
        fit <- object
        if(class(fit)!="lspec")
          stop("fit is not an lspec object")
   aa <- " Logspline Spectral Estimation\n"
   aa <- paste(aa,"=============================\n")
   aa <- paste(aa,"The fit was obtained by the command:\n ")
   cat(aa)
   print(fit$call)
   aic <- round(-2*fit$logl+fit$penalty*(fit$nknots+fit$natoms),2)
   logl <- round(fit$logl,2)
   ns <- fit$natoms
   nk <- fit$nknots
   nd <- ns + nk
   if(ns==0 && nk==1)
   aa <- paste(" Only 1 basis function, a constant, was fitted.\n")
   if(ns==0 && nk>1)
   aa <- paste(" A spline with",nk,"knots, was fitted;",
      "there were no lines in the model.\n")
   if(ns>0 && nk>1)
   aa <- paste(" A spline with",nk,"knots, was fitted;",
      "there were also",ns,"lines in the model.\n")
   if(ns>0 && nk==1)
   aa <- paste(" There were",nd,"basisfunctions, a constant and",
      ns,"lines, in the model.\n")
   aa <- paste(aa,"The log-likelihood of the model was",logl,
      "which corresponds to an AIC\n value of",aic,".\n\n")
   aa <- paste(aa,"The program went though",abs(fit$updown))
   if(fit$updown>0)
   aa <-paste(aa,"updown cycles, and reached a stable solution.\n")
   if(fit$updown<0)
   aa <-paste(aa,"updown cycles, and did not reach a stable solution.\n")
   p1 <- round(fit$penalty,2)
   n1 <- round(fit$minmass,4)
   nn <- floor(fit$sample/2)
   p2 <- round(log(nn),2)
   uu <- plspec(pi,fit)
   n2 <- round(uu*(-log(1-0.95^(1/nn))-1)/fit$sample,4)
   p3 <- (p1==p2)
   p4 <- TRUE
   if(n1/n2 > 1.2 || n2/n1 > 1.2) p4 <- FALSE
   if(p3==TRUE && p4==TRUE)aa<-paste(aa,
"Both penalty (AIC) and minmass were the default values. For penalty this\n",
"was log(n)=log(",nn,")=",p1," (as in BIC) and for minmass this was",n1,".\n")
   if(p3==TRUE && p4==FALSE)aa<-paste(aa,
   "Penalty (AIC) had the default values",
   "log(n)=log(",nn,")=",p1," (as in BIC).\n Minmass was",n1,
   ", the default would have been",n2,".\n")
   if(p3==FALSE && p4==FALSE)aa<-paste(aa,
   "Penalty was",p1,", the default would have been",
   "log(n)=log(",nn,")=",p2,"\n(as in BIC). Minmass was",n1,
   ", the default would have been",n2,".\n")
   if(p3==FALSE && p4==TRUE)aa<-paste(aa,
   "Penalty was",p1,", the default would have been, log(n)=log(",nn,")=",
        p2,"\n (as in BIC). Minmass was the default",n1,".\n\n")
   if(nk>1){aa<-paste(aa,"The locations of the knots were:")
             for(i in 1:nk)aa<-paste(aa,round(fit$knots[i],3))
      aa<-paste(aa,"\n")
   }
        if(ns>0){
           aa<-paste(aa,"The locations and the mass in each line were:\n")
   bb <- matrix(0,ncol=4,nrow=ns)
   for(i in 1:ns){
   bb[i,1]<-round(fit$atoms[i],3)
        bb[i,2]<-2*pi/(fit$atoms[i])
        bb[i,2]<-round(bb[i,2],2)
   bb[i,3]<-round(fit$mass[i],5)
   bb[i,4]<- round(100*fit$mass[i]/uu,2)
   }
   dimnames(bb) <- list(rep("",ns),c("angular frequency","period","mass",
        "% of total mass"))
   }
   cat(aa)
   if(ns>0)print(bb)
   invisible()
   
}

polymars <- function(responses, predictors, maxsize, gcv = 4., additive = FALSE, startmodel,
   weights, no.interact, knots, knot.space = 3, ts.resp, ts.pred, 
   ts.weights, classify, factors, tolerance = 1e-06, verbose = FALSE)
{
   #responses  - a vector (or matrix) of responses. (Can be a a vector of characters for classification)
   #predictors - a matrix of predictors with same number of cases as response. Columns are predictors.
   #OPTIONAL ARGUEMENTS
   #maxsize    - maximum number of basis function the model can contain 
   #gcv        - parameter for overall best model seletion
   #additive   - boolean, is the model to be additive
   #startmodel - either a matrix (m*4 or m*5) or a polymars object from a previous call to polymars 
   #             an initial model the procedure should start with in model selection
   #weights    - a vector of length equal to the number of cases
   #no.interact- a 2*l matrix of columns numbers of the predictor matrix( each row pair cannot 
   #              have interaction terms)
   #knots      - a vector specifying many knots per predictor are wanted (with -1 for categorical 
   #             variables) ncol(predictors)==length(knots), or a matrix with ncol(predictors) == 
   #             ncol(knots) with actual knot specified and filled out with NA's.
   #             Can also be a single number - "knots" number of knots per predictor                    
   #knot.space - minimum number of order statistics between knots
   #ts.resp    - testset reponses, same format as responses
   #ts.pred    - testset predictors, same format as predictors
   #ts.weights - testset weights, same format as weights
   #classify   - whether classification is to be done, set = TRUE if the response vector is integer, if 
   #             if character classify is automatically true
   #factors    - a vector of column numbers of the predictor matrix of categorical variables
   #tolerance  - a numerical parameter which may need to be made smaller if the program crashes
   #store the call to the polymars function
   call <- match.call()
   ism0 <- missing(classify)
   ism1 <- missing(ts.resp)
   ism2 <- missing(maxsize)
   ism3 <- missing(ts.pred)
   ism4 <- missing(ts.weights)
   ism5 <- missing(knots)
   ism6 <- missing(factors)
   ism7 <- missing(startmodel)
   ism8 <- missing(weights)
   ism9 <- missing(no.interact)
   if(!missing(responses))
      responses <- unstrip(responses)
   if(!missing(predictors))
      predictors <- unstrip(predictors)
   if(!missing(weights))
      weights <- unstrip(weights)
   if(!missing(no.interact))
      no.interact <- unstrip(no.interact)
   if(!missing(knots))
      knots <- unstrip(knots)
   if(!missing(ts.resp))
      ts.resp <- unstrip(ts.resp)
   if(!missing(ts.pred))
      ts.pred <- unstrip(ts.pred)
   if(!missing(ts.weights))
      ts.weights <- unstrip(ts.weights)
   if(!missing(factors))
      factors <- unstrip(factors)
   responses <- as.matrix(responses)
   predictors <- data.matrix(predictors)
   nresponses <- ncol(responses)
   npredictors <- ncol(predictors)
   ncases <- nrow(predictors)
   if(ism0)
      classify <- FALSE
   if(mode(responses) == "character" || classify == TRUE) {
      if(ncol(responses) > 1) {
         stop("When using character responses  or classify = TRUE only 1 response per case is allowed\n"
            )
      }
      char.responses <- responses
      int.responses <- as.integer(as.factor(responses))
      nresponses <- length(unique(responses))
      responses <- matrix(ncol = nresponses, nrow = ncases, data = 
         int.responses)
      for(i in 1:nresponses) {
         responses[, i] <- (responses[, i] == (unique(
            int.responses)[i]))
      }
      conversion <- matrix(ncol = 2, nrow = nresponses, c(unique(
         char.responses), unique(int.responses)))
      classify <- TRUE
      if(!ism1) {
         char.responses.test <- ts.resp
         ts.resp <- matrix(ncol = nresponses, nrow = length(
            char.responses.test), data = 0)
         for(i in 1:nresponses) {
            ts.resp[, i] <- as.integer(char.responses.test ==
               conversion[i, 1])
         }
      }
   }
   else {
      conversion <- FALSE
      classify <- FALSE
   }
   #maxsize that the model can grow to
   if(ism2) maxsize <- ceiling(min(6 * (ncases^(1/3)), ncases/4, 100))
   #if a testset is to be used in model selection
   if(!ism1 || !ism3) {
      if(ism1 || ism3) {
         stop("Both ts.resp (testsets responses) and ts.pred (testset predictors) should be specified\n"
            )
      }
      if(!is.matrix(ts.resp))
         ts.resp <- as.matrix(ts.resp)
      if(!is.matrix(ts.pred))
         ts.pred <- as.matrix(ts.pred)
      if(ncol(ts.resp) != nresponses) {
         stop("Testset should have the same number of responses as the training set\n "
            )
      }
      if(ncol(ts.pred) != npredictors) {
         stop("Testset should have the same number of predictors as the training set\n "
            )
      }
      if(nrow(ts.resp) != nrow(ts.pred)) {
         stop("Testset ts.pred and ts.resp should have the same number of cases (rows)"
            )
      }
      testsetmatrix <- cbind(ts.resp, ts.pred)
      testsetcases <- nrow(testsetmatrix)
      testset <- TRUE
      if(!ism4) {
         if(length(ts.weights) != testsetcases) {
            stop("length of testset weights misspecified\n"
               )
         }
         testset.weighted <- TRUE
         testsetmatrix <- cbind(ts.resp * ts.weights, ts.pred)
      }
      else {
         testset.weighted <- FALSE
         ts.weights <- 0
      }
   }
   else {
      testsetmatrix <- 0
      testsetcases <- 0
      testset <- FALSE
      testset.weighted <- FALSE
      ts.weights <- 0
   }
   #If the mesh is specified by the knots arguement this will be changed to
   #true later
   mesh.specified <- FALSE
   mesh.vector <- 0
   if(nrow(responses) != nrow(predictors)) {
      stop("The number of rows (cases) of the response and predictor matricies should be the same"
         )
   }
   if(!ism5 && !is.matrix(knots) && length(knots) != npredictors && length(
      knots) != 1) {
      stop("Length of vector of `knots per predictor' should be equal to number of predictors or 1\n"
         )
   }
   if(!ism5) {
      if(!is.matrix(knots)) {
         #if knots is specified as a single number it is expanded to a vector 
         #length npredictors
         if(length(knots) == 1) {
            knots <- rep(knots, npredictors)
            if(!ism6) {
               for(i in 1:length(factors)) {
                  if(!is.vector(factors)) {
                     stop("`factors' should be a vector whose elements are indicies of predictors that are factors\n"
                        )
                  }
                  # in knots the number of knots(per predictor) is specified
                  # or -1 if the predictor is a factor and all it values are levels  
                  knots[factors[i]] <- -1
               }
            }
         }
      }
      else {
         mesh <- knots
         mesh.vector <- vector(length = ncol(mesh) * nrow(mesh),
            mode = "double")
         knots <- vector(length = npredictors, mode = "integer")
         k <- 0
         for(i in 1:npredictors) {
            knots[i] <- length(unique(mesh[is.na(mesh[
               , i]) == FALSE, i]))
            for(j in 1:knots[i]) {
               k <- k + 1
               mesh.vector[k] <- unique(mesh[!is.na(
                  mesh[, i]), i])[j]
            }
         }
         if(!ism6) {
            for(i in 1:length(factors)) {
               if(!is.vector(factors)) {
                  stop("`factors' should be a vector whose elements are indicies of predictors that are factors\n"
                     )
               }
               # in knots the number of knots(per predictor) is specified
               # or -1 if the predictor is a factor and all it values are levels  
               knots[factors[i]] <- -1
            }
         }
         mesh.specified <- TRUE
      }
   }
   if(ism5) {
      knots <- rep(min(20, round(ncases/4)), npredictors)
      if(!ism6) {
         for(i in 1:length(factors)) {
            if(!is.vector(factors)) {
               stop("`factors' should be a vector whose elements are indicies of predictors that are factors\n"
                  )
            }
            # in knots the number of knots(per predictor) is specified
            # or -1 if the predictor is a factor and all it values are levels  
            knots[factors[i]] <- -1
         }
      }
   }
   startmodelsize <- 1
   #A starting model must be specified as a object of class polymars
   #or a matrix with 4 or 5 columns
   no.remove <- 0
   no.remove.size <- 0
   if(!ism7) {
      if(is.vector(startmodel))
         startmodel <- t(as.matrix(startmodel))
      v1 <- (class(startmodel) == "polymars")
      if(length(v1) == 0)
         v1 <- FALSE
      if(!(is.matrix(startmodel) || v1) || (is.matrix(startmodel) &&
         (ncol(startmodel) != 4 && (ncol(startmodel) != 5)))) {
         stop(paste(
            "startmodel should be a matrix with each row corresponding to",
            "a function with number of columns = 4 (or 5 for extra boolean\n",
            "column specifying predictors which cannot be removed)",
            "or startmodel should be a polymars object\n"))
      }
      if(is.matrix(startmodel)) {
         #Fifth column denotes which basis functions must remain in the model at 
         #all times
         if(ncol(startmodel) == 5) {
            no.remove <- vector(length = (nrow(startmodel))
               )
            j <- 0
            for(i in 1:nrow(startmodel)) {
               if(startmodel[i, 5] == TRUE) {
                  j <- j + 1
                  no.remove[j] <- i
               }
            }
            no.remove.size <- j
         }
         #The startknots are taken from the startmodel and put into a vector
         #The startmodel becomes a 4*n matrix with a "1" in the 2nd and 4th 
         #columns where knots appear
         startknots <- as.vector(t(cbind(startmodel[, 2], 
            startmodel[, 4])))
         startknots[is.na(startknots)] <- 0.
         startmodel <- matrix(startmodel[, 1:4], ncol = 4)
         startmodel[!is.na(startmodel[, 2]), 2] <- 1
         startmodel[is.na(startmodel[, 2]), 2] <- 0
         startmodel[is.na(startmodel[, 3]), 3] <- 0
         startmodel[startmodel[, 3] == 0, 4] <- 0
         for(i in 1:nrow(startmodel)) {
            if((!is.na(startmodel[i, 4])) && startmodel[
               i, 3] != 0)
               startmodel[i, 4] <- 1
         }
         startmodel[is.na(startmodel[, 4]), 4] <- 0
         startmodelsize <- nrow(startmodel) + 1
      }
      else {
         startmodelsize <- startmodel$model.size
         startmodel <- startmodel$model[-1,  ]
         startknots1 <- startmodel$knot1
         startknots2 <- startmodel$knot2
         L1 <- FALSE
         if(!is.null(startmodel$level1)) {
            L1 <- TRUE
            level1 <- startmodel$level1
         }
         if(L1) {
            startmodel$knot1[!is.na(level1)] <- 1
            startknots1[!is.na(level1)] <- level1[!is.na(
               level1)]
         }
         startknots <- cbind(startknots1, startknots2)
         startknots <- as.vector(t(startknots))
         startknots[is.na(startknots)] <- 0.
         startmodel <- cbind(startmodel[, "pred1"], startmodel[
            , "knot1"], startmodel[, "pred2"], startmodel[
            , "knot2"])
         startmodel[, 2] <- !is.na(startmodel[, 2])
         startmodel[, 4] <- !is.na(startmodel[, 4])
      }
   }
   else {
      startmodel <- 0
      startknots <- 0
   }
   if(!ism8) {
      if(length(weights) != ncases) {
         stop("Number of weights not equal to the numnber of cases\n"
            )
      }
      weighted <- TRUE
      responses <- responses * weights
   }
   else {
      weighted <- FALSE
      weights <- 0
   }
   datamatrix <- cbind(responses, predictors)
   #Predictors which cannot interact together in the model are specified 
   #by a 2*n matrix of predictor indicies
   if(!ism9) {
      if(!is.matrix(no.interact) || ncol(no.interact) != 2) {
         stop("list of interactions disallowed has been misspecified,must be a 2*n matrix"
            )
      }
      no.interact <- t(no.interact)
      no.interact.size <- ncol(no.interact)
   }
   else {
      no.interact.size <- 0
      no.interact <- 0
   }
   if(startmodelsize > maxsize) {
      stop("start model should not be of greater size than the max model size\n"
         )
   }
   #Some error checking on the startmodel
   if(startmodelsize != 1) {
      for(i in 1:(startmodelsize - 1)) {
         if(startmodel[i, 1] == 0) {
            stop("first column of startmodel cannot be zero\n"
               )
         }
         if(startmodel[i, 2] == 1) {
            if(startknots[(i * 2) - 1] < min(predictors[
               , startmodel[i, 1]]) || startknots[
               (i * 2) - 1] > max(predictors[, 
               startmodel[i, 1]])) {
               stop("Knot out of range of its predictor \n"
                  )
            }
         }
         if(startmodel[i, 4] == 1) {
            if(startknots[(i * 2)] <= min(predictors[,
               startmodel[i, 3]]) || startknots[(
               i * 2)] >= max(predictors[, startmodel[
               i, 3]])) {
               stop("Knot out of range of its predictor\n"
                  )
            }
         }
      }
      if(max(startmodel[, c(1, 3)] > npredictors)) {
         stop("Initial model misspecified on input\n")
      }
   }
   startmodel <- t(startmodel)
   resultmodelsize <- 0
   end.state <- 0
   step.count <- 0
   z <- .C("polymars",
      as.integer(npredictors),    
      as.integer(nresponses),    
      as.integer(ncases),   
      as.double(datamatrix),
      as.integer(knots),
      as.double(mesh.vector),
      as.integer(mesh.specified),
      as.integer(maxsize),             
      as.double(gcv),                   
      as.integer(additive),              
      as.integer(startmodelsize),         
      start.model = as.integer(startmodel),
      start.knots = as.double(startknots),
      as.integer(weighted),
      as.double(weights),
      as.integer(no.interact.size),        
      as.integer(no.interact),
      as.integer(no.remove.size), 
      as.integer(no.remove),
      as.integer(knot.space),      
      as.integer(testset),          
      as.double(testsetmatrix),
      as.integer(testsetcases),
      as.integer(testset.weighted),  
      as.double(ts.weights),
      as.integer(classify),  
      as.double(tolerance),   
      as.integer(verbose),     
      best.model = as.integer(matrix(nrow = maxsize, ncol = 4, data
          = rep(0, maxsize * 4))),
      coefficients = as.double(matrix(nrow = maxsize, ncol = 
         nresponses, data = rep(0., maxsize * nresponses))),

         steps = as.integer(matrix(nrow = maxsize * 2, ncol = 2,
         data = rep(0, maxsize * 4))),
      rss.gcv = as.double(matrix(nrow = maxsize * 2, ncol = 
         nresponses + 1, data = rep(0., maxsize * 2 * (
         nresponses + 1)))),
      modelsize = as.integer(resultmodelsize),
      modelknots = as.double(matrix(nrow = maxsize, ncol = 2, data = 
         rep(0., maxsize * 2))),
      coefficient.se.term = as.double(rep(0., maxsize)),
      end.state = as.integer(end.state),
      step.count = as.integer(step.count),
      PACKAGE = "polspline")
   #The C function returns information about how it ended
   if(z$end.state != 0 && z$end.state != 5) {
      switch(z$end.state,
         stop("Mis-specification of initial model\n"),
         stop("Initial model with non-linear function must contain the corresponding linear function\n"
            ),
         stop("Initial model contains two-predictor functions that require prerequisite functions\n"
            ))
   }
   else {
      model <- matrix(z$best.model[1:((z$modelsize - 1) * 4)], ncol
          = 4, byrow = TRUE)
      knot.values <- matrix(z$modelknots[1:((z$modelsize - 1) * 2)],
         ncol = 2, byrow = TRUE)
      for(i in 1:nrow(model)) {
         if(model[i, 2] != 0) {
            model[i, 2] <- knot.values[i, 1]
         }
         else {
            model[i, 2] <- NA
         }
         if(model[i, 4] != 0) {
            model[i, 4] <- knot.values[i, 2]
         }
         else {
            model[i, 4] <- NA
         }
      }
      if(length(knots[model[, 1]]) != 0 && min(knots[model[, 1]]) <
         0) {
         factor1 <- TRUE
         levels1 <- rep(NA, z$modelsize - 1)
         factor.variables <- unique(model[knots[model[, 1]] <
            0, 1])
         for(i in 1:length(factor.variables)) {
            for(j in 1:length(model[, 1])) {
               if(model[j, 1] == factor.variables[
                  i]) {
                  levels1[j] <- model[j, 2]
               }
            }
            model[model[, 1] == factor.variables[i], 2] <-
               NA
         }
         levels1 <- c(NA, levels1)
      }
      else {
         factor1 <- FALSE
      }
      coefs <- matrix(z$coefficients[1:(z$modelsize * nresponses)],
         ncol = nresponses)
      #The model that the C-function returns does not explicitly contain an intercept
      #so in formatting the output one is added
      if(z$modelsize > 1) {
         if(factor1 == FALSE) {
            model <- rbind(c(0, NA, 0, NA), model)
            model <- data.frame(model, coefs)
            if(nresponses == 1) {
               dimnames(model) <- list(1:z$modelsize,
                  c("pred1", "knot1", "pred2",
                  "knot2", "coefs"))
            }
            else {
               dimnames(model) <- list(1:z$modelsize,
                  c("pred1", "knot1", "pred2",
                  "knot2", paste("Coefs", 1:
                  nresponses)))
            }
         }
         if(factor1 == TRUE) {
            model[(knots[model[, 1]] < 0), 2] <- NA
            model <- rbind(c(0, NA, 0, NA), model)
            model <- data.frame(model[, 1:2], levels1,
               model[, 3:4], coefs)
            if(nresponses == 1) {
               dimnames(model) <- list(1:z$modelsize,
                  c("pred1", "knot1", "level1",
                  "pred2", "knot2", "coefs"))
            }
            else {
               dimnames(model) <- list(1:z$modelsize,
                  c("pred1", "knot1", "level1",
                  "pred2", "knot2", paste("Coefs",
                  1:nresponses)))
            }
         }
      }
      else {
         model <- data.frame(0, NA, 0, NA, coefs)
         if(nresponses == 1) {
            dimnames(model) <- list(1:z$modelsize, c(
               "pred1", "knot1", "pred2", "knot2",
               "coefs"))
         }
         else {
            dimnames(model) <- list(1:z$modelsize, c(
               "pred1", "knot1", "pred2", "knot2",
               paste("Coefs", 1:nresponses)))
         }
      }
      #for later plotting the ranges and medians of the predictors are stored
      ranges.and.medians <- matrix(ncol = npredictors, nrow = 3,
         data = 0)
      for(i in 1:npredictors) {
         ranges.and.medians[1, i] <- min(predictors[, i])
      }
      for(i in 1:npredictors) {
         ranges.and.medians[2, i] <- max(predictors[, i])
      }
      for(i in 1:npredictors) {
         ranges.and.medians[3, i] <- median(predictors[, i])
      }
      # A table with information from the fitting is formatted here
      steps <- matrix(z$steps[1:(2 * (z$step.count + 1))], ncol = 2,
         byrow = TRUE)
      rss.gcv <- matrix(z$rss.gcv[1:((nresponses + 1) * (z$step.count +
         1))], ncol = nresponses + 1, byrow = TRUE)
      fitting <- data.frame(steps, rss.gcv)
      if(testset == FALSE) {
         if(nresponses == 1) {
            dimnames(fitting) <- list(1:(nrow(fitting)),
               c("0/1", "size", "RSS", "GCV"))
         }
         else {
            dimnames(fitting) <- list(1:nrow(fitting),
               c("0/1", "size", paste("RSS", 1:
               nresponses), "GCV"))
         }
      }
      else {
         if(classify == FALSE) {
            if(nresponses == 1) {
               dimnames(fitting) <- list(1:(nrow(
                  fitting)), c("0/1", "size",
                  "RSS", "T.S. RSS"))
            }
            else {
               dimnames(fitting) <- list(1:nrow(
                  fitting), c("0/1", "size",
                  paste("RSS", 1:nresponses),
                  "T.S. RSS"))
            }
         }
         else {
            if(nresponses == 1) {
               dimnames(fitting) <- list(1:(nrow(
                  fitting)), c("0/1", "size",
                  "RSS", "T.S.M.C."))
            }
            else {
               dimnames(fitting) <- list(1:nrow(
                  fitting), c("0/1", "size",
                  paste("RSS", 1:nresponses),
                  "T.S.M.C."))
            }
         }
      }
      # if their are factors present in the model the factors must be stored for use during plotting
      if(factor1 == TRUE) {
         model2 <- model[-1,  ]
         factors.in.model <- unique(model2[knots[model2[, 1]] <
            0, 1])
         maxfactors <- 0
         for(i in 1:length(factors.in.model)) {
            maxfactors <- max(maxfactors, length(unique(
               predictors[, factors.in.model[i]])))
         }
         factor.matrix <- matrix(ncol = length(factors.in.model),
            nrow = maxfactors + 2, data = NA)
         for(i in 1:length(factors.in.model)) {
            factor.matrix[1, i] <- factors.in.model[i]
            factor.matrix[2, i] <- length(unique(predictors[
               , factors.in.model[i]]))
            for(j in 3:(length(unique(predictors[, 
               factors.in.model[i]])) + 2)) {
               factor.matrix[j, i] <- unique(
                  predictors[, factors.in.model[
                  i]])[j - 2]
            }
         }
      }
      else {
         factor.matrix <- 0
      }
      if(nresponses == 1) {
         model <- cbind(model, model[,1])
         dimnames(model)[[2]][length(dimnames(model)[[2]])] <-
            "SE"
      }
      else {
         for(i in 1:nresponses) {
            model <- cbind(model, model[,1])
            dimnames(model)[[2]][length(dimnames(model)[[
               2]])] <- paste("SE", i)
         }
      }
      result <- list(model = model, fitting = fitting, model.size = z$
         modelsize,  responses = nresponses, ranges.and.medians = 
         ranges.and.medians, call = call, conversion = 
         conversion, factor.matrix = factor.matrix)
      class(result) <- "polymars"
      #refit the coefficients
      dd <- design.polymars(result,predictors)
      model2 <- result$model
      rsquared2 <- rep(0,nresponses)
      for(i in 1:nresponses){
        if(z$modelsize>1) mm <- summary(lm(responses[, i] ~ dd[, -1]))
        else mm <- summary(lm(responses[, i] ~ 1 ))
        rsquared2[i] <- mm$r.squared
        mm <- mm$coefficients
        model2[,i+factor1+4] <- mm[,1]
        model2[,i+factor1+4+nresponses] <- mm[,2]
      }
      result$model <- model2
      result$Rsquared <- rsquared2
      #calculates fitted values and residual of the data according to the
      #model returned 
      if(z$modelsize > 1) {
         fitted <- predict.polymars(result, x = predictors)
         residuals <- responses - fitted
      }
      else {
         fitted <- matrix(ncol = nresponses, nrow = ncases,
            data = coefs[1, 1])
         residuals <- matrix(ncol = nresponses, nrow = ncases,
            data = responses - coefs[1, 1])
      }
      result$residuals = residuals
      result$fitted = fitted
      return(result)
   }
}

################################################################################################

predict.polymars<-function(object,x,classify=FALSE,intercept,...)
{
 # Produces predicted values for a polymars object
 # pmars.model  an object returned from a call to polymars
 # x           a matrix with number of columns equal to number of columns of predictor matrix in 
 #             original call to polymars and predictor values in the corresponding columns. Can 
 #             also be a matrix with number of column equal to the number of predictors in the 
 #             model, in the order of the original dataset.
 # classify    If the original call to polymars was for classification setting  classify=TRUE will 
 #             the new data otherwise it will return the multi-response fitted values.
 # intercept   By default TRUE. The full intercept is included; or when FALSE the intercept is left out.
 #             Can also be given a numerical value
 
 if(missing(intercept))
  {
   intercept<-TRUE
  }
 if(!missing(x))x <- unstrip(x)
 # some error checking
                if(class(object)!="polymars")
         stop("object is not a polymars object")
  pmars.model <- object
 # The x matrix number of columns can be of length equal to the number of 
 # predictors in the original model or shorten to the number of predictors in 
 # the model in `pmars.model'
 if(!(is.matrix(x))) 
  {
   if(length(unique(pmars.model$model[, "pred1"]))== 1 ||  ncol(pmars.model$ranges.and.medians)== 1  )
    {
   x<-matrix(data=x,ncol=1)
    }
  }
 if((is.matrix(x) && ncol(x) 
    != length(unique(pmars.model$model[,"pred1"]))))
  {
   if(ncol(x) != ncol(pmars.model$ranges.and.medians))
    {
     
     
     stop("Input should be a matrix with number of columns equal to either number of original predictors or number of predictors in model\n")
     
    }
  }
 
 # If the number of columns of the matrix is not length equal to number of 
 # predictors it is expanded to that size.
 if(is.matrix(x) && ncol(x) == length(unique(pmars.model$model[, "pred1"])) && ncol(x) != ncol(pmars.model$ranges.and.medians))
  {
   tempmatrix<-x
   
   x<-matrix(nrow=nrow(tempmatrix),ncol=ncol(pmars.model$ranges.and.medians),data = 0)
   for(i in 1:length(unique(pmars.model$model[, "pred1"]))) 
    {
     for(j in 1:nrow(tempmatrix))
      {
       x[j,sort(unique(pmars.model$model[,"pred1"]))[i]]<-x[j]
      }
    }   
  }
 # If x is a vector put it into matrix form expanding it if it is of 
 # length equal to only the number of predictors in the model in `pmars.model'
 if(!(is.matrix(x))) 
  {
   if(!(length(x) == ncol(pmars.model$ranges.and.medians) || length(x) == unique(pmars.model$model[, "pred1"]))) 
    {
     stop("The vector of values must be equal in length to either the number of original predictors or predictors in the model\n")
    
    }
   if(length(x) == unique(pmars.model$model[, "pred1"]) && length(x) != ncol(pmars.model$ranges.and.medians)) 
    {
     x <- rep(0, ncol(pmars.model$ranges.and.medians))
     for(i in 1:length(unique(pmars.model$model[, "pred1"]))) 
      {
       x[sort(unique(pmars.model$model[, "pred1"]))[i]]<-x[i]
      }
    }
   x <- t(as.matrix(x))
  }
 
 # Checking to see if there are factor variables in the model
 if(dimnames(pmars.model$model)[[2]][3] == "level1")
  {
   level1<-TRUE
   pmars.model$model<-pmars.model$model[,c(1:(5+pmars.model$responses))]


   #if(dimnames(pmars.model$model)[[2]][6] == "level2"){level2<-TRUE}else{level2<-FALSE}
  }
 else
  {
   level1<-FALSE
   pmars.model$model<-pmars.model$model[,c(1:(4+pmars.model$responses))]
   #if(dimnames(pmars.model$model)[[2]][5] == "level2")
   # {level2<-TRUE}else{level2<-FALSE}
  }
 # Setting up the fitted responses matrix
 responses<-pmars.model$responses
 Y <- matrix(ncol = responses, nrow = nrow(x), data = rep(0, nrow(x)))
 Y1 <- matrix(ncol = 1, nrow = nrow(x), data = rep(0, nrow(x)))
 Y2 <- matrix(ncol = 1, nrow = nrow(x), data = rep(0, nrow(x)))
 if(is.logical(intercept))
  {
   if(intercept==TRUE)
    {
     for(i in 1:responses)Y[,i] <- pmars.model$model[1,ncol(pmars.model$model)-responses+i]
     }
   else
     {
      if(intercept==FALSE)
       {
        for(i in 1:responses)Y[,i] <- 0.0
        
       }
     }
  }
 else
  {
   if(is.numeric(intercept))
    {
     if(length(intercept)==responses)
      {
       for(i in 1:responses)Y[,i] <- intercept[i]
      }
     else
      {
       if(length(intercept) != 1)
        {
                stop("Intercept arguement mispecified \n")
                
        }
       for(i in 1:responses)Y[,i] <- intercept
      }

    }
  }
 # Computing fitted values
 if(pmars.model$model.size>1)
  {
   for(i in 2:pmars.model$model.size) 
    {
  
     Y2[] <- 1
   
     Y1[] <- x[,pmars.model$model[i, "pred1"]]
     if(!is.na(pmars.model$model[i, "knot1"])) 
      {
       Y1 <- Y1 - pmars.model$model[i,"knot1"]
       Y1[Y1 < 0,] <- 0
      }
     if(level1)
      {

     if(!is.na(pmars.model$model[i, "level1"]))
      {
       Y1<- (Y1 == pmars.model$model[i, "level1"])
      }
    }
   if(!is.na(pmars.model$model[i, "pred2"]) & pmars.model$model[i, "pred2"] != 0) 
    {
     Y2[] <- x[,pmars.model$model[i,"pred2"]]
     if(!is.na(pmars.model$model[i,"knot2" ])) 
      {
       Y2 <- Y2 - pmars.model$model[i,"knot2"]
       Y2[Y2 < 0,] <- 0
      }
     #if(level2)
     #{
     #  if(!is.na(pmars.model$model[i, "level2"]))
     #        {
     #   Y2<- (Y2 == pmars.model$model[i, "level2"])
     #  }
     #}
    }
   
    for(j in 1:responses){Y[,j]<-Y[,j]+(Y1 * Y2 * pmars.model$model[i,ncol(pmars.model$model)-responses+j])}
   }
  }
 # If classification is to be used the original polymars fitting expanded the 
 # response into a vector of indicator variables. The largest of the responses 
 # correspondes to the fitted class for each case.
 if(classify == TRUE)
  {
   for(i in 1:nrow(Y))
    {
     Y[i,]<-Y[i,]==max(Y[i,])
    }
   if(is.matrix(pmars.model$conversion))
   Z<-Y
   
   Y<-matrix(ncol=1,nrow=nrow(Z))
   for(i in 1:nrow(Y))
    {
     for(j in 1:ncol(Z))
      {
       
       if(Z[i,j] == 1) Y[i,] <- pmars.model$conversion[j]
      }     
    }
   }
# else
#  {
   # if classification was used but the full multiple response fitted response 
   # matrix is requested the response names (corresponding to the classes) are 
   # added.
  # if(is.matrix(pmars.model$conversion))
  #  {
  #   dimnames(Y)[[2]]<-list(pmars.model$conversion[,1])
#
#    }
#  }
        
 return(Y)
}
################################################################################################
print.polymars<-function(x,...)
{
        summary.polymars(x)
}
################################################################################################
summary.polymars<-function(object,...)
{
       if(class(object)!="polymars")
        stop("object is not a polymars object")
        pmars.model <- object
        cat("Call:\n")
        print(pmars.model$call)
        cat("\nModel fitting\n\n")
        print(pmars.model$fitting)
        cat("\n\nModel produced\n\n")
        print(pmars.model$model)
        if(pmars.model$responses != 1)
                cat("\nRESPONSES :", pmars.model$responses, "\n")
        if(!is.null(pmars.model$Rsquared))
                cat("\nRsquared :",round(pmars.model$Rsquared,3),"\n")
        invisible()
}


plot.polymars<-function(x,predictor1,response,predictor2,xx,add=FALSE,n,xyz=FALSE,contour.polymars=FALSE,xlim,ylim,intercept,...)
{
 # pmars.model      a polymars object
 # predictor1      the column number in the original predictor matrix of the predictor of interest
 # response        with multi-response polymars the  column number in the original response matrix.
 #                 Default is 1
 # predictor2      the second predictor for a contour or persp plot. For single response data 
 #                 plot(pmars1,2,6) is understood as 3-d plot of predictors 2 and 6.
 # xx              Values for the other predictors can be given, using the same format as for the 
 #                 predict fuhnction. By default median values are used.
 # add             should the plot be added to another. (for 2-d plots only)
 # n               For 2-d plot the number of points the function is interploted over. For 3-d plots
 #                 the a n*n mesh is interploted over. Default 2-d: 100, 3-d 33.
 # xyz             sometimes a call can be ambiguous: plot(pmars1,6,2) a 2-d plot with 2nd response 
 #                 or a 3-d plot. Use xyz=TRUE for 3-d.
 #contour.polymars By default a 3-d if a `persp' plot. contour.polymars=TRUE asks for a `contour' plot.
 #intercept        same as for predict function. =TRUE intercepr is included =FALSE it is left out, or can
 #                 be given a numerical value.

       if(class(x)!="polymars")
        stop("x is not a polymars object")
 pmars.model <- x
 if(missing(xx))xx<-pmars.model$ranges.and.medians[3,]
 if(length(xx) != ncol(pmars.model$ranges.and.medians))
  {
   stop("xx should be of length equal to the number of predictors in original data\n")
  
  }
  x <- xx
 if(!missing(predictor2))xyz <- TRUE
 if(missing(predictor2) && (!missing(response)) && pmars.model$responses == 1)
  {
   if(missing(predictor1) && xyz == TRUE)
    {

     stop("You must specify 2 predictor numbers")
     

    }
   xyz<-TRUE
   predictor2<-response
   response<-1
   
  } 
 if(contour.polymars == TRUE)
  {
   xyz<-TRUE
  }
 
 if(missing(intercept))
  {
   intercept<-TRUE
  }
 
 if(xyz==TRUE)
  {
   
   
   if(missing(n))n<-33
   if(missing(response))
    {
     
     if(missing(xlim))
      {
       persp.polymars(pmars.model,
                  predictor1,
                  predictor2,
                  n=n,
                  
                  contour.polymars=contour.polymars,
                  intercept=intercept,
                  ...)
       
      }
     else
      {
      
       persp.polymars(pmars.model,
                  predictor1,
                  predictor2,
                  n=n,
                  xlim=xlim,
                  
                  contour.polymars=contour.polymars,
                  intercept=intercept,
                  ...)
      }
    }
   else
    {
     
     if(missing(xlim))
      {
       persp.polymars(pmars.model,
                  predictor1,
                  predictor2,
                  response,
                  n=n,
                  
                  contour.polymars=contour.polymars,
                  intercept=intercept,
                  ...)
      }
     else
      {
      persp.polymars(pmars.model,
                 predictor1,
                 predictor2,
                 response,
                 n=n,
                 xlim=xlim,
                 contour.polymars=contour.polymars,
                 intercept=intercept,
                 ...)
      }
    }
   invisible(return())
   }
 else
  {
   
   if(missing(predictor1))
    {
     cat("predictor should be specified \n")
    }
  if(pmars.model$responses != 1 && missing(response)&& missing(predictor2))
   {
    cat("Response should be specified  (default: response =1)\n")
   }
  
  #check to see that the predictor is in the model
  inmodel<-FALSE
  for(i in 2:pmars.model$model.size)
   {
    if(pmars.model$model[i,"pred1"] == predictor1)inmodel<-TRUE
   }
  #check to see if the predictor is a factor

  if(is.matrix(pmars.model$factor.matrix))
   {
    if(length(pmars.model$factor.matrix[1,pmars.model$factor.matrix[1,]== predictor1]) != 0)
     {
      isfactor<-TRUE
     }
    else
     { 
      isfactor<-FALSE
     }
   }
  else
   {
    isfactor<-FALSE
   }

  if(isfactor == TRUE)
   { 
   

    pred.values <- matrix(nrow = pmars.model$factor.matrix[2,pmars.model$factor.matrix[1,]==predictor1], ncol = ncol(pmars.model$ranges.and.medians),data = x, byrow = TRUE)
    factors<-pmars.model$factor.matrix[-c(1,2),pmars.model$factor.matrix[1,]==predictor1]
    pred.values[,predictor1]<- factors[!is.na(factors)]

   
   
   mesh<-factors[!is.na(factors)]        
   }
  else
   {
    if(missing(n))n<-100
    if(missing(xlim))xlim<-c(pmars.model$ranges.and.medians[1,predictor1],pmars.model$ranges.and.medians[2,predictor1])
    
    pred.values <- matrix(nrow = n, ncol = ncol(pmars.model$ranges.and.medians), 
    data = x, byrow = TRUE)
    
    mesh <- matrix(seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/(n-1)),nrow=1)
    pred.values[,predictor1]<-mesh
   }
  if(missing(response) && missing(predictor2))response<-1
  if(missing(response))response <- 1
  if(response > pmars.model$responses || response < 0)
   {
    stop("response arguement = ",response,"is out of range\n")
    
   }

  model<-pmars.model$model

 Y<-predict.polymars(pmars.model,pred.values,intercept=intercept)
 
 if(isfactor == FALSE)
  {
   if(add == FALSE)
    {
     if(pmars.model$responses == 1)
      {
       plot(mesh,Y,...,type="l",xlab=paste("Predictor ",predictor1),ylab="Response")
        
      }
     else
      {
        
       plot(mesh,
            Y[,response],
            type="l",
            xlab=paste("Predictor ",predictor1),
            ylab=paste("Response ",response),
            ...)
      }
    }
  else
   {
        
    points(mesh,
           Y,
           type="l")
   }
  }

 if(isfactor == TRUE)
  {
         
   if(add == FALSE)
    {
     if(pmars.model$responses == 1)
      {
       plot(mesh,Y,...,xlab=paste("Predictor ",predictor1),ylab="Response")
     }  
    else
     { 
       plot(mesh,
           Y[,response],
           type="l",
           xlab=paste("Predictor ",predictor1),
           ylab=paste("Response ",response),
           ...)
     }
    } 
   else
    {
     points(mesh,
            Y,
            type="l")
    }
  }





        
  invisible()
 }
}


################################################################################################
persp.polymars<-function(x, predictor1, predictor2, response, n= 33,xlim,ylim,xx,contour.polymars,main,intercept,...)
{
       if(class(x)!="polymars")
        stop("x is not a polymars object")
  pmars.model <- x
 # used by the plot.polymars function
 # not designed for stand alone use.

 if(missing(xx))xx<-pmars.model$ranges.and.medians[3,]
 if(missing(xlim))xlim<-c(pmars.model$ranges.and.medians[1,predictor1],pmars.model$ranges.and.medians[2,predictor1])
 if(missing(ylim))ylim<-c(pmars.model$ranges.and.medians[1,predictor2],pmars.model$ranges.and.medians[2,predictor2])
 if(missing(predictor1) || missing(predictor2)) 
  {
   stop("You must specify 2 predictor numbers\n")
   
  }
 if(pmars.model$responses != 1 && missing(response)) 
  {
   cat("Response should be specified  (default: response =1)\n")
  }
 if(missing(response))response <- 1
 if(response > pmars.model$responses || response < 0)
  {
   stop("response arguement = ",response,"is out of range\n")
   
  }

 if(sum(as.integer(predictor1==pmars.model$model[,1])) == 0)
  {
   stop("Predictor 1 not in model\n")
   
  }
 if(sum(as.integer(predictor2==pmars.model$model[,1])) == 0)
  {
   stop("Predictor 2 not in model\n")
  
  }
 

 X <- seq(xlim[1],xlim[2],(xlim[2] - xlim[1])/(n-1))

 y <- seq(ylim[1],ylim[2],(ylim[2] - ylim[1])/(n-1)) 
 meshX <- rep(X, n)
 meshY <- rep(y, n)
 meshY <- sort(meshY)
 
 pred.values <- matrix(nrow = n^2, ncol = ncol(pmars.model$ranges.and.medians), 
 data = xx, byrow = TRUE)
 
 for(i in 1:(n^2))pred.values[i, predictor1] <- meshX[i]
 for(i in 1:(n^2))pred.values[i, predictor2] <- meshY[i]
 Z <- predict.polymars(pmars.model, pred.values,intercept=intercept)[,response]
 Z <- matrix(Z, ncol = n, byrow = FALSE)
 
 
 xtitle<-paste("Predictor", predictor1)
 ytitle<-paste("Predictor", predictor2)
 if(pmars.model$responses > 1) 
  {
 
   if(missing(main) && (!contour.polymars))
    {
     ztitle <- paste("Response", response)
    }
     
   if(missing(main) && (contour.polymars))
     {
      ztitle <- paste("Contour of response",response)
     }
  }
 else
  {
   if(missing(main) && (!contour.polymars))ztitle <- "Response"
   if(missing(main) && contour.polymars)ztitle <- paste("Contour of response")
  }
 if(!contour.polymars)
  {
   persp(X, y, Z, xlab = xtitle, ylab= ytitle, zlab = ztitle, ...)
  }
 else
  {
   
     contour(X, y, Z, xlab = xtitle, ylab = ytitle , main = ztitle, ...)
   
   }
 invisible()
}
################################################################################################
design.polymars<-function(object,x)
{
 # Produces predicted values for a polymars object
 # pmars.model  an object returned from a call to polymars
 # x           a matrix with number of columns equal to number of columns of predictor matrix in 
 #             original call to polymars and predictor values in the corresponding columns. Can 
 #             also be a matrix with number of column equal to the number of predictors in the 
 #             model, in the order of the original dataset.
 
 if(!missing(x))x <- unstrip(x)
 # some error checking
                if(class(object)!="polymars")
         stop("object is not a polymars object")
  pmars.model <- object
 # The x matrix number of columns can be of length equal to the number of 
 # predictors in the original model or shorten to the number of predictors in 
 # the model in `pmars.model'
 if(!(is.matrix(x))) 
  {
   if(length(unique(pmars.model$model[, "pred1"]))== 1 ||  ncol(pmars.model$ranges.and.medians)== 1  )
    {
   x<-matrix(data=x,ncol=1)
    }
  }
 if((is.matrix(x) && ncol(x) 
    != length(unique(pmars.model$model[,"pred1"]))))
  {
   if(ncol(x) != ncol(pmars.model$ranges.and.medians))
    {
     
     
     stop("Input should be a matrix with number of columns equal to either number of original predictors or number of predictors in model\n")
     
    }
  }
 
 # If the number of columns of the matrix is not length equal to number of 
 # predictors it is expanded to that size.
 if(is.matrix(x) && ncol(x) == length(unique(pmars.model$model[, "pred1"])) && ncol(x) != ncol(pmars.model$ranges.and.medians))
  {
   tempmatrix<-x
   
   x<-matrix(nrow=nrow(tempmatrix),ncol=ncol(pmars.model$ranges.and.medians),data = 0)
   for(i in 1:length(unique(pmars.model$model[, "pred1"]))) 
    {
     for(j in 1:nrow(tempmatrix))
      {
       x[j,sort(unique(pmars.model$model[,"pred1"]))[i]]<-x[j]
      }
    }   
  }
 # If x is a vector put it into matrix form expanding it if it is of 
 # length equal to only the number of predictors in the model in `pmars.model'
 if(!(is.matrix(x))) 
  {
   if(!(length(x) == ncol(pmars.model$ranges.and.medians) || length(x) == unique(pmars.model$model[, "pred1"]))) 
    {
     stop("The vector of values must be equal in length to either the number of original predictors or predictors in the model\n")
    
    }
   if(length(x) == unique(pmars.model$model[, "pred1"]) && length(x) != ncol(pmars.model$ranges.and.medians)) 
    {
     x <- rep(0, ncol(pmars.model$ranges.and.medians))
     for(i in 1:length(unique(pmars.model$model[, "pred1"]))) 
      {
       x[sort(unique(pmars.model$model[, "pred1"]))[i]]<-x[i]
      }
    }
   x <- t(as.matrix(x))
  }
 
 # Checking to see if there are factor variables in the model
 if(dimnames(pmars.model$model)[[2]][3] == "level1")
  {
   level1<-TRUE
   pmars.model$model<-pmars.model$model[,c(1:(5+pmars.model$responses))]


   #if(dimnames(pmars.model$model)[[2]][6] == "level2"){level2<-TRUE}else{level2<-FALSE}
  }
 else
  {
   level1<-FALSE
   pmars.model$model<-pmars.model$model[,c(1:(4+pmars.model$responses))]
   #if(dimnames(pmars.model$model)[[2]][5] == "level2")
   # {level2<-TRUE}else{level2<-FALSE}
  }
 # Setting up the fitted responses matrix
 responses<-pmars.model$responses
 Y <- matrix(ncol = 1, nrow = nrow(x), data = rep(1, nrow(x)))
 Y1 <- matrix(ncol = 1, nrow = nrow(x), data = rep(0, nrow(x)))
 Y2 <- matrix(ncol = 1, nrow = nrow(x), data = rep(0, nrow(x)))
 # Computing fitted values
 if(pmars.model$model.size>1)
  {
   for(i in 2:pmars.model$model.size) 
    {
  
     Y2[] <- 1
   
     Y1[] <- x[,pmars.model$model[i, "pred1"]]
     if(!is.na(pmars.model$model[i, "knot1"])) 
      {
       Y1 <- Y1 - pmars.model$model[i,"knot1"]
       Y1[Y1 < 0,] <- 0
      }
     if(level1)
      {

     if(!is.na(pmars.model$model[i, "level1"]))
      {
       Y1<- (Y1 == pmars.model$model[i, "level1"])
      }
    }
   if(!is.na(pmars.model$model[i, "pred2"]) & pmars.model$model[i, "pred2"] != 0) 
    {
     Y2[] <- x[,pmars.model$model[i,"pred2"]]
     if(!is.na(pmars.model$model[i,"knot2" ])) 
      {
       Y2 <- Y2 - pmars.model$model[i,"knot2"]
       Y2[Y2 < 0,] <- 0
      }
     #if(level2)
     #{
     #  if(!is.na(pmars.model$model[i, "level2"]))
     #        {
     #   Y2<- (Y2 == pmars.model$model[i, "level2"])
     #  }
     #}
    }
   
    Y<-cbind(Y,Y1 * Y2)
   }
  }
 return(Y)
}
################################################################################################

logspline <- function(x, lbound, ubound, maxknots=0, knots, nknots=0,
   penalty= -1, silent = TRUE,mind= -1, error.action=2)
{
   call <- match.call()
   if(!missing(x))x <- unstrip(x)
   data <- x
   ilx <- 0; iux <- 0
   if(!missing(lbound)){ilx <- 1;jlx <- lbound}
   if(!missing(ubound)){iux <- 1;jux <- ubound}


   u2 <- length(data)
   data <- data[!is.na(data)]
   nsample <- length(data)
   if(nsample<10)stop("not enough data")
   if(u2 !=nsample) print(paste("***", u2-nsample, " NAs ignored in data"))
   data <- sort(data)

   # data can not be beyond the boundaries of the density
   if(!missing(lbound)) if(data[1] < lbound) stop("data below lbound")
   if(!missing(ubound)) if(data[nsample] > ubound) stop("data above ubound")
   mm <- range(data)
   if(!missing(lbound)) mm <- range(c(mm, lbound))
   if(!missing(ubound)) mm <- range(c(mm, ubound))

   # boundaries
   ilow <- (!missing(lbound)) * 1
   iupp <- (!missing(ubound)) * 1
   low <- 0
   upp <- 0
   if(ilow == 1) low <- lbound
   if(iupp == 1) upp <- ubound

   # get the maximal dimension
   intpars <- c(-100, rep(0, 9))
   z <- .C("nlogcensorx", z = as.integer(intpars),
      PACKAGE = "polspline")
   maxp <- z$z[1]

   # organize knots
   kts <- vector(mode = "double", length = max(maxp))
   if(maxknots > maxp - 5) warning(paste("maxknots reduced to", maxp))
   nknots <- -nknots
   if(!missing(knots)) {
      nknots <- length(knots)
      knots <- sort(knots)
      if(!missing(lbound)) if(min(knots) < lbound)
         stop("data (knots) below lbound")
      if(!missing(ubound)) if(max(knots) > ubound)
         stop("data (knots) above ubound")
      if(nknots < 3) stop("need at least three starting knots")
      if(nknots > maxp - 5) stop(paste("at most", maxp - 5, "knots possible"))
      kts[1:nknots] <- knots
   }

   silent <- (silent == FALSE)

   # group parameters
   intpars <- c(nsample, maxknots, nknots, silent, 1-ilow, 1-iupp,mind)
   dpars <- c(penalty, low, upp)
   data <- c(data, rep(0, maxp))

   # do it
   z <- .C("nlogcensor",
      ip = as.integer(intpars),
      coef = as.double(data),
      dp = as.double(dpars),
      logl = as.double(rep(0, maxp)),
      ad = as.integer(rep(0, maxp)),
      kts = as.double(kts),
      PACKAGE = "polspline")

   # error messages
   if(z$ip[1] != 0 && z$ip[1]<100) {
      if(z$ip[1] == 17) warning("too many knots beyond data")
      if(z$ip[1] == 18) warning("too many knots before data")
      if(z$ip[1] == 39) warning("too much data close together")
      if(z$ip[1] == 40) warning("no model could be fitted")
      if(z$ip[1] == 2) warning("error while solving system")
      if(z$ip[1] == 8) warning("too much step-halving")
      if(z$ip[1] == 5) warning("too much step-halving")
      if(z$ip[1] == 7) 
         warning("numerical problems, likely tail related. Try lbound/ubound")
      if(z$ip[1] == 1) warning("no convergence")
      i <- 0
      if(missing(knots))i<- 1     
      if(z$ip[1] == 3 && i==1) 
        warning("right tail extremely heavy, try running with ubound")
      if(z$ip[1] == 4 && i==1) 
        warning("left tail extremely heavy, try running with lbound")
      if(z$ip[1] == 6 && i==1) 
        warning("both tails extremely heavy, try running with lbound and ubound")
      if(z$ip[1] == 3 && i==0) 
        warning("right tail too heavy or not enough knots in right tail")
      if(z$ip[1] == 4 && i==0) 
        warning("left tail too heavy or not enough knots in left tail")
      if(z$ip[1] == 6 && i==0) 
        warning("both tails too heavy or not enough knots in both tail")
      if(error.action==0) stop("fatal error")
      if(error.action==1) {
         print("no object returned")
         invisible()
      }
      if(error.action==2) {
          if(ilx==0 && iux==0)z <- oldlogspline(x)
          if(ilx==0 && iux==1)z <- oldlogspline(x,ubound=jux)
          if(ilx==1 && iux==0)z <- oldlogspline(x,lbound=jlx)
          if(ilx==1 && iux==1)z <- oldlogspline(x,lbound=jlx,ubound=jux)
          z <- oldlogspline.to.logspline(z,x)
          z$call <- call
          warning("re-ran with oldlogspline")
          z
      }
   }
   else{
   if(z$ip[1]>100) {
      warning(" Not all models could be fitted")
   }
   # organize logl
   logl <- cbind(z$ad, z$logl)
   logl <- cbind(2+(1:z$ip[3]),logl[1+(1:z$ip[3]),  ])
   kk <- (1:length(logl[,1]))
   kk <- kk[logl[, 2] == 0 ]
   if(length(kk)>0)logl <- logl[-kk,]
   # bye bye
   fit <- list(call = call, nknots = z$ip[2], coef.pol = z$coef[1:2], coef.kts = 
      z$coef[2 + (1:z$ip[2])], knots = z$kts[1:z$ip[2]], maxknots = z$ip[3]+2,
      penalty = z$dp[1], bound = c(ilow, low, iupp, upp), samples = nsample,
      logl = logl, range = mm, mind = z$ip[7])
   class(fit) <- "logspline"
   fit}
}
plogspline <- function(q, fit)
{
    if(class(fit)!="logspline")
       stop("fit is not a logspline object")
   if(!missing(q))q <- unstrip(q)
    sq <- rank(q)
    q <- sort(q)
    z <- .C("rpqlsd",
        as.double(c(fit$coef.pol, fit$coef.kts)),
        as.double(fit$knots),
        as.double(fit$bound),
        as.integer(1),
        pp = as.double(q),
        as.integer(length(fit$knots)),
        as.integer(length(q)),
      PACKAGE = "polspline")
    zz <- z$pp[sq]
    if(fit$bound[1] > 0) zz[q<fit$bound[2]] <- 0
    if(fit$bound[3] > 0) zz[q>fit$bound[4]] <- 1
    zz
}
qlogspline <- function(p, fit)
{
    if(class(fit)!="logspline")
       stop("fit is not a logspline object")
   if(!missing(p))p <- unstrip(p)
    sp <- rank(p)
    p <- sort(p)
    z <- .C("rpqlsd",
        as.double(c(fit$coef.pol, fit$coef.kts)),
        as.double(fit$knots),
        as.double(fit$bound),
        as.integer(0),
        qq = as.double(p),
        as.integer(length(fit$knots)),
        as.integer(length(p)),
      PACKAGE = "polspline")
    zz <- z$qq[sp]
    zz[p<0] <- NA
    zz[p>1] <- NA
    zz
}
rlogspline <- function(n, fit)
{
    if(class(fit)!="logspline")
       stop("fit is not a logspline object")
    pp <- runif(n)
    qlogspline(pp, fit)
}
dlogspline <- function(q, fit)
{
    if(class(fit)!="logspline")
       stop("fit is not a logspline object")
   if(!missing(q))q <- unstrip(q)
    x <- q
    y <- fit$coef.pol[1] + x * fit$coef.pol[2]
    for(i in 1:length(fit$knots)) 
       y <- y + fit$coef.kts[i] * ((abs(x - fit$knots[i]) +x- fit$knots[i])/2)^3
    y <- exp(y)
    if(fit$bound[1] > 0) y[x < fit$bound[2]] <- 0
    if(fit$bound[3] > 0) y[x > fit$bound[4]] <- 0
    y
}
plot.logspline <-function(x, n = 100, what = "d", add = FALSE, xlim, xlab = "", ylab = "", type = "l", ...)
{
        fit <- x
    if(class(fit)!="logspline")
       stop("fit is not a logspline object")
        if(add){
                plim <- (par()$usr)[1:2]
                u4 <- plim[1]
                u3 <- plim[2]
                if(!missing(xlim)) {
                        u4 <- max(xlim[1], plim[1])
                        u3 <- min(xlim[2], plim[2])
                }
        }
        else{
        if(missing(xlim)) {
                u1 <- qlogspline(0.01, fit)
                u2 <- qlogspline(0.99, fit)
                u3 <- 1.1 * u1 - 0.1 * u2
                u4 <- 1.1 * u2 - 0.1 * u1
        }
        else {
                u3 <- xlim[1]
                u4 <- xlim[2]
        }}
        xx <- (0:(n - 1))/(n - 1) * (u4 - u3) + u3
        if(what == "d" || what == "D") yy <- dlogspline(xx, fit)
        if(what == "f" || what == "F" || what == "p" || what == "P")
                yy <- plogspline(xx, fit)
        if(what == "s" || what == "S") yy <- 1 - plogspline(xx, fit)
        if(what == "h" || what == "H") yy <- dlogspline(xx, fit)/(1 - plogspline(xx, fit))
        if(missing(xlab)) xlab <- ""
        if(missing(ylab)) ylab <- ""
        if(missing(type)) type <- "l"
        if(add)lines(xx,yy, ...)
        else plot(xx, yy, xlab = xlab, ylab = ylab, type = type, ...)
        invisible()
}
print.logspline <- function(x,...)
{
      summary.logspline(x)
}
summary.logspline <- function(object,...)
{
        fit <- object
    if(class(fit)!="logspline")
       stop("fit is not a logspline object")
   ul <- fit$penalty
   um <- fit$samples[1]
   if(length(fit$samples)>1)
   um <- fit$samples[1]+ fit$samples[4]
   else
   um <- fit$samples
   kk <- fit$logl[fit$logl[,2] != 0,1]
   ad <- fit$logl[fit$logl[,2] != 0,2]
   ll <- fit$logl[fit$logl[,2] != 0,3]
   bb <- -2 * ll + ul * (kk-1)
   cc1 <- bb
   cc2 <- bb
   cc2[1] <- Inf
   cc1[length(bb)] <- 0
   if(length(bb) > 1) {
      for(i in 1:(length(bb) - 1)) {
         cc1[i] <- max((ll[(i + 1):(length(bb))] - ll[i])/(kk[(i + 1):
                (length(bb))] - kk[i]))
         cc2[i + 1] <- min((ll[1:i] - ll[i + 1])/(kk[1:i] - kk[i + 1]))
      }
   }
   c3 <- cc2 - cc1
   cc1[c3 < 0] <- NA
   cc2[c3 < 0] <- NA
   uu <- cbind(kk, ad, ll, bb, 2 * cc1, 2 * cc2)
   ww <- rep("", length(bb))
   dimnames(uu) <- list(ww, c("knots", "A(1)/D(2)", "loglik", "AIC",
      "minimum penalty", "maximum penalty"))
   print(round(uu, 2))
   cat(paste("the present optimal number of knots is ",kk[bb== min(bb)],"\n"))
   if(ul == log(um))
      cat(paste("penalty(AIC) was the default: BIC=log(samplesize): log(",
         um, ")=", round(ul, 2), "\n"))
   else cat(paste("penalty(AIC) was ", round(ul, 2),
         ", the default (BIC) ", "would have been", round(log(um), 2), "\n"))
   invisible()
}
polyclass <- function(data, cov, weight, penalty, maxdim, exclude, include,
   additive = FALSE, linear, delete=2, fit, silent = TRUE, normweight = TRUE, tdata, 
   tcov, tweight, cv, select=0, loss, seed)
{
   call <- match.call()
   if(!missing(cov))cov <- unstrip(cov)
   if(!missing(exclude))exclude <- unstrip(exclude)
   if(!missing(include))include <- unstrip(include)
   if(!missing(data))data <- unstrip(data)
   if(!missing(weight))weight <- unstrip(weight)
   if(!missing(tdata))tdata <- unstrip(tdata)
   if(!missing(tweight))tweight <- unstrip(tweight)
   if(!missing(tcov))tcov <- unstrip(tcov)
   it <- 0
   ntdata <- 0
   if(!missing(cv)) it <- 2
   if(!missing(tdata))it <- 1
        if(missing(cv)) cv <- 0
        if(it==1||it==0) cv <- 0
   if(it==2){
      if(!missing(seed)){
         if(sum(seed)!=0){
            if(length(seed)>11) assign(".Random.seed",  seed[1:12], envir=.GlobalEnv)
            else set.seed(seed[1])
            seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
         }
      }
      else{
         if(!missing(fit)){
            if(fit$method==2) assign(".Random.seed",  fit$seef, envir=.GlobalEnv)
         }
         seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
      }
   }
   z <- .C("spolyx", mk = as.integer(rep(-3,13)),
      PACKAGE = "polspline")
   MAXKNOTS <- z$mk[1]
   MAXSPACE <- z$mk[2]    
   if(missing(data)) stop("there has to be data")
   if(length(data) < 25) stop("not enough data")
   if(is.integer(data) == FALSE){
      if(max(abs(as.integer(data) - data)) < 0.001)
         data <- as.integer(data)
      else stop("not-integer data")
   }
   if(it == 1) {
      if(is.integer(tdata) == FALSE){
         if(max(abs(as.integer(tdata) - tdata)) < 0.001)
            tdata <- as.integer(tdata)
         else stop("not-integer test data")
      }
      alldata <- c(data,tdata)
      if(min(alldata)<0) stop("negative data")
      clss <- min(alldata):max(alldata)
      if(min(alldata) == 1){
         data <- data - 1
         tdata <- tdata - 1
      }
      ntdata <- length(tdata)
      if(missing(tweight)) tweight <- rep(1,ntdata)
      if(length(tweight)!=ntdata)stop("length tweight is incorrect")
      if(normweight == TRUE)tweight <- tweight*ntdata/sum(tweight)
   }
   else{
      if(min(data)<0) stop("negative data")
      clss <- min(data):max(data)
      if(min(data) == 1) data <- data - 1
   }
   nclass <- length(clss)
   ndata <- length(data)   
   nu <-  exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
   if(nu) xx <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
   yy <- sample(ndata)
   if(nu)assign(".Random.seed",  xx, envir=.GlobalEnv)
   if(missing(weight)) weight <- rep(1,ndata)
   if(it==2){
      if(sum(abs(seed[1]))==0) myord <- 1:ndata
      else myord <- sample(ndata)
      data <- data[myord]
      weight <- weight[myord]
   }
   if(length(weight)!=ndata)stop("length weight is incorrect")
   if(normweight == TRUE)weight <- weight*ndata/sum(weight)
   if(missing(cov)) {
      stop("covariates required")
   }
   else {
      if(length(cov) == ndata)
         cov <- matrix(cov, ncol = 1, nrow = ndata)
      if(length(cov[, 1]) != ndata)
         stop("covariates not ndata * ncov matrix")
      if(it==2)cov <- cov[myord,]
      ncov <- length(cov[1,  ])
      nms <- 1:ncov
      if(is.matrix(cov))
         nms <- dimnames(cov)[[2]]
      if(length(nms) != ncov)
         nms <- 1:ncov
   }
   if(missing(penalty) && it ==0)
      penalty <- log(ndata)
   if(missing(penalty) && it >0)
      penalty <- 0
   il <- 1
   if(select==1) il <- 0
   if(select==2) il <- 2
   if(delete!=0 && delete !=1) delete <- 2
   iml <- missing(loss)
   if(iml) loss <- 1 - diag(rep(1,nclass))
   if(il!=1 && !iml)
      stop("if loss is specified, select has to be 0")
   if((it == 0) && !iml)
      warning("loss only has effect when there is a test-set or CV is used")
   if(it == 1){
      if(missing(tcov)) {
         if(ncov!=0)stop("missing tcov")
         tcov <- 0
      }
      else {
         if(length(tcov) == ntdata)
         tcov <- matrix(tcov, ncol = 1, nrow = ntdata)
         if(length(tcov[, 1]) != ntdata)
         stop("test-covariates not ntdata * ncov matrix")
         ntcov <- length(cov[1,  ])
         if(ntcov!=ncov) stop("wrong number of test-covariates")
      }
   }
   naction <- nclass
   if(it>0){
      if(is.matrix(loss)==FALSE)stop("loss is not a matrix")
      if(length(loss[1,])!=nclass)stop("loss has not nclass columns")
      naction <- length(loss[,1])
   }
   if(additive) {
      if(!missing(exclude)) stop("cannot have exclude and additive")
      if(!missing(include)) stop("cannot have include and additive")
      include <- c(0, 0)
   }
   if(missing(exclude) + missing(include) == 0)
      stop("only 1 from exclude and include allowed")
   vexclude <- 0   
   if(missing(exclude) == FALSE) {
      if(length(exclude) == 2)
         exclude <- matrix(exclude, ncol = 2, nrow = 1)
      if(length(exclude[1,  ]) != 2) stop("exclude has wrong shape")
      if(min(exclude) < 0 || max(exclude) > ncov)
         stop("exclude has wrong values")
      vexclude <- as.vector(t(exclude))
      vexclude <- c(length(vexclude)/2, vexclude)   
   }
   if(missing(include) == FALSE || additive) {
      if(length(include) == 2)
         include <- matrix(include, ncol = 2, nrow = 1)
      if(length(include[1,  ]) != 2)
         stop("include has wrong shape")
      if(min(include) < 0 || max(include) > ncov)
         stop("include has wrong values")
      include <- t(apply(include, 1, sort))
      if(length(include) == 2)
         include <- matrix(include, ncol = 2, nrow = 1)
      vexclude <- as.vector(t(include))
      vexclude <- c( - length(vexclude)/2, vexclude)
   }
   if(missing(maxdim)) {
      maxdim <- floor(4 * (ndata)^(1/3))+1
      maxdim <- min(ndata/2, MAXSPACE-1, (nclass-1)*maxdim)
      maxdim <-  - maxdim
   }
   if(maxdim > MAXSPACE - 1) {
      maxdim <- MAXSPACE - 1
      print(paste("maximum dimension reduced to", maxdim))
   }
   lins <- rep(0, MAXSPACE)
   if(!missing(linear)) {
      linear[linear <= 0] <- ncov + 1
      linear[linear > ncov + 1] <- ncov + 1
      lins[linear] <- 1
   }
   if(additive)
      vexclude <- c(-1, 0, 0)   # do it
   fitter <- 0
   bbtt <- matrix(0, ncol = 4 + max(data), nrow = abs(maxdim))
   cckk <- matrix(0, ncol = (MAXKNOTS + 1), nrow = ncov+1)
   if(!missing(fit)) {
      if(class(fit)!="polyclass")stop("fit is not a polyclass object")
      fitter <- (fit$nclass-1)*(fit$nbas)
      if(fit$ncov != ncov)
         stop("ncov and fit's ncov are different")
      if(fit$nclass != nclass)
         stop("nclass and fit's nclass are different")
                a1 <- length(fit$fcts[1,])
      bbtt[1:fit$nbas,  ] <- fit$fcts[,-a1]
      bbtt <- as.vector(t(bbtt))
      bbtt[is.na(bbtt)] <- -1
      a1 <- length(fit$knots[1,  ])
                a2 <- as.vector(t(fit$knots))
      cckk <- as.vector(cckk)
                cckk <- c(a1,a2,cckk)
      cckk[is.na(cckk)] <- -1
   }
   mindist <- 3*nclass
        if(missing(tdata)){
           tdata<-0
      tcov <-0
           tweight <- 0
        }
   ranges <- NA
   if(ncov == 1)
      ranges <- matrix(range(cov), ncol = 1, nrow = 2)
   if(ncov > 1)
      ranges <- apply(cov, 2, range)   # done
   cov <- as.single(t(cov))
   aicx <- as.single(rep(0,1000))
   intpars <-c(ndata,nclass,ncov,mindist,maxdim,silent,fitter,cv,it,ntdata,
             naction,il,delete)
   anova <- loss
   if(length(anova)<abs(maxdim)*4)anova<-c(anova,rep(0,abs(4*maxdim)))
   z <- .C("spoly",
      intpars = as.integer(intpars),
      as.integer(data),
      as.single(cov),
      anova = as.double(anova),
      as.double(penalty),
      bbtt = as.double(bbtt),
      cckk = as.double(cckk),
      as.integer(vexclude),
      as.integer(lins),
      logl = as.double(rep(0, 11*MAXSPACE+1)),
      as.double(weight),
      as.integer(tdata),
      as.single(t(tcov)),
      as.double(tweight),
      bbb = as.double(rep(0, MAXSPACE*nclass)),
      aicx=as.single(aicx),
      PACKAGE = "polspline")
   ndim <- z$intpars[1]
   aicx <- z$aicx[1:4]
   aicy <- 0
   if(it==2){
      aicy <- z$aicx[6:(z$aicx[5])]
      aicy <- matrix(aicy,ncol=3,byrow=TRUE)
      if(z$aicx[5]<995)aicy[length(aicy[,1]),2]<- Inf
      dimnames(aicy) <- list(NULL,c("pen-min","pen-max","cv-loss"))
   }
   nclass <- z$intpars[2]+1
   nbas <- z$intpars[3]
   maxdim <- abs(maxdim)
   z$bbtt <- matrix(z$bbtt, nrow = maxdim, ncol = 3 + nclass, byrow =  TRUE)
   z$bbtt <- z$bbtt[1:nbas,  ]
   z$cckk <- matrix(z$cckk, nrow = ncov+1, ncol = MAXKNOTS + 1, byrow = TRUE)
   z$cckk <- z$cckk[1:ncov,]
   z$cckk <- matrix(z$cckk, nrow = ncov)
   z$cckk <- z$cckk[, 1:(1 + max(z$cckk[, 1]))]
   z$cckk <- matrix(z$cckk, nrow = ncov)
   l1 <- max(z$cckk[, 1])
   for(i in 1:(ncov))
      if(z$cckk[i, 1] != l1) z$cckk[i, (z$cckk[i, 1] + 2):(l1 + 1)] <- 
            NA
   if(l1 > 0)
      dimnames(z$cckk) <- list(nms, c("K", 1:l1))
   if(l1 == 0)
      dimnames(z$cckk) <- list(nms, "K")
   z$bbtt <- matrix(z$bbtt, ncol = 3 + nclass)
   z$bbtt <- cbind(z$bbtt, 0)
   dimnames(z$bbtt) <- list(1:nbas, c("dim1", "knot1", 
      "dim2", "knot2", as.character(clss)))
   z$bbtt[z$bbtt[, 3] == -1, 3:4] <- NA
   z$bbtt[z$bbtt[, 4] == 0, 4] <- NA
   z$bbtt[1,1] <- NA
   i <- z$logl[1]
   z$logl <- matrix(z$logl[2:(11*i+1)],ncol=11,byrow=TRUE)
   z$logl[z$logl[,10]<0,10] <- NA
   z$logl[z$logl[,11]<0,11] <- NA
   z$logl[1,11] <- Inf
   dimnames(z$logl) <- list(NULL, c("dim","loss","l-lik-trn","loss-trn",
      "sq-err-trn","l-lik-test","loss-tst","sq-err-tst"
              , "A/D","pen-min","pen-max"))
   if(it!=1){
      dimnames(z$logl)[[2]][2] <- "AIC"
      z$logl <- z$logl[,-(6:8)]
   }
   anova <- z$anova[2:(1+z$anova[1])]
   anova[anova<0] <- NA
   anova <- matrix(anova,ncol=3,byrow=TRUE)
   dimnames(ranges) <- list(c("min", "max"), nms)
   z$bbtt[0, 0] <- NA
   z$bbtt[z$bbtt[, 2] == 0, 2] <- NA
   z$bbtt[z$bbtt[, 2] == 0, 4] <- NA
   if(nclass==naction)
      yyy <- clss
   else
      yyy <- 1:naction
   if(it!=0)dimnames(loss) <- list(as.character(yyy),clss)
   if(il!=1)loss <- -1
   bbb <- z$bbb
   bbb <- bbb[1:(nbas*nclass)]
   bbb <- matrix(bbb,nrow=nbas,byrow=TRUE)
   if(it==0){
      nfit <- list(call = call, ncov = ncov, ndim = ndim, nclass = nclass,
      nbas = nbas, fcts = z$bbtt, knots = z$cckk, penalty = penalty,
      method = it, ranges = ranges, logl= z$logl, 
      sample = ndata, wgtsum = sum(weight), covnames = nms,
      classnames = clss, beta = bbb, delete = delete, anova = anova)
   }
   else{
      if(it==1)
         nfit <- list(call = call, ncov = ncov, ndim = ndim, nclass = nclass,
         nbas = nbas, naction = naction, fcts = z$bbtt, knots = z$cckk,
         loss = loss, penalty = penalty, method = it, ranges = ranges,
         logl= z$logl, sample = ndata, tsample = ntdata,
         wgtsum = sum(weight), covnames = nms, classnames = clss, beta = bbb,
         delete = delete, anova = anova,
         select = select, twgtsum = sum(tweight))
      else
         nfit <- list(call = call, ncov = ncov, ndim = ndim, nclass = nclass,
         nbas = nbas, naction = naction, fcts = z$bbtt, knots = z$cckk,
         cv = cv, loss = loss, penalty = penalty, method = it, ranges = ranges,
         logl= z$logl, sample = ndata, wgtsum = sum(weight), covnames = nms,
         classnames = clss, cv.aic = aicx, cv.tab = aicy, seed = seed, 
         beta = bbb, delete = delete, anova = anova, select = select)
   }
   class(nfit) <- "polyclass"
   nfit
}
cpolyclass <- function(cov, fit)
{
      if(class(fit)!="polyclass")stop("fit is not a polyclass object")
      if(!missing(cov))cov <- unstrip(cov)
   xxx <- ppolyclass(cov, fit)
   yyy <- fit$classnames
   if(fit$method!=0){
   if(length(fit$loss)==1)
          fit$loss <- 1 - diag(rep(1,fit$nclass))
   xxx <- t(-fit$loss%*%t(xxx))
   if(fit$nclass==fit$naction)
      yyy <- fit$classnames
   else
      yyy <- 1:fit$naction
   }
   zzz <- xxx[, 1]
   www <- rep(yyy[1], length(zzz))
   for(i in 2:length(yyy)) {
      www[zzz < xxx[, i]] <- yyy[i]
      zzz[zzz < xxx[, i]] <- xxx[zzz < xxx[, i], i]
   }
   www
}
ppolyclass <- function(data, cov, fit)
{
   imf <- missing(fit)
   if(imf) {
      fit <- cov
      cov <- data
   }
      if(!missing(cov))cov <- unstrip(cov)
      if(!missing(data))data <- unstrip(data)
      if(class(fit)!="polyclass")stop("fit is not a polyclass object")
   if(is.matrix(cov) == FALSE)
      cov <- matrix(cov, ncol = 1)
   if(length(cov[1,  ]) != fit$ncov) {
      if(length(cov[1,  ]) == 1 && length(cov[, 1]) == fit$ncov)
         cov <- t(cov)
      else stop("incorrect number of covariates")
   }
   ncase <- length(cov[, 1])
   nclass <- fit$nclass
   nbas <- length(fit$fcts[, 1])
   if(imf || missing(data))
      data <- rep(-1, ncase)
   if(length(data) == 1)
      data <- rep(data, ncase)
   if(is.integer(data) == FALSE)
      if(max(abs(as.integer(data) - data)) < 0.001)
         data <- as.integer(data)
      else stop("not-integer data")
   w2 <- fit$classnames
   if(data[1] != -1 && (min(w2) > min(data) || max(w2) < max(data)))
      stop("data has wrong range")
   if(min(data) == 0)
      data <- data + 1
   ppp <- matrix(0, ncol = nclass, nrow = ncase)
   for(i in 1:(nclass - 1))
      ppp[, i] <- (fit$fcts[1, (4 + i)])
   if(nbas > 1)
      for(j in 2:nbas) {
         uuu <- cov[, fit$fcts[j, 1]]
         if(is.na(fit$fcts[j, 2]) == FALSE) {
            uuu <- uuu - fit$knots[fit$fcts[j, 1], fit$fcts[
              j, 2] + 1]
            uuu[uuu < 0] <- 0
         }
         vvv <- rep(1, ncase)
         if(is.na(fit$fcts[j, 3]) == FALSE) {
            vvv <- cov[, fit$fcts[j, 3]]
            if(is.na(fit$fcts[j, 4]) == FALSE) {
              vvv <- vvv - fit$knots[fit$fcts[j, 3], fit$
                fcts[j, 4] + 1]
              vvv[vvv < 0] <- 0
            }
         }
         uuu <- uuu * vvv
         for(i in 1:(nclass - 1))
            ppp[, i] <- ppp[, i] + uuu * fit$fcts[j, (4 + i
              )]
      }
   ppp <- ppp-apply(ppp,1,max)
   ppp <- exp(ppp)
   zzz <- ppp[, nclass]
   for(i in 1:(nclass - 1))
      zzz <- zzz + ppp[, i]
   for(i in 1:nclass)
      ppp[, i] <- ppp[, i]/zzz
   if(data[1] == -1)
      dimnames(ppp) <- list(NULL, fit$classnames)
   else ppp <- ppp[cbind(1:ncase, data)]
   ppp
}
plot.polyclass <- function(x,cov,  which, lims, what, data, n, xlab, ylab, zlab, ...)
{
      if(class(x)!="polyclass")stop("x is not a polyclass object")
      if(!missing(cov))cov <- unstrip(cov)
     fit <- x
   here <- c(-1, -1)
   if(length(which) == 1 || length(which) == 2)
      here[1] <- length(which)
   here[2] <- as.integer(what)
   if(here[2]< 1||here[2]>8) stop("what is wrong")
   if(min(here) < 0)
      stop("which is wrong")
   if(here[2] < 5.5 && here[1] == 1)
      stop("which and what contradict")
   if(here[2] > 5.5 && here[1] == 2)
      stop("which and what contradict")
   if(length(cov) != fit$ncov)
      stop("length of cov is wrong")
   clbs <- fit$covnames
   ww <- fit$classnames
   w1 <- (1:fit$ncov)
   if(missing(lims))
      lims <- NULL
   if(length(lims) != 0 && length(lims) != (here[1] * 2))
      stop("lims is wrong")
   wa <- 0
   for(i in 1:length(which)){
   if(is.numeric(which) == FALSE)
      wa <- c(wa,w1[which[i] == clbs])
   else wa <- c(wa,w1[w1 == which[i]])
   }
   wa <- wa[-1]
   if(length(wa) != here[1])
      stop("which is wrong")
   wb <- clbs[wa]
   if(here[2] < 3.5 || here[2] > 7.5) {
      if(missing(data))
         stop("data is missing")
      if(length(data) > 1)
         stop("only one class (data) allowed")
   }
   if(length(lims) == 0) {
      if(here[1] == 1)
         lims <- fit$ranges[, wa]
      else lims <- c(fit$ranges[, wa[1]], fit$ranges[, wa[2]])
   }
   if(missing(xlab))
      xlab <- as.character(wb[1])
   if(missing(ylab)) {
      if(here[1] == 2)
         ylab <- as.character(wb[2])
      if(here[2] > 6.5)
         ylab <- "probability"
      if(here[2] == 6)
         ylab <- "class"
   }
   if(missing(zlab) && here[2] == 2)
      zlab <- "probability"
   if(missing(n) && here[1] == 1)
      n <- 250
   if(missing(n) && here[1] == 2)
      n <- 50
   if(here[1] == 1) {
      cov <- matrix(cov, byrow = TRUE, nrow = n, ncol = fit$ncov)
      c1 <- lims[1] + ((lims[2] - lims[1]) * (0:(n - 1)))/(n - 1)
      cov[, wa] <- c1
   }
   if(here[1] == 2) {
      cov <- matrix(cov, byrow = TRUE, nrow = n * n, ncol = fit$ncov)
      c1 <- lims[1] + ((lims[2] - lims[1]) * (0:(n - 1)))/(n - 1)
      c11 <- (rep(c1, n))
      cov[, wa[1]] <- c11
      c2 <- lims[3] + ((lims[4] - lims[3]) * (0:(n - 1)))/(n - 1)
      c22 <- sort(rep(c2, n))
      cov[, wa[2]] <- c22
   }
   if(here[2] <= 3) {
      v1 <- ppolyclass(data, cov, fit)
      v1 <- matrix(v1, n, n)
      if(here[2] == 1)
         contour(c1, c2, v1, xlab = xlab, ylab = ylab, ...)
      if(here[2] == 2)
         persp(c1, c2, v1, xlab = xlab, ylab = ylab, zlab = zlab,
            ...)
      if(here[2] == 3)
         image(c1, c2, v1, xlab = xlab, ylab = ylab, ...)
   }
   if(here[2] == 6) {
      v1 <- cpolyclass(cov, fit)
      plot(c1, v1, type = "l", ylim = range(ww), xlab = xlab, ylab = 
         ylab, ...)
   }
   if(here[2] == 8) {
      v1 <- ppolyclass(data, cov, fit)
      plot(c1, v1, type = "l", xlab = xlab, ylab = ylab, ...)
   }
   if(here[2] == 4 || here[2] == 5) {
      v1 <- cpolyclass(cov, fit)
      v1 <- matrix(v1, n, n)
      if(here[2] == 5)
         image(c1, c2, v1, xlab = xlab, ylab = ylab, ...)
   }
   if(here[2] == 4) {
      zz <- range(v1)
      z1 <- 1 * (v1 < zz[1] + 0.5)
      contour(c1, c2, z1, xlab = xlab, ylab = ylab, levels = 0.5, 
         labex = 0, ...)
      if(zz[2] - zz[1] > 1)
         for(i in (zz[1] + 1):(zz[2] - 1)) {
            z1 <- 1 * (v1 < i + 0.5)
            contour(c1, c2, z1, labex = 0, levels = 0.5, 
              add = TRUE, ...)
         }
   }
   if(here[2] == 7) {
      v1 <- ppolyclass(cov, fit)
      plot(c1, v1[, 1], type = "l", xlab = xlab, ylab = ylab, ylim = 
         c(0, 1), ...)
      abline(h = c(0, 1))
      zz <- length(v1[1,  ])
      if(zz > 2)
         for(i in 2:zz) {
            v1[, 1] <- v1[, 1] + v1[, i]
            lines(c1, v1[, 1])
         }
   }
   invisible()
}
rpolyclass <- function(n, cov, fit)
{
      if(class(fit)!="polyclass")stop("fit is not a polyclass object")
     if(!missing(cov))cov <- unstrip(cov)
   if(n < 1)
      stop("n is wrong")
   if(is.matrix(cov) == FALSE)
      cov <- matrix(cov, nrow = 1)
   if(length(cov[1,  ]) != fit$ncov)
      stop("wrong number of covariates")
   if(n > 1 && length(cov[, 1]) == 1)
      cov <- matrix(cov, nrow = n, ncol = fit$ncov, byrow = TRUE)
   if(n != length(cov[, 1]))
      stop("cov has wrong number of rows")
   vv <- ppolyclass(cov, fit)
   ww <- runif(n)
   zz <- rep(fit$nclass, n)
   for(i in 2:fit$nclass)
      vv[, i] <- vv[, i] + vv[, (i - 1)]
   for(i in fit$nclass:1)
      zz[ww < vv[, i]] <- i
   if(min(fit$classnames) == 0)
      zz <- zz - 1
   zz
}
print.polyclass <- function(x,...)
{
   summary.polyclass(x)
}
summary.polyclass <- function(object,...)
{
      if(class(object)!="polyclass")stop("object is not a polyclass object")
    fit <- object
   it <- fit$method
   cat("========================POLYCLASS summary=======================\n")
   cat(paste("The fit was obtained with\n  "))
   cat("\b\b")
   print(fit$call)
   cat(paste("There were",fit$nclass,"classes and", fit$ncov,"covariates.\n"))
   if(it == 1) cat(paste("There were", fit$sample, "trial cases and",
      fit$tsample, "test cases.\n\n"))
   else cat(paste("There were",fit$sample,"cases.\n\n"))
   if(it == 0) {
      cat("The model selection was carried out using AIC.\n")
      if(0.99 < fit$penalty/log(fit$wgtsum)&&fit$penalty/log(fit$wgtsum) <1.01){
         cat(paste("The penalty was the default, log("))
         cat(paste(round(fit$wgtsum), "\b)="))
         cat(paste(round(log(fit$wgtsum), 2), "\b.\n"))
      }
      else{
          cat(paste("The penalty was", round(fit$penalty, 2), 
            "\b, the default would have been log("))
          cat(paste("\b",round(fit$wgtsum), "\b)="))
          cat(paste(round(log(fit$wgtsum), 2), "\b.\n"))
      }
   }
   if(it == 1) cat("The model selection was carried out using a test set.\n")
   if(it == 2) cat(paste("The model selection was carried out using", fit$cv,
      "\b-fold cross-validation.\n"))
   if((it == 1 || it == 2) && fit$select==0) {
      a2 <- range(fit$loss + diag(rep(1, fit$nclass)))
      if(a2[2] == 1 && a2[1] == 1) cat("The standard loss-matrix was used.\n")
      else cat("A loss matrix was provided.\n")
   }
   if((it == 1 || it == 2) && fit$select ==2) {
      cat("The sum of squared probabilities was used for the loss.\n")
   }
   if((it == 1 || it == 2) && fit$select ==1) {
      cat("Minus the test set log likelihood was used for the loss.\n")
   }
   if(length(fit$logl)<12)fit$logl<-t(as.matrix(fit$logl))
   a2 <- fit$logl[fit$logl[,1]==fit$ndim,  ]
   cat(paste("The model had dimension", fit$ndim, 
       "\b, log-likelihood",round(a2[3], 2)))
   if(it == 0) cat(paste(" and AIC", round( a2[2], 2), "\b.\n\n"))
   if(it == 1) cat(paste(" and loss", round( a2[2], 2), "\b.\n\n"))
   if(it == 2) cat(paste(" and AIC", round( a2[2], 2), "\b.\n\n"))
   if(it == 2){
      cat(paste("The penalty was cross-validated between", 
         round(fit$cv.aic[1],2)))
      if(fit$cv.aic[2]<0) cat(paste(" and Inf to",round(fit$cv.aic[4],2)))
      else cat(paste(" and",round(fit$cv.aic[2],2),"to", 
        round(fit$cv.aic[4],2)))
      cat(paste(" (loss",round(fit$cv.aic[3],2),"\b).\n"))
         cat(paste("The default penalty would have been log(",
         round(fit$wgtsum), "\b)=", round(log(fit$wgtsum), 2), "\b.\n"))
   }
   cat("The locations of the knots:\n")
   dimnames(fit$knots)[[2]][1] <- "Number"
   print(round(fit$knots, 3))
   cat(paste("\n  There are", fit$nbas, "basis functions, summarized below:\n"))
   a3 <- length(dimnames(fit$fcts)[[2]])
   for(i in 5:a3)
      dimnames(fit$fcts)[[2]][i] <- paste("Class",dimnames(fit$fcts)[[2]][i])
   print(round(fit$fcts, 3))
   cat("The first basis function is the constant function. For all others,\n")
   cat("the first column and the third column indicate on which covariates\n")
   cat("that basis function depends. If the third column is NA, the basis\n")
   cat("function depends on only one covariate.\n")
   cat("For the nonconstant basis functions the second and the fourth column\n")
   cat("indicate on which knot the function depend. If these columns are NA,\n")
   cat("the basis function is linear in this covariate.\n")
   cat("The remaining columns give the coefficients.\n")
   cat("\n")
   cat("================================================================\n")
   if(fit$method==0)
   cat("The influence of the penalty parameter is summarized below:\n")
   if(fit$method==1)
   cat("The effect of the penalty in the final run is summarized below:\n")
   if(fit$method==2)
   cat("The equivalence of the penalty parameter is summarized below:\n")
   dimnames(fit$logl)[[1]] <- rep("",length(fit$logl[,1]))
   fit$logl[,3:5] <- fit$logl[,3:5]/fit$wgtsum
   if(fit$method==1)fit$logl[,6:8] <- fit$logl[,6:8]/fit$twgtsum
   fit$logl[,"AIC"] <- (-fit$logl[,"AIC"])
   print(round(fit$logl, 3))
   fit$logl[,"AIC"] <- (-fit$logl[,"AIC"])
   if(fit$method==2){
  cat("The relation between the CV-loss and the penalty is summarized below:\n")
   dimnames(fit$cv.tab)[[1]] <- rep("",length(fit$cv.tab[,1]))
   print(round(fit$cv.tab, 3))}
   cat("================================================================\n")
   cat("The importance-anova decomposition is:\n")
   anova <- fit$anova
   anova[,3] <- anova[,3]*100
   dimnames(anova) <- list(rep("",length(anova[,1])),
        c("Cov-1","Cov-2","Percentage"))
   print(round(anova,2))
   cat("================================================================\n")
   invisible()
}
beta.polyclass <- function(fit, which, xsp = 0.4, cex)
{
      if(class(fit)!="polyclass")stop("fit is not a polyclass object")
 plot(c(0, 1), c(0, 1), axes = FALSE, xlab = "", xlim = c(0.1, 0.9), ylim =
                c(0.1, 0.9), ylab = "", type = "n")
        lines(c(0, 1, 1, 0, 0), c(1, 1, 0, 0, 1))
        xsp <- xsp/4
        if(missing(which))
                which <- fit$classnames
        if(fit$classnames[1] == 0)
                which <- which + 1
        if(missing(cex))
                cex <- par()$cex
        nb <- fit$nbas
        lines( c(4 * xsp,4*xsp),c(0,1))
   b11 <- fit$beta
   for(i in 1:nb) {
      b1 <- fit$fcts[i,  ]
      y1 <- 1 - (i - 0.5)/nb
      y0 <- 1 - (i - 1)/nb
      y2 <- 1 - (i - 0)/nb
      lines(c(0, 4 * xsp), c(y0, y0))
      lines(c(2 * xsp, 2 * xsp), c(y0, y2))
      aa <- fit$covnames[b1[1]]
      if(is.na(b1[2])) aa <- paste(aa, "linear")
      else aa <- paste(aa, "at", signif(fit$knots[b1[1], b1[2] + 1], 2))
      if(i==1)aa <- "constant"
      text(xsp, y1, aa, cex = cex)
      if(is.na(b1[3]) == FALSE) {
         aa <- fit$covnames[b1[3]]
         if(is.na(b1[4])) aa <- paste(aa, "linear")
         else aa <- paste(aa, "at", signif(fit$knots[b1[3], b1[4] + 1], 2))
         text(xsp * 3, y1, aa, cex = cex)
      }
      lines(c(4 * xsp + 0.03, 0.97), c(1 - (i - 0.1)/nb, 1 - (i - 0.1)/nb))
   }
   b2 <- range(b11)
   for(i in 1:nb) {
      b1 <- round(((0.92-4*xsp)*(b11[i,]-b2[1]))/(b2[2]-b2[1])+4*xsp+0.04,2)
      aa <- rep(1, length(b1))
      for(j in 1:length(b1)) aa[j] <- sum(abs(b1[1:j] - b1[j]) < 0.01)
      bb <- max(aa)
      if(bb > 1) bb <- 0.7/(nb * bb)
      for(j in which) text(b1[j], 1 - (i - 0.2)/nb + bb * (aa[j] - 1), 
         as.character(fit$classnames[j]))
   }
   invisible()
}
testhare <- c(4.974595958,0,1,2.456985,8,38,5.229125,1,3.422498434,0,0,2.177377,7,49,5.277500,0
,4.290693972,1,0,4.381446,20,54,5.485566,0,11.301950208,0,0,3.526174,10,65,4.621450,0
,10.683645663,0,1,1.150400,5,34,4.766442,1,3.741203855,1,0,5.087841,13,52,6.405083,0
,7.141522554,0,1,1.056958,8,46,4.682535,1,2.563535609,1,1,1.278860,6,25,4.556451,1
,3.701746380,0,1,3.999343,10,54,6.084539,0,6.395697579,0,1,1.336799,7,59,4.631800,0
,0.275924575,1,0,4.426891,10,45,5.141796,1,7.993160854,0,0,1.512389,8,53,4.976703,0
,10.650698724,0,1,2.227674,3,56,4.770898,0,1.015110143,1,0,5.693455,13,71,4.648958,0
,3.805403838,1,0,2.315779,5,45,4.921255,1,8.068892808,1,1,3.921555,14,29,4.820110,1
,0.944656017,1,0,6.564750,10,46,4.194352,1,1.320377850,1,1,2.505369,30,53,5.652503,0
,4.858707158,0,1,3.818449,8,24,5.283514,1,12.207398556,0,1,3.033311,3,50,4.506939,0
,10.981959783,0,0,2.896733,17,62,4.637291,0,3.407096607,0,1,1.175291,13,58,5.060192,0
,7.508765234,0,0,3.347511,9,44,5.096031,1,5.665519855,0,1,1.961776,3,30,5.067910,1
,11.655373133,0,0,3.555977,3,44,4.810457,0,1.961668982,1,1,7.169299,5,28,5.025885,1
,9.198057574,1,1,2.027242,7,40,5.077524,1,8.781112429,0,0,3.017898,19,61,4.660392,0
,1.093624486,1,1,6.143122,8,39,4.956558,1,2.924913855,1,1,4.724365,13,67,4.627196,0
,10.315301712,0,0,3.139070,15,70,6.654374,0,3.512635454,0,1,1.635826,3,50,6.034860,0
,1.883750176,1,1,2.651898,7,51,5.034317,1,4.690818787,1,1,5.931697,13,29,5.220239,0
,4.525531470,0,1,1.882305,4,48,6.011017,0,1.643812226,1,1,4.041397,17,39,5.052686,1
,4.777447362,0,1,1.906920,5,34,5.695211,1,2.127848760,0,0,2.436000,7,42,4.745345,1
,8.193705169,0,0,5.602968,15,39,5.161291,1,10.871217535,0,1,3.307412,3,50,4.506939,0
,3.646383136,0,1,3.742517,4,42,5.357143,1,1.809633167,1,1,3.949245,5,28,5.015292,1
,9.403031955,0,1,1.650917,9,66,4.374088,0,2.598629459,0,0,5.005954,20,54,5.270361,0
,0.716941135,0,0,2.367035,15,29,5.151093,1,1.149237689,1,0,3.917028,3,66,5.229125,0
,6.634921268,1,0,3.619424,2,39,4.753973,1,6.677475551,0,0,2.680901,3,59,4.619330,0
,1.853423959,1,0,7.280685,13,39,5.374839,1,1.977927943,1,1,1.479475,5,40,5.329045,1
,3.071427712,1,1,6.376802,17,61,4.761905,0,5.729473555,0,1,2.209530,4,70,5.251685,0
,12.999040825,0,1,3.899504,7,45,4.800717,1,4.059714382,0,1,3.377176,5,32,4.535342,1
,0.744850674,1,0,4.984910,20,55,6.196016,0,0.855878893,1,0,7.315716,7,47,5.851493,1
,9.801334057,0,0,1.829557,5,54,5.921052,1,8.941585984,0,0,3.592722,4,33,4.594660,1
,2.148160196,1,0,2.452756,7,41,4.230605,1,3.481338340,1,1,4.809930,7,34,5.247021,1
,6.696794779,0,1,5.234634,17,46,4.693797,0,3.835564836,1,0,4.091902,13,32,5.455447,1
,1.005947051,1,0,8.155644,17,56,5.352583,0,3.254906373,0,0,2.894001,13,38,4.860499,1
,3.673054257,0,1,3.269990,8,56,6.070261,0,6.852679664,0,0,2.728159,3,27,4.902511,1
,3.766346315,0,0,1.107997,6,52,4.599488,0,0.722085805,1,1,2.561833,13,31,4.705882,1
,5.808865352,0,1,4.302919,12,48,4.778376,0,0.758610999,1,1,4.873266,13,68,4.781478,0
,8.604434532,0,0,3.738226,10,50,4.980119,0,0.649993805,1,0,1.755723,5,62,5.077524,0
,8.846562896,0,0,2.704749,7,33,5.344762,1,8.872918703,0,1,1.049267,13,57,4.848485,0
,4.903475126,0,0,4.319373,3,52,5.294117,0,11.336698376,0,0,2.686271,4,60,4.475359,0
,3.752937025,0,1,6.675384,8,52,5.418258,1,1.177802659,1,1,1.113398,12,35,5.425139,1
,8.781171469,0,1,1.224230,8,43,5.010377,1,9.027689550,0,1,2.545398,6,23,5.728220,1
,3.495951230,0,1,4.830449,13,23,4.115462,1,3.514034237,0,1,4.004309,17,54,4.984073,0
,9.832671190,0,0,2.545287,5,22,3.697551,1,4.314590803,1,0,3.193582,10,49,4.848485,1
,1.502567605,1,0,4.507756,6,41,5.131558,1,5.588331033,0,1,4.173111,8,38,4.328138,1
,4.473941263,1,1,3.461476,17,51,5.038911,0,3.918024806,1,0,2.239662,3,43,4.666667,1
,9.091762674,0,1,1.049056,5,48,4.836185,0,1.080213129,1,1,4.159584,30,43,5.518136,1
,7.434034034,0,1,3.408190,5,43,4.984073,1,4.965729778,0,1,3.387252,8,32,5.659616,1
,4.086529910,0,1,2.015208,5,52,4.860499,0,2.178294984,1,1,4.064211,5,39,5.221878,1
,3.600221681,0,1,1.750524,10,50,4.355976,1,0.633003287,1,0,4.965685,17,22,4.923234,1
,6.274058768,0,1,1.757575,7,28,4.841229,1,0.808337851,1,1,1.089425,7,39,5.014839,1
,8.438751752,0,1,1.066097,3,53,5.454607,0,8.904313667,1,1,3.176962,13,43,5.180603,1
,5.557537136,1,1,2.930177,10,29,5.474375,1,5.199344839,0,1,1.768274,10,29,5.959141,1
,3.468550196,0,1,3.505014,5,32,4.535342,1,6.120901706,0,1,1.958902,8,46,4.682535,1
,2.638955051,1,0,5.593618,10,46,5.555451,1,5.793515954,0,1,3.807752,7,49,4.944132,1
,0.002876836,0,1,2.647479,7,22,4.960784,1,5.017167941,0,1,3.003154,7,37,4.974027,1
,3.084513249,0,1,1.570307,13,45,4.847189,1,10.665055115,0,1,3.320491,8,41,6.014000,1
,10.370884446,0,0,1.159649,5,50,5.252364,0,3.721920684,0,1,3.817728,4,42,5.357143,1
,3.608864926,1,0,2.325495,10,44,4.724556,1,4.270983923,1,1,2.324736,7,53,4.250432,0
,3.312833467,1,0,3.841955,8,58,4.827945,0,4.252580879,1,1,1.212867,8,32,5.516086,1
,1.926961555,1,1,4.700265,12,68,5.286123,0,2.337763801,1,0,3.896992,3,55,5.374839,0
,3.856019937,0,0,2.141238,8,52,4.819277,0,7.267598202,0,1,2.835711,6,33,4.508021,1
,2.815494852,1,0,3.391577,5,34,5.421687,1,6.646758005,0,1,7.461381,7,44,5.115846,1
,5.757489061,0,1,1.978295,27,67,4.908459,0,7.072637549,0,0,1.542867,7,55,4.781461,1
,0.242647827,0,1,4.614868,8,24,5.329045,1,5.721986249,0,0,3.803286,27,53,4.607373,0
,1.869202546,1,1,1.941625,9,44,5.059026,1,7.093576279,0,0,1.198225,5,43,4.660392,1
,4.924369386,1,0,7.218409,10,56,4.414404,0,5.541972634,0,1,2.106218,15,48,5.223193,1
,4.777928196,1,1,3.838427,10,34,4.902511,1,2.288150712,1,1,4.782008,11,52,5.386785,0
,3.875663734,0,1,3.145217,5,40,5.038871,1,6.435290254,0,0,1.437421,11,43,4.915615,1
,3.651424411,0,0,4.634806,3,57,5.128117,0,3.437407842,1,1,3.025336,10,53,5.583828,1
,5.942381638,0,1,1.767708,7,45,5.386785,1,1.844332249,1,1,1.122816,7,48,5.038911,1
,6.159103114,0,0,1.327351,13,64,5.153882,0,0.864839180,0,0,3.753201,8,73,4.666667,0
,4.012976033,0,0,3.832005,10,43,5.257000,1,1.308869072,1,0,2.659949,16,56,4.921529,0
,4.437839179,0,1,4.888316,10,45,4.615931,1,5.572560916,0,1,3.144377,3,44,5.951397,1
,2.299757972,1,1,4.471766,13,36,4.850811,1,1.677105018,1,1,4.428328,8,29,5.453168,1
,6.498067211,0,0,2.614110,5,53,4.861484,0,4.267064982,1,0,3.097657,3,55,5.052686,0
,1.625555118,1,1,1.600333,7,30,4.548680,1,8.675170724,0,0,2.526882,8,35,5.315730,1
,2.231627990,1,1,3.222574,15,47,5.000000,1,4.314665431,0,1,4.291517,10,59,4.652018,0
,4.281169110,0,1,4.733467,10,51,4.913402,0,3.290729834,0,1,2.427684,9,60,4.493949,0
,2.659039388,1,1,2.069452,9,47,4.276304,1,4.318971205,0,0,4.470250,8,25,4.694526,1
,5.860576271,0,0,2.539947,17,52,4.668973,0,2.315594708,1,1,4.020647,23,29,5.115846,1
,2.804584226,1,1,5.534834,7,36,5.553775,1,7.975260185,0,0,3.792803,5,62,4.660265,0
,1.106782961,1,1,9.358279,28,37,5.333333,1,2.501710933,1,0,2.566242,10,58,5.552011,0
,7.529826116,1,1,1.469639,17,51,4.672253,0,7.114662981,0,0,3.256092,10,60,5.295317,0
,4.663097138,1,1,3.591701,7,48,5.770498,0,2.589886223,1,1,3.334416,10,76,4.668973,0
,4.424446586,0,0,1.197006,8,49,4.948717,1,0.602173486,1,0,5.747700,9,37,5.517594,1
,2.288680116,0,1,6.681697,15,25,5.313040,1,4.138863511,0,0,1.686351,10,28,6.073310,1
,0.143355979,0,1,4.700276,17,67,4.686909,0,1.620006848,1,1,3.728193,13,66,4.550068,0
,4.152298128,1,1,5.478927,7,23,4.908459,1,6.251817060,0,1,4.321571,13,61,4.886249,0
,9.441347218,0,1,1.159361,27,53,4.478945,0,4.664831930,0,1,1.100395,4,33,4.594265,1
,7.672083922,0,1,3.318092,15,65,5.225269,0,2.932772797,1,0,4.330532,5,45,4.753750,1
,3.730763147,1,0,1.196599,9,75,4.292613,0,8.669274496,0,1,4.078849,10,54,4.733728,0
,3.210499902,1,0,4.826839,10,57,4.997703,0,5.740347834,0,1,3.048919,11,58,4.652018,0
,1.818485626,1,1,4.413885,22,61,6.091449,0,6.508071678,1,1,3.734937,10,55,4.892449,0
,8.548549536,0,1,1.650514,5,30,4.886216,1,4.437758745,0,0,6.218407,10,63,4.491464,0
,2.424876189,1,0,3.361363,3,55,5.374839,0,1.102067806,1,1,4.190865,7,46,6.214974,1
,4.385735762,1,1,4.993822,3,30,4.607373,1,2.014738577,1,1,6.646280,15,72,4.615620,0
,6.180637882,0,1,3.544051,10,62,5.223193,0,4.254613903,0,1,2.863435,7,50,4.503865,0
,8.403602001,0,1,3.592442,14,36,5.345836,1,4.647634726,0,1,2.318024,4,31,4.966996,1
,5.131910086,0,1,2.351699,11,40,5.277500,1,1.344906684,1,0,1.606278,8,46,5.404638,1
,1.532573575,1,0,6.670015,19,68,5.116169,0,1.232745480,0,1,5.068406,13,51,4.965363,0
,1.037380417,1,0,6.248003,8,49,4.891389,0,8.589759928,0,1,1.391212,7,62,4.145781,0
,1.790415350,1,1,2.591412,8,46,4.952207,1,8.620250368,0,0,2.692689,10,39,5.423261,1
,0.767846852,1,0,8.315848,10,43,5.219121,1,9.112477010,0,0,2.741637,8,35,5.315730,1
,0.931669122,1,1,2.236075,13,42,4.789794,1,12.337309722,0,0,1.643076,13,44,4.923659,1
,5.651648457,0,1,2.655148,7,67,5.625326,0,9.887523533,0,1,3.089551,17,51,5.336655,0
,0.527471159,1,0,5.521132,9,37,5.517594,1,2.635730106,1,1,2.267928,10,64,4.408289,0
,0.684619288,1,1,1.799920,12,30,4.272742,1,3.696458325,0,0,1.079333,6,30,5.070603,1
,1.492219376,0,1,2.218436,10,31,5.324759,1,6.978924470,0,1,1.167364,13,53,5.045599,0
,4.287556533,0,1,3.195852,20,61,4.694526,0,9.563446520,0,1,1.345802,7,53,4.789546,0
,6.543785368,1,1,3.793917,12,36,5.006571,1,2.673429072,0,0,6.431115,3,31,4.652324,1
,3.497507432,1,1,3.026543,14,37,6.250000,1,7.274422542,0,1,1.216880,5,63,4.731417,0
,11.423145428,0,1,3.599164,10,55,4.921529,0,9.941126805,0,1,3.158260,3,48,4.649801,0
,5.630463122,0,1,4.641624,10,54,4.535342,0,4.441394618,0,1,1.021535,5,53,4.521807,0
,6.712453194,0,1,1.310118,13,44,5.340002,1,0.220304507,0,0,1.260606,10,62,5.526557,0
,6.194856035,1,1,3.475943,7,39,4.075414,1,6.909127758,0,1,7.057882,7,65,4.798963,0
,6.178783402,0,1,4.908034,8,33,5.257357,0,3.377913012,0,1,6.624742,13,53,4.633481,0
,7.414119272,0,0,5.547483,13,33,4.953681,1,6.751937838,0,1,2.774286,8,39,5.521156,1
,5.701612672,0,0,7.545380,11,34,5.052686,1,2.066030486,1,1,3.576639,8,55,5.370431,0
,9.254740104,0,1,1.728768,5,22,4.881406,1,9.515433653,0,1,3.845054,17,51,5.336655,0
,4.151351126,1,0,4.953528,7,49,5.263158,1,4.188501767,1,1,4.221468,15,45,4.741448,0
,3.896630265,1,1,4.916828,8,44,5.235233,1,1.822954456,1,0,1.203984,12,49,4.864693,1
,3.568449743,0,0,5.435036,15,41,5.120809,1,1.912548820,1,1,7.668028,13,41,4.704970,1
,11.219042551,0,0,2.074424,2,61,4.635125,0,11.056621394,0,0,2.503133,17,65,6.250000,0
,7.418174278,0,1,3.182933,5,46,6.531973,1,3.800188802,0,1,6.740669,4,53,5.155131,0
,5.602112811,0,0,1.622664,13,58,4.535342,0,6.257003117,0,1,3.171081,13,42,4.430379,1
,0.518431071,1,1,5.500709,23,24,5.401257,1,4.727847199,0,1,4.545202,10,59,4.652018,0
,0.724555281,1,1,4.653086,10,49,4.991069,1,1.903682072,1,1,3.880554,10,30,4.759858,1
,2.109078189,1,1,4.012484,12,68,5.286123,0,4.259027131,1,1,3.932303,7,48,5.770498,0
,2.744858570,1,1,4.085707,10,33,4.850712,1,2.220550942,1,1,4.523721,11,52,5.386785,0
,7.380341967,1,0,1.163121,13,53,5.154913,0,1.774944681,1,1,3.547665,5,28,5.015292,1
,3.751020721,0,1,2.283102,3,33,4.419417,1,6.856077620,0,1,7.096555,7,65,4.798963,0
,1.547230496,1,0,3.234160,4,48,4.991342,0,11.157422549,0,1,5.202307,18,52,4.991342,0
,11.707553546,0,0,3.998700,7,57,4.980119,0,7.696203645,1,0,2.601323,5,54,5.247021,0
,5.096846550,1,0,1.701784,6,46,5.869379,1,0.723263842,1,1,3.526346,7,53,4.808812,0
,5.397526946,0,0,1.536233,4,56,4.631800,0,6.440626935,0,1,2.975475,6,45,4.570437,1
,1.409737724,0,0,7.226188,13,55,4.761905,0,6.641189913,0,1,4.236774,10,64,5.952871,0
,6.508748090,0,1,3.301017,12,38,6.141898,1,4.921195146,0,0,1.177132,9,33,5.142595,1
,0.664008194,1,1,3.267832,12,56,5.796012,0,2.199608074,1,0,3.221255,7,46,5.229125,1
,7.090401998,1,1,1.632593,8,50,5.270361,0,1.412069554,1,1,4.108582,13,50,4.778376,1
,7.200055097,0,0,1.507909,6,58,4.453618,0,4.685260559,0,1,8.766184,9,59,4.642308,0
,4.333250147,1,1,4.984817,13,28,5.241720,1,7.054771580,0,1,2.586440,8,31,5.624463,1
,2.085870276,1,1,1.130246,20,33,4.176713,1,2.760954702,1,1,2.942864,7,21,5.052686,0
,8.557328138,0,1,5.510820,10,52,4.466325,0,4.062899508,1,1,1.980304,23,51,5.439283,1
,4.257764630,0,1,3.032253,2,29,5.006571,1,6.023658070,0,1,3.500897,5,44,6.110101,1
,1.665516820,1,0,5.236119,10,49,5.370431,1,4.713091096,1,1,3.199575,10,49,4.784644,0
,1.136867995,1,0,6.025526,17,70,5.580490,0,1.489985563,1,1,3.291423,33,56,5.514099,0
,3.156369948,1,0,2.237694,10,44,4.724556,1,4.897130273,0,1,3.183937,5,71,4.439968,0
,2.909972143,1,1,3.821681,10,76,4.668973,0,0.473195174,1,1,4.441651,20,39,5.015292,1
,1.846312487,1,1,3.836800,8,32,5.168114,1,5.763980772,0,1,5.612299,17,33,5.014839,1
,5.833846860,1,0,1.199922,13,54,4.870861,0,1.191873836,1,1,3.212073,8,44,5.355851,1
,3.438199779,1,1,2.004780,12,67,6.013071,0,4.176913133,1,0,3.683610,17,51,5.038911,0
,3.245481586,1,0,4.344132,13,39,5.504342,1,6.260273716,0,1,2.581937,5,55,4.631800,0
,9.268476483,0,1,3.358700,7,51,4.887685,0,10.335815530,0,1,7.678052,14,52,4.701095,0
,2.965332277,1,1,5.646492,11,69,5.112992,0,2.341294041,1,1,5.263336,10,38,5.376453,1
,6.822608427,0,1,2.639990,7,28,4.808812,1,2.250602160,1,1,2.721885,5,69,5.138322,0
,6.300978500,0,0,2.470716,3,59,4.619330,0,2.465998212,0,0,6.831051,3,31,4.652324,1
,1.066376267,1,1,4.962165,17,49,4.718646,0,4.101858701,0,1,1.495405,7,57,6.959705,0
,1.693029265,0,1,3.683914,5,29,5.407597,1,1.007218322,1,0,6.338842,12,45,5.142595,1
,8.734197112,0,0,2.838642,17,68,5.090253,0,4.708322109,0,0,5.810662,12,50,5.094267,0
,5.008139408,1,1,1.529218,25,66,5.656162,0,7.197577796,0,0,3.648851,8,34,5.040121,1
,8.259847282,0,1,3.615229,4,33,5.120764,1,5.956809496,0,1,2.401528,18,57,4.960819,0
,0.135921434,1,1,4.131625,17,41,5.661268,1,4.811381111,0,1,4.171850,12,32,5.661270,1
,11.183490087,0,0,2.992922,12,63,5.554637,0,7.198701062,0,0,5.353519,17,46,4.693797,0
,6.782307836,1,1,3.679246,2,39,4.753973,1,5.698595772,0,0,7.207650,11,34,5.052686,1
,3.805061920,0,1,4.685538,16,54,4.572111,0,8.245176787,0,1,2.535451,8,53,5.549887,0
,2.845214672,0,1,3.296065,3,64,5.185781,0,9.193522879,0,0,3.012419,18,56,4.577911,0
,6.649138882,0,1,3.782074,5,31,5.303301,1,6.293812206,1,0,2.499770,6,43,4.508264,1
,1.355687434,1,0,4.915542,12,58,5.221878,0,2.837860446,1,0,1.412286,8,37,5.720019,1
,9.479260786,0,0,2.993065,5,62,4.567398,0,7.091882950,0,1,2.566071,17,56,4.860499,1
,2.213798429,1,1,3.622597,3,62,5.023578,0,4.572762277,0,0,2.816184,12,40,5.115846,1
,7.170860199,0,0,3.612320,3,55,5.796012,0,4.790314174,0,1,2.316505,7,49,4.516129,0
,1.607172815,1,0,4.952683,12,36,5.366974,1,0.896199928,1,0,5.897702,17,30,5.043083,1
,3.449437714,1,1,3.232331,28,65,4.362469,0,9.753663756,0,1,7.526090,7,53,5.416025,1
,6.601340091,0,1,1.681674,17,39,5.500175,1,6.854971010,0,0,3.997224,13,40,5.169417,1
,1.568935028,1,1,7.717937,9,33,4.575657,1,2.867387874,1,0,4.569772,13,39,5.504342,1
,2.090968802,1,0,3.980132,10,69,6.748466,0,9.157326991,0,1,3.296270,7,58,4.736275,0
,9.830822584,0,1,1.160683,5,22,4.881406,1,4.054926218,0,1,3.163021,20,61,4.694526,0
,1.440696348,1,0,3.394164,13,60,4.533199,0,8.996651002,0,1,3.415186,7,42,4.362469,1
,1.297020996,1,0,3.878364,8,54,5.120809,0,1.378073953,1,1,3.211672,13,75,5.229125,0
,2.098558667,1,1,5.720832,6,66,5.111615,0,3.759938275,0,1,4.306056,11,50,6.312191,1
,4.893984218,0,0,3.603907,7,61,5.063291,0,7.869546174,0,0,1.644335,3,42,4.322629,1
,2.344838072,1,1,4.683118,12,53,5.078968,1,4.960971610,0,1,1.089539,5,41,4.694526,1
,9.330519779,0,0,1.536900,7,53,4.536092,0,1.105600935,1,0,4.549566,7,37,5.220239,0
,0.560492445,0,0,3.006726,18,50,4.766442,1,3.049369954,1,1,1.514977,4,48,5.318160,1
,4.341731608,0,1,4.888017,7,40,5.488113,1,5.456516720,0,0,5.344640,7,50,4.529359,0
,5.750005841,0,1,1.312609,3,30,5.067910,1,1.564712708,1,0,5.955641,20,57,5.544314,0
,5.426550372,0,1,2.767213,8,49,4.603557,1,2.779461860,1,1,6.596323,12,53,5.261336,1
,7.205733380,0,0,5.545706,5,38,5.697535,1,5.048819190,1,0,1.322403,13,23,4.203487,1
,0.545143450,1,0,4.658551,20,53,4.655240,0,5.995114114,0,1,2.703438,6,45,4.570437,1
,1.765002397,0,1,6.481649,8,33,4.577839,1,3.507343644,1,1,1.281901,8,61,5.733408,0
,0.532532162,1,0,1.495249,15,55,5.248639,0,5.797894983,0,0,5.659140,13,50,4.810457,0
,6.604511283,1,1,3.800932,8,29,5.164568,1,7.346587367,0,1,4.530767,20,55,5.481173,0
,5.069804537,0,1,4.181787,12,44,5.488113,1,0.339388948,0,0,9.448597,4,42,4.698308,1
,8.358854636,0,1,2.142333,10,62,5.352583,0,6.467726356,0,1,1.634026,10,53,5.015566,0
,6.643724944,1,1,8.089670,17,41,4.733485,1,6.919238833,0,1,2.022032,7,32,4.549815,1
,1.939941117,1,1,2.050028,10,52,4.655240,0,2.982583506,1,0,5.139349,12,28,4.741448,1
,1.810415111,0,1,3.409266,5,29,5.407597,1,8.570213281,0,1,4.862987,23,36,4.930935,1
,1.855761762,1,1,3.593343,7,38,5.733508,1,0.301597154,1,0,5.599797,8,45,5.832464,1
,0.231423197,0,0,1.062920,10,62,5.526557,0,0.345455383,1,0,7.131616,3,45,5.174546,1
,4.800155273,1,1,3.416740,11,55,4.980119,0,8.094522422,0,0,3.268683,7,27,6.061189,1
,4.578820375,0,1,2.946149,3,46,5.328577,1,4.411192012,1,0,2.870890,7,43,5.075993,1
,1.998407097,1,0,1.101122,7,48,5.038911,1,1.920281072,1,0,4.796375,15,44,4.624277,1
,6.057648934,0,0,1.113554,8,64,4.374088,0,6.048239898,0,0,5.933510,5,58,4.302066,1
,7.355324373,0,0,1.706783,9,58,4.706487,0,4.214826654,0,1,1.175137,7,57,6.959705,0
,6.013621460,1,1,2.599346,17,38,4.466325,1,2.043685884,1,1,4.903934,10,61,4.631770,0
,5.407607486,0,1,4.750799,12,44,5.488113,1,4.125843671,1,1,4.593218,7,45,4.814913,1
,10.592281430,0,1,1.023707,17,64,4.533199,1,8.030761938,0,0,5.385319,13,72,4.879078,0
,0.430427654,1,1,5.721438,20,40,4.650769,1,9.111283940,1,1,2.968081,7,40,5.077524,1
,7.515296381,1,0,2.305218,4,49,5.474375,0,3.820811482,0,0,5.717769,10,50,4.893999,1
,6.483626840,0,0,1.974287,6,58,4.453618,0,4.508714258,1,1,6.328627,17,50,5.390110,1
,2.375836544,1,1,1.053986,17,31,4.885980,1,1.999463538,1,0,4.978882,13,52,4.839637,0
,3.040389022,0,1,2.748946,13,55,4.166667,0,7.211924252,0,1,1.815023,8,46,4.682535,1
,4.881413052,1,0,2.811838,5,57,4.327874,0,2.005410373,1,0,3.951752,10,53,5.685352,0
,4.609222177,1,1,5.766891,3,33,5.161291,1,6.752495877,0,1,1.183893,4,37,4.784644,1
,0.208722173,0,1,6.427892,8,30,5.055576,1,2.176248739,1,1,4.656861,22,47,5.376453,1
,3.760232022,0,0,5.711707,15,55,4.701095,0,3.925154710,1,1,1.224083,18,42,5.096089,1
,5.312882297,0,1,2.241892,4,38,5.241315,1,10.110561068,0,1,2.274060,10,56,4.706487,0
,3.453600597,1,1,3.411477,5,36,4.696845,1,5.103860480,0,1,1.733784,3,37,5.504342,1
,1.394725925,1,1,4.341496,10,27,5.015566,1,4.579546124,1,1,5.020559,18,49,5.370431,0
,1.792864341,1,0,3.038276,10,36,4.080358,1,9.071814150,0,1,2.779851,3,33,5.185781,1
,6.091612196,0,1,3.155432,12,40,4.558028,1,8.145444102,0,0,3.384657,12,42,4.960784,1
,4.180701950,0,1,2.975262,3,33,4.419417,1,11.523924806,0,1,1.243254,13,58,4.672253,0
,12.248223385,0,1,3.902886,8,45,5.164568,1,4.068336384,0,0,3.626877,6,53,5.326697,0
,3.242764985,1,0,1.518876,9,75,4.292613,0,6.572059410,0,1,1.648025,12,56,4.938272,0
,1.561962705,1,0,3.583774,4,48,4.991342,0,3.597322147,1,1,3.786477,10,33,5.299465,1
,8.625950781,1,1,1.283505,13,58,5.006571,0,1.107881892,1,1,5.897928,8,39,4.272742,1
,6.834428988,1,1,3.316387,8,29,5.164568,1,1.457088764,1,0,3.600972,15,62,4.793944,0
,9.031411884,0,1,1.342553,8,34,5.832464,1,6.058416348,0,0,3.868683,3,56,4.109609,0
,6.172218574,1,1,3.529930,13,45,4.766596,1,7.229191107,0,0,1.194385,9,58,4.706487,0
,1.655939195,1,1,2.396522,10,46,4.535342,0,4.750297681,0,0,2.868717,4,70,5.251685,0
,6.528564183,1,1,3.688580,10,44,4.800717,1,1.669746362,0,0,4.164697,17,64,5.237828,0
,3.979072929,0,0,3.743189,3,66,4.723693,0,1.083100761,1,1,5.689226,12,44,5.904718,1
,7.408603030,0,1,2.350274,8,34,5.941006,1,1.797228265,1,0,6.354184,7,46,4.827945,1
,8.524258140,0,1,1.573125,8,66,4.637013,0,2.638208892,1,1,2.877399,9,35,4.493895,1
,0.800734659,1,0,6.762334,13,65,4.921255,0,2.972883404,1,1,4.039496,5,65,4.593059,0
,3.090467846,1,1,3.326468,4,50,5.318160,0,10.907723850,0,0,7.988268,7,53,5.416025,1
,7.952081973,0,1,5.389587,13,33,4.953681,1,9.147780870,1,1,3.650056,13,43,5.180603,1
,8.456363376,0,1,3.100534,7,50,5.000000,0,1.751124609,1,1,5.671827,13,23,4.705882,1
,0.778723502,1,1,4.347169,8,51,5.153882,0,0.432523381,1,0,4.831612,3,60,5.561514,0
,2.802435199,1,1,5.764255,16,31,5.391265,1,9.200843995,0,1,1.090543,7,46,4.778846,0
,0.989988651,0,1,3.249323,10,35,4.683626,1,8.125830337,0,1,1.916480,7,52,5.201327,0
,3.900602465,0,1,2.245333,8,60,5.207717,0,6.148142134,1,1,4.151725,10,49,4.766442,1
,3.924550284,0,1,1.006315,6,62,5.669801,0,3.992218117,0,1,3.059786,10,45,5.024872,1
,9.364093883,0,0,3.155882,10,50,5.132883,0,2.041921540,1,0,6.337418,7,46,4.827945,1
,1.058448307,1,0,1.079272,12,35,5.425139,1,2.561950667,1,0,3.701132,16,60,5.242941,0
,2.891827123,1,0,4.975455,10,34,5.554567,1,1.065517223,1,0,4.755102,28,68,6.280743,0
,2.235826966,1,1,2.710385,10,30,4.561979,1,4.603598643,0,0,1.202609,5,70,4.921255,0
,3.359902494,0,0,1.549379,8,60,5.094267,0,6.284424408,0,0,4.048934,7,28,4.563989,1
,6.993184066,0,1,4.821479,13,25,5.363205,1,7.456607707,0,1,4.032171,10,54,4.733728,0
,0.323357455,1,0,5.590777,16,52,4.648111,0,10.908473437,0,0,4.361976,17,53,5.115846,0
,4.245837517,0,1,1.965004,5,32,4.750900,1,6.287846860,0,1,3.153801,10,36,5.624385,1
,0.919237043,1,0,1.130161,10,68,4.548680,0,6.168837323,0,1,1.878300,6,39,5.784654,1
,8.086577793,0,0,5.794484,12,44,5.929093,1,6.249093265,0,0,1.603284,13,37,5.359112,1
,8.853233659,0,1,4.791674,12,57,5.000000,0,0.808009138,1,0,1.653349,9,48,6.154446,1
,2.050633455,1,1,7.609329,13,41,4.704970,1,11.241442774,0,0,3.096521,15,35,5.094267,1
,7.042880279,0,1,4.289976,8,31,5.164568,1,8.608023099,0,1,4.688038,10,54,4.733728,0
,6.093956216,0,1,3.094646,5,50,4.960819,1,2.201470875,1,1,5.949804,15,39,5.687042,1
,3.846729506,0,1,3.157237,5,32,4.535342,1,10.008533299,0,1,1.044579,5,59,4.562997,0
,1.394567123,1,0,1.831151,7,55,4.624277,0,5.958579279,0,1,9.670487,10,60,5.182124,0
,7.237163870,0,0,3.682433,10,50,4.798963,1,7.344749431,0,1,3.719710,13,28,5.386379,1
,4.124655849,0,1,1.828047,14,35,5.370431,1,6.491587475,0,1,4.974500,13,61,4.886249,0
,3.461048257,1,0,3.288262,3,56,4.680553,0,6.569245920,0,1,1.003871,12,56,4.938272,0
,5.387570287,0,1,3.880464,3,62,4.364066,0,1.864536497,1,0,5.825757,13,66,4.800717,0
,5.276904936,1,1,2.449857,28,38,5.580232,1,0.180977227,0,0,6.761964,8,30,5.055576,1
,0.002821072,0,1,3.259859,6,37,4.976703,1,6.273446445,0,1,2.089837,7,28,4.808812,1
,5.248314074,1,0,4.016090,15,54,5.326697,0,4.074713586,0,1,1.614378,3,37,4.766442,1
,6.782409465,0,1,1.454855,6,33,5.474375,1,3.537634320,0,1,2.668229,7,37,6.389871,0
,4.227414512,0,0,1.983044,8,41,5.266344,1,0.306393552,1,0,4.095328,15,64,5.033223,0
,5.763045197,0,1,1.500819,6,47,6.994941,1,1.378936682,1,0,3.451202,13,51,4.771733,0
,3.336393456,0,1,1.219664,3,50,6.034860,0,2.636784901,1,1,1.505303,7,39,5.219121,1
,6.289302730,0,1,1.367362,8,50,5.270361,0,1.525525784,1,1,7.044859,13,59,4.904786,0
,1.034546670,1,1,5.103517,13,41,5.732484,1,3.084538560,1,1,3.894454,7,42,5.295317,1
,1.680179850,1,0,4.599548,10,56,4.614682,0,7.051316381,0,1,3.207134,13,38,5.178184,1
,5.178364865,0,0,2.461112,10,49,4.741448,1,10.209324904,0,1,1.720407,13,58,4.672253,0
,3.462754478,0,0,4.651890,13,30,5.474375,1,1.943983868,1,1,1.979412,9,44,5.059026,1
,0.333382673,0,1,1.829802,7,39,6.185896,1,7.365314347,0,1,3.992740,3,35,4.835737,1
,6.456711297,0,0,1.955456,10,48,5.391265,1,5.352734390,0,1,1.381658,12,60,4.603557,0
,3.117490032,1,1,4.862531,8,39,4.827945,1,10.464199781,0,1,2.654715,3,56,4.770898,0
,6.648622390,0,0,1.418481,6,58,4.453618,0,0.504983812,1,1,2.871429,10,33,4.681194,1
,1.132197103,1,1,5.547047,20,68,4.619330,0,8.204686747,0,0,2.118565,7,33,5.344762,1
,0.548629641,1,0,8.382483,12,47,6.641995,1,6.545835318,0,1,2.851951,16,61,4.830680,0
,1.614850887,1,0,5.732372,20,41,4.907975,1,4.326508009,1,1,1.128477,10,51,5.877699,1
,9.631693548,0,1,1.589332,7,53,4.987757,0,3.743417556,0,0,4.826593,10,41,5.021689,1
,0.711919688,1,0,5.311243,8,34,6.344507,1,1.406760777,1,0,7.458299,8,67,4.723693,0
,5.347828937,0,1,1.984332,12,60,4.603557,0,4.511347593,0,1,4.925597,3,45,5.590170,1
,5.807706633,0,1,3.428959,8,34,5.422386,1,4.077919875,0,1,6.562409,30,55,4.556451,0
,6.051981912,0,1,1.919870,8,58,5.257000,0,5.732195844,0,1,4.652636,8,65,4.467590,0
,6.597855218,0,1,1.732938,13,44,5.340002,1,10.118567072,0,1,3.887820,4,45,6.558120,1
,6.110741125,0,0,2.134272,5,53,4.861484,0,0.231828608,0,0,4.035851,17,60,5.882353,0
,0.777022351,1,0,3.151760,7,30,5.295317,1,2.341648323,0,0,2.016336,4,40,5.103104,1
,6.279425508,0,0,1.674932,17,39,5.555556,1,10.721694673,0,1,2.376518,6,53,4.899540,0
,9.511700864,0,1,4.601852,15,51,5.106757,0,6.266444418,0,0,2.530300,18,57,4.960819,0
,7.547995835,0,1,4.582309,5,32,5.625000,1,1.458256258,1,1,4.566258,25,31,6.240738,1
,1.733922908,1,1,4.160118,7,26,5.763505,1,10.296185095,0,0,6.443537,5,62,5.588507,1
,0.616789585,1,0,6.848180,33,51,4.901409,0,0.190548415,0,1,1.433576,3,36,5.096031,1
,4.597437441,0,1,4.877607,7,37,5.624713,0,0.487907532,0,0,2.016946,8,37,5.237828,1
,4.573749060,1,1,4.648474,13,28,5.241720,0,6.351624030,1,1,1.223701,7,48,5.416645,0
,5.444808172,0,1,5.103662,18,37,4.252083,1,1.510551130,1,1,1.072253,2,47,4.330127,1
,1.584217669,1,1,7.082270,13,59,4.904786,0,7.480105808,0,1,4.115575,17,45,4.933737,1
,0.264522374,0,1,3.556272,17,31,6.419274,1,2.422858835,1,1,3.409903,5,32,4.760953,1
,5.450188944,1,0,1.790452,6,46,5.869379,1,6.809018450,0,1,2.416089,6,43,4.508264,1
,2.337273737,1,1,4.383944,8,48,4.841229,1,2.079221042,1,1,5.173593,15,68,5.039189,0
,8.267420054,0,0,1.442905,3,42,4.322629,0,6.906518314,0,1,4.727592,8,31,5.164568,1
,6.078153858,1,1,4.333743,6,51,5.767761,0,7.729096555,0,1,3.453160,13,54,5.266344,0
,6.317741962,0,0,1.894594,17,36,4.713139,1,3.696723353,0,0,2.437806,5,59,5.142595,0
,9.527512819,0,1,3.724191,8,55,5.015292,0,5.980060589,0,1,3.536102,5,50,4.960819,1
,5.877707151,0,0,2.010150,10,49,4.741448,1,7.932894615,0,1,3.430362,4,33,5.120764,1
,2.042875027,1,1,4.943407,12,68,5.286123,0,6.313682701,0,1,3.987807,27,38,5.095541,1
,3.943776047,0,0,4.250438,11,48,4.869480,0,4.190370232,0,1,3.602082,10,42,6.369427,0
,0.713752810,1,0,4.484808,5,59,5.642155,0,2.061156533,1,0,6.993872,10,58,4.997703,0
,7.265972062,0,1,2.504398,8,34,5.941006,1,5.897122381,0,1,5.522065,14,54,4.563989,0
,1.117037792,1,1,2.476285,15,37,5.296764,1,2.519232332,1,1,3.615621,13,40,4.956558,1
,6.722139594,1,1,2.970835,5,46,5.735394,1,2.939263624,1,0,4.783918,9,44,5.553775,1
,1.263742230,1,0,5.732667,2,41,5.327739,1,3.554809385,0,0,2.372134,8,52,4.819277,0
,2.005697320,1,0,3.737095,10,69,6.748466,0,9.113218144,0,0,1.729485,3,54,4.810457,0
,1.903591075,1,0,3.448379,8,66,4.916011,1,6.852265320,0,1,3.989755,10,50,5.811836,0
,1.148551185,1,0,7.176261,30,50,5.486694,1,11.838745909,0,0,2.502532,4,60,4.475359,0
,9.465308131,0,1,3.676296,8,55,5.015292,0,3.652735761,1,0,3.855869,4,59,5.521156,0
,7.709119867,0,1,5.872373,8,37,4.879415,1,11.227722208,0,1,1.768116,3,35,5.304117,1
,4.514278936,0,1,3.831382,3,43,5.050762,1,8.168211337,0,1,3.540189,14,36,5.345836,1
,2.110967002,1,1,4.250007,10,61,4.631770,0,0.170603433,0,1,2.554719,11,55,5.151093,0
,6.142549669,0,0,4.511429,12,37,4.976703,1,7.829874802,0,0,1.761346,10,27,5.290843,1
,4.041970755,0,0,3.166400,7,56,4.615620,0,0.881428488,1,0,4.963674,15,38,7.096774,1
,3.895394592,0,1,2.966590,7,50,4.503865,0,4.734134560,0,1,1.159845,8,42,4.921529,1
,5.549556897,0,1,2.515525,12,42,6.321264,1,5.957744865,1,1,4.282083,10,29,5.096089,1
,2.612217756,1,1,3.283575,5,32,4.760953,1,2.351318893,1,0,9.433906,7,38,5.474375,1
,3.929335919,0,1,1.969351,3,52,4.938272,0,6.057858733,0,1,3.963494,5,50,4.960819,1
,11.505287025,0,1,2.650312,9,45,4.705882,1,2.425573164,1,1,3.449912,12,40,4.870861,1
,8.537079653,0,0,1.919414,10,66,4.686909,0,10.535590867,0,0,5.666118,7,69,4.506939,0
,2.136745227,1,1,6.834885,15,25,5.313040,1,1.926380538,1,1,3.226138,7,38,5.733508,1
,8.115219096,0,0,1.465175,17,54,4.364066,0,5.472815064,0,1,1.782311,8,37,5.856070,1
,1.130849975,1,0,4.216702,28,68,6.280743,0,1.163170736,1,0,7.256262,30,50,5.486694,1
,3.365204161,0,1,4.654139,33,56,4.904786,0,1.162584948,1,1,3.640609,10,36,5.333006,1
,5.522293804,0,1,3.506413,10,34,4.881406,1,4.280830294,0,1,2.784282,6,56,4.540842,1
,0.536465463,1,1,6.009152,8,23,5.519851,1,5.372209310,0,1,4.339617,3,23,4.655240,1
,8.092996790,0,0,4.227226,10,44,5.476925,1,1.660042695,1,1,4.237032,10,43,4.773922,1
,5.046438225,0,1,1.891330,27,67,4.908459,0,0.558458981,1,1,2.598031,4,32,5.487283,1
,2.431132009,1,0,4.795160,12,42,4.423004,1,3.911608432,0,1,3.265681,7,56,4.615620,0
,1.728596360,1,0,5.031483,27,70,4.604683,0,6.485542586,0,1,2.892669,18,57,4.960819,0
,4.455368793,1,0,1.558971,25,66,5.656162,0,3.356534447,0,0,3.694339,15,69,4.902490,0
,3.764739669,1,1,3.220591,12,32,5.806452,1,8.603190010,1,0,3.052693,21,67,4.610694,0
,11.808447247,0,1,2.969251,17,54,4.851086,0,0.937717720,0,1,3.799739,10,35,4.683626,1
,0.356154239,0,0,1.705729,10,41,4.631800,1,7.749694261,0,0,1.176630,8,60,4.778846,0
,0.377782457,1,1,3.880753,3,33,4.944419,1,1.442296832,1,0,3.641154,13,60,4.533199,0
,2.966253607,1,1,1.941936,7,38,4.672253,1,4.243549545,0,0,4.388806,17,51,4.668973,0
,5.085674026,0,1,4.854721,8,65,4.467590,0,6.941571441,0,0,3.838433,12,45,4.422167,1
,1.894093992,1,1,3.204225,3,45,5.520686,1,1.899673674,0,0,4.416458,3,53,4.835737,0
,3.350516263,1,1,3.591717,4,59,5.521156,0,0.412278028,1,0,5.988977,18,51,5.081007,0
,2.308422808,1,1,4.487907,17,42,4.800717,1,4.680489773,1,0,2.724866,3,64,5.201457,0
,2.165319102,1,0,4.384848,13,52,4.839637,0,11.550136209,0,1,1.629800,7,50,4.860499,1
,7.761382067,0,0,1.970176,8,66,4.686909,0,2.709093688,1,1,3.535508,13,63,5.252364,0
,2.795778440,1,1,4.713039,10,33,4.850712,1,3.655220622,0,1,2.232516,7,50,4.503865,0
,5.131761460,0,0,3.191143,10,41,5.096031,1,1.933563557,1,1,1.074989,8,44,5.247021,1
,2.386244466,1,0,3.105901,12,57,5.652957,0,6.579558296,1,1,4.911388,10,52,4.686909,0
,1.653981343,1,1,3.754309,7,49,4.693797,0,3.616694513,0,0,3.340925,3,52,4.590991,0
,3.348087451,1,1,3.596255,7,31,6.202187,1,3.937152011,0,1,7.048824,11,26,4.871677,1
,4.448652618,0,0,1.899213,5,25,5.007613,1,3.314438255,1,0,3.568579,10,47,4.921255,0
,4.378647088,0,1,3.082401,8,32,5.659616,1,2.894753047,1,1,4.472487,10,27,5.454546,1
,1.143836483,1,0,4.906165,20,59,4.583412,0,2.078924038,1,1,5.038161,8,41,4.759858,1
,7.632391183,1,1,3.766607,8,29,5.164568,1,2.532663288,1,0,2.514690,27,65,5.000000,0
,7.292980488,0,1,3.598738,13,38,5.178184,1,7.161655578,0,1,3.531391,13,38,5.178184,1
,7.995158576,0,1,1.195915,27,61,4.631800,0,5.356725785,0,1,2.807563,3,66,4.983549,0
,5.864581357,0,1,1.653911,7,59,4.631800,0,3.815578056,0,1,3.490062,4,42,5.357143,1
,2.088531185,1,1,4.074264,15,58,5.646925,0,5.654649446,0,1,2.197629,4,70,5.251685,0
,8.609041606,1,0,4.729127,12,50,5.163978,0,2.017008686,1,0,1.169704,10,51,5.929271,0
,1.592183384,1,1,2.791577,5,53,5.521473,0,11.340743492,0,1,3.460046,8,45,5.164568,1
,2.497296371,1,0,3.562934,10,69,4.637013,0,0.924902449,1,0,5.116869,17,40,5.094267,1
,0.877101983,1,0,6.776094,13,65,4.921255,0,2.914522704,1,1,3.430647,7,48,5.266344,1
,7.691193425,0,0,1.529326,6,51,5.516086,0,3.500237395,1,1,1.383206,5,43,5.420771,1
,6.667544048,0,1,3.596455,8,42,5.120432,1,5.717699502,0,1,3.406005,5,44,5.070667,1
,1.343640841,1,1,6.402983,5,46,4.819277,1,5.281411490,0,0,3.603305,9,52,4.923659,0
,8.514367758,0,0,2.699022,5,62,4.567398,0,0.345293636,1,0,7.424966,3,45,5.174546,1
,3.256129004,1,1,4.639832,13,66,3.875617,0,3.943650016,0,0,1.172036,10,28,6.073310,1
,5.750799191,0,1,3.048840,7,38,5.242941,1,0.837870163,1,1,4.224013,13,50,6.273158,0
,1.888124836,1,1,4.047443,8,29,5.453168,1,5.908433161,0,0,1.920029,20,41,4.454354,1
,1.049472063,1,0,5.003171,12,44,5.904718,1,3.586541725,0,0,3.288140,7,34,5.095541,1
,6.635303664,0,0,5.943127,14,54,4.563989,0,5.947873756,1,1,1.189624,6,39,5.784654,1
,7.800592523,0,1,5.629727,14,52,5.514311,0,1.939952102,1,1,3.039541,15,43,6.136303,1
,1.348567277,1,0,3.943967,15,61,5.455447,0,3.894577310,1,1,2.075025,7,48,4.843221,0
,3.193398435,1,1,3.658782,28,65,4.362469,0,4.817679250,0,0,1.025435,9,33,5.142595,0
,2.860194480,1,1,1.165846,23,45,4.733485,1,1.352776048,1,1,1.676759,27,30,5.454546,1
,4.912300494,1,1,3.153114,7,32,5.381357,1,4.890521981,0,1,2.582964,13,57,4.463000,0
,5.764706540,0,1,2.609888,5,30,4.907975,1,2.767113746,1,0,2.928604,7,57,4.798963,0
,4.431233163,0,0,1.345869,10,41,5.454546,1,1.623662089,1,1,3.463678,27,45,4.827945,1
,2.181049038,1,1,3.567826,3,62,5.023578,0,5.721996131,0,1,3.195724,6,42,5.452375,1
,1.627023852,1,1,5.494793,12,51,5.421687,1,2.751551764,1,1,1.687466,5,41,6.545970,1
,7.271223955,0,1,1.774302,10,25,4.980119,1,5.997232191,1,0,4.188447,7,57,5.617264,0
,2.815038151,1,1,4.098625,8,39,4.827945,1,6.538407808,0,1,3.898660,20,39,4.493895,1
,3.586533365,1,1,2.288991,13,31,5.397807,1,1.168339111,1,1,3.037992,4,59,5.014839,0
,1.131764089,1,0,6.213348,5,42,4.535342,1,5.773740162,1,1,5.733903,10,32,4.886216,0
,5.724049911,0,1,3.742550,6,46,5.034317,1,8.209486959,0,1,1.318034,7,43,4.960819,0
,5.905065192,0,0,1.436360,4,39,5.034317,0,1.312351945,0,1,2.809238,7,44,6.481796,1
,5.327245863,0,0,1.840784,20,41,4.454354,1,2.016230894,1,0,2.050312,10,43,5.386785,1
,7.273675921,0,1,3.319160,8,33,5.303301,1,0.679058643,1,1,4.789431,10,49,4.991069,1
,6.208641571,0,1,1.570053,7,28,4.841229,0,10.193101052,0,1,2.174288,10,49,5.081007,0
,7.983991138,0,1,1.496593,7,23,5.228350,1,2.707835686,0,0,1.239440,12,65,5.219121,0
,1.505188818,1,1,7.607318,13,59,4.904786,0,0.849979285,1,0,5.965647,8,50,5.625326,0
,2.102019053,1,1,4.984487,13,39,4.567398,1,3.535572626,0,1,3.852810,13,50,5.059026,0
,2.768660076,0,1,3.661696,18,44,5.257000,1,10.673161184,1,0,3.207696,8,42,5.589223,1
,4.175695912,0,1,1.367948,7,53,4.960819,0,1.419040855,0,0,1.312853,10,48,4.311626,0
,1.093027247,1,1,3.994409,5,53,5.359112,1,1.966217633,1,0,5.150978,9,51,4.593059,0
,2.128991491,1,0,3.787845,10,53,5.685352,0,8.002784875,0,1,1.808558,5,27,4.950651,1
,6.583644943,0,0,1.654456,17,39,5.555556,1,1.113495390,1,1,4.798297,8,49,4.652018,1
,0.550321056,1,1,6.502562,8,23,5.519851,1,0.723488188,1,1,2.155395,13,31,4.705882,1
,0.954077783,1,0,3.535933,3,66,4.525292,0,4.975119414,1,0,1.484172,7,43,4.733485,1
,6.969411931,0,0,1.264260,6,58,4.453618,0,6.712535366,0,0,2.394560,12,75,4.742505,0
,2.433923883,1,0,2.138953,4,59,5.576548,0,10.577546171,0,0,2.936957,17,65,6.250000,0
,4.503652270,0,1,1.208192,23,54,4.712121,1,2.298074342,0,0,1.138595,7,43,5.104738,1
,3.400328945,0,1,1.293391,13,32,4.694526,1,1.984224372,1,0,2.070662,7,42,4.745345,1
,5.792484433,0,0,1.958561,10,46,5.940885,1,5.048172060,0,1,4.938344,10,53,5.009940,0
,3.724314833,0,1,4.003840,10,60,5.359738,0,0.703710942,1,1,1.118981,7,39,5.014839,1
,3.068358640,1,0,4.306813,12,59,6.130060,0,2.793717710,1,1,1.469506,10,41,4.800717,1
,0.924667925,1,0,4.032832,7,37,5.220239,1,6.724213596,0,0,7.648622,33,59,6.595520,0
,5.042581582,1,1,3.873687,17,45,5.423261,1,6.244645554,0,1,1.163991,13,40,5.484352,1
,1.069207876,1,0,5.898909,12,44,5.904718,0,9.685962357,0,1,7.351507,7,53,5.416025,1
,0.602467925,1,1,4.164720,17,61,4.385608,0,6.735501806,0,1,3.681254,12,40,4.558028,1
,2.416085327,1,1,2.998528,13,41,4.991342,1,0.883339159,1,0,4.666377,22,25,4.960784,1
,1.029780683,1,1,7.469264,30,50,5.486694,1,2.678863910,1,1,6.210031,12,53,5.261336,1
,0.218703943,0,0,5.567017,5,52,5.920780,0,1.135487874,1,1,4.155227,20,49,4.923659,0
,9.563382952,0,1,3.367165,6,35,5.514311,1,9.073131862,0,0,2.514044,8,47,5.217020,1
,4.171939803,0,0,3.681111,7,56,4.615620,0,1.519229517,1,1,3.118656,29,67,4.563989,0
,10.778655855,0,1,4.842721,5,45,5.329681,0,4.337965989,0,1,3.301727,10,57,4.985775,0
,7.280470398,0,1,4.031752,17,51,4.841229,0,1.007990146,1,0,5.502734,15,45,4.988877,1
,2.495596502,1,1,3.492195,11,44,5.163978,1,7.608110567,0,1,1.963279,12,65,4.590991,0
,2.118993806,1,1,3.483719,10,68,4.472136,0,3.558037083,1,1,1.582273,5,55,4.519892,0
,6.017214503,0,1,7.847679,17,60,5.052686,0,11.540110101,0,0,1.896460,12,52,5.416025,0
,1.668611958,1,1,3.480821,10,63,4.499433,0,3.584244189,1,1,4.929268,7,34,5.247021,1
,2.127077850,1,1,2.093565,7,64,5.488114,0,9.608429844,0,1,3.034624,8,55,5.015292,0
,3.864569215,0,0,5.000620,10,50,4.893999,1,12.183432204,0,1,1.583632,12,44,4.830909,1
,1.336573082,1,1,3.780255,8,54,5.120809,0,2.681474118,1,1,2.129729,14,36,5.863020,0
,3.537433598,1,1,3.872712,12,32,5.806452,1,0.588255172,1,0,8.332843,12,47,6.641995,1
,0.454187936,0,1,3.252996,3,39,4.758241,1,6.768400653,0,0,1.001791,11,43,4.915615,1
,1.388905490,1,0,3.142672,17,63,4.328138,0,7.651766392,1,1,6.356338,2,33,5.962848,1
,3.187847518,0,1,3.994416,8,27,5.583828,1,3.393697460,0,1,3.739596,3,53,5.115846,0
,9.269941007,0,1,3.593269,8,50,4.632703,0,1.861529935,1,0,4.709698,22,47,5.376453,1
,8.697841282,0,0,4.726717,8,55,4.385608,0,5.412574003,0,1,3.330418,10,57,5.153882,0
,6.995501779,0,0,2.478400,3,59,4.619330,0,4.330310494,0,0,3.917717,3,66,4.723693,0
,2.576165090,1,1,6.066838,18,54,5.661270,1,1.145485281,1,0,5.285216,2,41,5.327739,1
,5.329422753,1,1,3.671430,7,62,4.766442,0,3.614254190,1,0,2.784692,7,54,5.677647,0
,4.165876197,1,1,6.332982,5,41,4.944419,1,4.256215571,1,0,3.280260,8,58,4.550068,1
,5.088455679,0,1,4.052454,13,39,6.069946,1,6.116862102,0,1,3.977267,12,38,6.141898,1
,12.683461723,0,0,1.517336,5,70,5.043558,0,4.453783981,1,0,3.871339,10,49,4.848485,1
,8.785755092,1,1,1.782309,12,56,4.752127,0,4.130454944,1,1,6.445827,17,50,5.390110,1
,0.404675946,1,0,3.814497,11,53,5.064476,0,4.072113862,0,0,5.897269,9,29,4.561979,1
,0.740543355,1,0,6.359513,8,51,4.902490,0,3.758998262,1,0,4.317449,17,40,5.201327,1
,6.319724316,1,1,5.784411,5,54,4.899540,1,3.303044922,1,1,7.979470,13,48,5.404638,1
,7.902598376,0,1,1.160953,10,59,4.682826,1,1.647431277,1,0,3.702247,13,57,5.096031,0
,5.017489872,0,1,3.201001,3,42,4.781478,1,9.130907734,0,1,4.749968,7,42,5.421048,1
,2.178538364,1,0,2.171074,17,59,4.796997,0,11.721630968,0,0,1.572511,12,52,5.416025,0
,4.982466788,0,1,1.210030,7,51,4.503865,0,0.989578057,0,0,7.625260,20,62,4.660392,0
,0.198328257,0,1,3.333207,5,52,5.318160,0,1.195690138,1,0,2.651865,17,52,5.229125,0
,1.917606493,1,1,2.632788,9,48,4.704970,0,3.186227065,1,1,1.684827,7,38,4.672253,1
,3.760231427,1,1,1.262469,10,51,5.877699,1,0.371193691,1,0,6.507787,17,52,6.180629,0
,6.708271608,0,1,2.340167,8,34,5.941006,1,3.092500668,1,0,4.816183,22,64,4.801516,0
,8.129999085,0,0,3.453010,15,49,4.520859,0,5.612371164,1,0,1.525717,10,44,5.160907,1
,11.718612459,0,1,3.267709,17,46,4.704970,0,1.803107106,1,0,1.905088,27,67,4.742505,0
,7.042425661,0,0,4.871394,10,54,5.229125,0,4.486400924,0,1,2.693333,7,60,4.991342,0
,2.458152764,1,1,4.166088,8,45,4.724556,1,5.081226855,1,1,3.079174,7,56,4.991069,0
,4.851088949,0,0,1.655767,5,70,4.921255,0,7.004463844,0,0,1.899629,8,32,4.412188,1
,8.060202675,0,1,2.705815,10,34,5.033223,1,1.011714856,1,1,3.617589,21,45,5.993707,1
,0.694334841,1,0,2.869502,12,52,4.933737,0,4.379653798,0,1,1.582491,5,33,5.115846,1
,1.992834936,1,0,3.526825,3,25,5.451704,1,1.692036651,1,0,4.520205,17,44,5.333006,1
,4.328743730,0,1,1.668471,7,47,4.594265,0,5.694383361,0,1,1.481445,3,30,5.067910,1
,1.597763919,1,1,5.364744,20,41,4.907975,1,3.566550944,1,1,2.005828,7,48,4.843221,0
,5.294380126,0,0,5.128683,7,50,4.529359,0,1.312443267,1,1,4.076913,15,72,4.896896,0
,2.146343958,1,1,3.147547,12,41,4.983549,1,0.499759341,1,0,2.347291,18,36,4.506939,1
,0.377100646,0,1,4.048678,7,46,5.412659,1,5.598202728,1,1,2.681853,20,35,5.842951,1
,11.677910276,0,0,4.933050,22,66,4.810457,0,1.700213530,1,0,4.104284,12,53,4.516129,0
,9.020224218,0,1,1.600587,5,43,4.902490,1,6.147254227,0,1,2.120357,8,49,4.603557,1
,7.100147957,0,1,1.412633,17,39,5.500175,1,3.997298399,1,1,5.518855,6,43,6.373774,1
,3.041361978,1,0,4.398945,8,46,5.128117,1,8.169657531,0,1,1.744003,7,23,5.228350,1
,6.510418634,0,1,1.718795,10,53,5.015566,0,6.828632147,1,1,2.016646,12,48,5.709323,0
,0.861517178,1,1,5.276851,8,40,4.902511,1,4.886619000,1,0,6.515116,17,45,4.648958,1
,9.676396170,0,1,3.117316,11,52,5.201327,0,3.716733232,0,1,2.108093,13,67,4.444445,0
,1.115358443,1,1,5.601878,20,68,4.619330,0,1.303703540,1,1,3.507072,12,42,5.685352,1
,0.196904947,0,1,3.018814,10,26,4.643764,1,1.976208169,1,1,3.824025,27,42,5.509923,1
,0.531891014,1,1,3.220245,20,48,6.090869,1,3.118361672,1,0,5.715606,17,45,4.413292,1
,9.422584610,1,0,2.823839,12,24,4.830680,1,0.002647261,0,1,2.844084,7,22,4.960784,1
,1.108267800,1,0,6.904618,8,49,4.891389,1,2.908533443,1,1,2.273455,14,36,5.863020,1
,6.778521025,1,1,2.785345,8,37,5.554637,1,6.173335308,0,1,3.364887,8,60,4.913402,0
,6.935856038,0,0,3.818273,6,62,5.615465,0,0.166028795,0,0,5.748957,15,38,6.033400,1
,2.725026051,1,0,4.288533,8,46,5.128117,1,4.612444617,1,1,3.030475,8,40,5.657501,1
,1.281861622,1,1,5.674807,2,41,5.327739,1,9.339026857,1,1,3.871680,13,43,5.180603,1
,10.735606323,0,1,3.281375,8,55,5.015292,1,2.365873897,1,0,2.535791,27,65,5.000000,0
,11.766978135,0,0,2.478327,18,57,5.420771,0,4.271065215,0,1,1.213050,11,41,5.412659,1
,0.616397869,1,1,4.264361,17,22,4.923234,1,4.027677286,0,0,1.190329,10,28,6.073310,0
,8.809499703,1,1,3.195532,7,56,4.242424,0,0.945373395,1,1,1.512612,5,55,4.933303,0
,4.457597521,0,1,4.558692,17,51,4.668973,0,10.968263637,0,0,1.326111,5,70,5.043558,0
,1.893454171,1,1,5.306797,7,32,4.839637,1,3.071826066,1,0,1.022251,10,34,5.266344,1
,3.332564283,0,1,4.634036,7,34,5.247021,1,5.498649241,0,0,1.586036,6,46,5.096031,1
,0.983960569,1,1,3.265702,8,46,4.830680,1,0.002802357,0,1,2.041719,7,22,4.960784,1
,0.146355799,0,1,3.405241,5,35,5.520686,1,10.276878528,0,1,1.884935,5,22,4.881406,1
,5.476632789,0,0,4.334094,17,69,5.071590,0,2.890838924,0,1,5.516277,16,31,5.391265,1
,0.740958745,1,0,8.021389,10,43,5.219121,1,5.918680277,0,0,2.598374,13,46,5.080005,1
,6.252165398,0,1,3.278764,7,63,4.374999,0,8.144544727,0,1,5.920142,14,52,5.514311,0
,7.272492165,0,1,3.646474,10,50,5.811836,0,6.113740440,0,1,1.835806,7,44,4.984073,1
,0.240082926,0,1,3.061189,17,31,6.419274,1,2.826153656,0,1,3.521198,18,44,5.257000,1
,6.600704658,0,1,1.072442,4,37,4.784644,1,0.589988940,1,0,1.087028,15,55,5.248639,0
,3.565092881,0,1,4.115564,11,50,6.312191,1,2.239053344,1,1,3.994855,12,40,4.870861,1
,3.933021427,0,1,1.398644,4,48,6.011017,0,6.621891045,0,1,2.730461,12,26,4.247670,1
,10.185764983,0,0,6.359992,5,62,5.588507,0,0.895861781,1,1,5.710522,8,40,4.902511,1
,5.411703008,0,1,1.056656,5,43,4.864693,1,8.565351894,0,1,3.348079,8,54,4.451705,0
,6.421887837,0,1,9.766629,10,60,5.182124,0,6.298610506,0,0,1.253611,8,47,4.191617,1
,0.296036847,0,1,4.810572,17,55,4.791564,0,10.851082559,0,1,3.814584,12,48,4.861484,0
,7.980208024,0,1,5.487632,10,52,4.466325,0,3.638660988,0,0,2.472651,8,40,5.326697,1
,7.442259017,0,1,3.339889,3,61,4.619330,0,7.968494005,0,1,2.078416,8,31,5.624463,1
,5.638054393,1,1,2.112804,10,29,5.474375,1,2.274404360,1,1,3.887508,11,44,5.163978,1
,2.806913498,1,1,3.938302,20,64,4.650930,0,4.014024757,0,1,4.271452,10,51,4.913402,0
,4.495647291,0,1,1.342853,17,29,4.752127,1,8.205679002,0,1,1.964848,7,52,5.201327,0
,10.934958368,0,1,3.501997,28,68,6.308775,0,1.075246807,1,0,3.377548,17,58,6.343058,0
,10.341693264,0,0,1.863250,12,44,4.830909,1,7.693085207,0,0,3.165204,7,54,5.481173,0
,5.395744543,0,1,4.457089,12,48,4.778376,0,5.382884698,0,1,3.773782,6,42,5.452375,1
,6.508818756,0,1,2.879026,7,35,5.219121,1,2.990790446,1,1,2.670362,10,61,4.677072,0
,2.069391739,1,1,3.632507,3,67,4.224999,0,1.827157268,0,1,3.541418,5,29,5.407597,1
,1.100326438,1,0,3.830103,33,30,5.132883,1,1.041723950,1,0,5.027721,10,42,5.174546,1
,4.441376641,0,1,2.877808,8,41,5.219121,1,3.330571019,1,1,5.685433,10,39,4.650930,1
,2.516389054,1,0,2.964775,4,54,4.549815,0,6.898805381,0,1,1.874730,4,37,4.784644,1
,11.377748124,0,1,4.265396,3,58,5.652957,0,5.471347402,0,1,5.358442,13,50,4.810457,0
,9.166552745,0,1,3.909758,7,64,4.820110,0,5.878316974,0,0,5.537839,4,24,4.572111,1
,0.786873728,0,1,4.422472,8,69,4.519892,0,4.164596816,0,1,1.600253,3,43,5.237828,1
,11.542036207,0,0,4.889596,6,53,4.952207,0,6.423826599,1,1,3.495273,10,67,4.419417,0
,3.506319462,0,0,1.130150,8,60,5.094267,0,9.498279253,1,1,3.506834,5,62,4.839637,0
,12.609272471,0,1,5.727082,13,36,5.344762,1,4.777893721,0,0,4.754430,4,52,4.829433,0
,2.623103825,1,1,4.749771,10,33,4.850712,1,2.883005385,1,0,2.714045,7,57,4.798963,0
,5.539691312,0,1,2.584135,3,43,4.984073,1,9.465759539,0,1,2.798065,10,62,5.352583,0
,4.067581461,1,1,3.001352,9,32,5.290592,1,3.096699745,0,1,3.091198,8,27,5.583828,1
,4.130818534,0,1,3.075352,10,42,6.369427,1,2.200995115,1,1,2.587967,8,43,5.194805,1
,0.844702499,0,0,3.591533,8,73,4.666667,0,2.302195126,1,0,5.316506,7,41,4.784644,1
,7.954934108,0,1,5.940959,17,25,5.132883,1,1.886116862,1,1,3.103121,12,45,4.983549,1
,4.210835201,0,0,1.733284,10,54,5.717564,0,2.131763536,0,1,4.552194,12,31,5.333333,1
,2.606818892,1,1,4.143072,13,53,4.983549,0,0.967519928,0,1,3.333475,20,30,4.631800,1
,1.577373909,1,0,6.969599,11,44,4.744147,1,0.928473702,1,0,5.932476,13,71,4.648958,0
,3.928922370,1,0,3.830729,7,35,4.771733,1,8.851494640,1,0,1.376508,17,59,5.419018,0
,11.587160723,0,1,2.125626,13,44,4.606335,1,9.169098320,0,1,3.465583,14,44,5.085716,1
,5.396104741,0,1,3.956372,8,28,5.770498,1,1.627566038,1,0,4.930407,20,69,4.319955,1
,2.878920068,1,0,2.669477,4,59,5.576548,0,8.402053876,0,0,1.943587,7,57,4.250432,0
,6.477433818,0,1,3.327807,15,60,5.254470,0,2.249856597,1,1,1.139254,8,49,5.303301,1
,6.756788662,0,1,3.181352,10,46,5.447472,1,4.565235126,0,1,3.432877,15,41,5.153882,1
,1.107124452,1,0,4.116839,7,54,4.841229,0,5.829129699,0,1,1.574529,4,37,4.784644,1
,3.344702727,1,1,3.640927,4,40,5.327739,1,4.106279312,0,0,1.307658,3,43,5.237828,1
,5.597583817,0,1,1.789552,5,33,5.063291,1,4.889016791,0,0,5.076378,22,62,5.309829,0
,1.614847931,1,0,1.776776,12,39,5.697535,1,5.700602099,1,1,3.228300,7,33,5.180268,1
,4.154081131,1,1,5.970791,8,64,5.303301,1,3.498311267,0,1,4.873115,16,54,4.572111,0
,1.725516742,1,1,3.473946,12,45,4.983549,1,3.088626334,0,1,1.506025,8,27,5.787719,1
,1.813557055,1,1,2.023611,5,41,4.886249,1,6.473210783,0,1,3.951521,5,45,5.585256,1
,11.904159777,0,1,3.296605,13,51,4.621678,0,4.618503669,1,1,4.554781,13,28,5.241720,1
,2.148192411,1,0,1.794118,12,46,5.109458,0,1.399641222,1,1,3.882329,5,46,6.495191,1
,6.524877740,0,1,4.092007,8,58,4.980119,0,4.257113222,0,1,2.469730,12,61,5.106757,0
,8.075994008,0,0,3.373190,7,53,4.521553,0,3.195923170,0,0,3.511813,17,54,5.196646,0
,2.868072448,1,0,3.325605,7,69,5.717564,0,0.134297954,1,0,4.100890,5,21,4.593059,1
,10.994402054,0,1,3.463139,10,52,4.650930,1,2.920086105,1,1,3.686064,10,76,4.668973,0
,7.114572805,0,1,2.859050,8,34,5.941006,1,6.318422168,0,0,3.211131,7,44,5.080005,1
,4.261228357,1,1,5.402267,17,43,5.416025,1,0.771896875,0,1,1.980577,7,43,5.966562,1
,5.793255777,0,0,7.124845,6,36,4.997398,1,1.139539301,1,1,2.473024,15,37,5.296764,1
,2.413839400,1,0,3.237792,20,64,4.650930,0,8.958817363,0,1,1.670993,10,59,4.682826,0
,2.822138019,1,0,2.963084,15,61,4.687360,0,3.422131727,0,1,3.251236,7,46,5.206833,1
,7.993613717,0,0,2.485564,6,69,4.621613,0,8.782852576,0,0,2.436971,8,35,5.315730,1
,5.397216464,0,1,3.467820,12,26,5.660932,1,2.223455060,1,0,3.055016,15,43,6.136303,1
,6.027392173,0,1,3.852236,10,34,4.881406,1,1.084586392,1,1,9.265588,28,37,5.333333,1
,1.243429923,1,1,4.679064,4,41,5.488114,1,5.668138399,0,1,1.388813,5,42,5.257000,1
,9.413060378,0,1,3.270537,11,52,5.201327,0,1.342867651,0,0,5.186753,2,41,5.327739,1
,6.845927230,0,1,1.199097,4,21,4.731417,1,0.002594154,0,1,2.510500,7,22,4.960784,1
,3.591301516,1,0,3.714610,3,43,5.735394,1,8.740492055,0,1,1.320928,5,48,4.836185,0
,4.955223384,0,1,3.600015,2,29,5.006571,1,6.020906848,0,1,3.806682,4,41,4.847189,1
,6.102748685,0,1,1.163931,3,61,5.229125,0,1.794673196,1,0,4.604695,13,36,5.077524,1
,0.961791287,1,0,1.897478,8,45,5.135196,1,4.509406570,0,1,4.138604,10,53,4.759858,0
,12.773241666,0,1,3.246316,8,46,4.459131,1,5.380294851,0,1,3.815470,5,44,5.070667,1
,2.745994109,1,1,1.887257,7,22,4.613830,1,1.070126274,1,0,4.309045,7,37,5.220239,1
,4.643592929,1,0,2.261823,5,57,4.327874,0,6.427625098,0,1,3.039427,7,39,4.605263,1
,4.353008821,0,1,2.223608,12,61,5.106757,0,9.791228457,0,1,3.375487,7,58,4.736275,0
,0.958278043,1,1,5.988311,17,46,6.027281,1,2.460672056,1,1,1.545435,8,49,5.303301,1
,3.814724112,1,1,4.389676,8,65,4.850712,0,6.472706164,1,1,4.784875,17,64,5.970874,0
,7.701176767,1,1,4.943282,8,56,4.587815,0,7.717132494,1,1,1.110856,5,34,5.350588,1
,5.253226996,1,1,7.734394,15,45,5.315730,1,2.861673863,0,1,3.736011,10,42,5.829612,1
,5.639338624,0,0,2.882569,9,50,5.423261,0,1.139180586,1,1,4.084861,10,23,4.850811,1
,1.520258869,1,1,4.470049,3,41,5.096031,1,2.082513273,1,1,4.215316,8,48,4.841229,1
,0.451145072,1,0,5.746813,23,45,4.408289,1,5.234386256,0,1,1.974097,4,53,5.359112,0
,1.842741723,1,1,3.575789,5,28,5.015292,1,2.634777668,1,0,1.380876,13,58,4.781478,0
,3.054968894,1,0,4.879149,13,32,5.344762,1,7.791163189,0,0,3.961538,13,60,4.494666,0
,7.626997128,1,1,2.394343,5,54,5.247021,0,5.744454066,0,0,3.634870,13,67,6.004324,0
,0.944391903,1,0,6.257458,13,65,4.921255,0,5.512751067,1,1,7.889996,7,55,5.698029,1
,1.473004997,1,0,3.968566,8,51,4.381244,0,3.695483268,1,1,3.887463,5,32,4.535342,1
,2.256091719,1,0,3.728736,7,56,5.366974,1,7.127344269,0,1,3.892185,4,54,4.408289,1
,5.495598537,0,1,3.465119,12,36,5.488114,1,2.024487003,1,0,3.465722,8,66,4.916011,0
,0.432877386,1,0,4.224798,2,43,5.810369,1,8.950332014,0,1,5.615929,17,54,4.577900,1
,3.217240913,0,1,3.111013,2,51,5.412659,0,7.391514326,0,0,2.190754,10,61,4.921529,0
,3.679468099,0,1,7.903886,11,26,4.871677,1,3.752650627,0,1,4.828012,13,50,5.416760,1
,1.359926640,1,1,5.950181,17,60,5.205962,0,1.779742809,1,1,2.673481,5,41,4.886249,1
,0.240097132,0,1,3.313608,17,31,6.419274,1,7.041773733,0,1,5.807772,12,56,4.984073,0
,5.274468466,0,0,4.617729,7,25,4.563989,1,8.570785937,1,1,4.003927,9,43,5.323971,1
,0.869031974,1,0,4.322735,20,53,4.843404,0,3.080014795,1,1,1.763183,8,30,6.171599,1
,3.117189568,1,1,5.752226,16,31,5.391265,1,2.512258594,1,0,4.975315,5,47,5.248639,1
,0.393262367,1,1,4.707015,20,45,5.584316,1,0.029258002,1,1,1.660248,19,22,4.983549,1
,10.925533555,0,1,1.099849,27,53,4.478945,0,11.226861385,0,0,2.292972,5,55,5.170100,0
,10.520582169,0,1,1.818042,12,44,4.830909,1,2.969630861,1,1,3.246843,4,40,5.327739,0
,5.093499793,1,1,3.942309,17,65,5.269940,0,7.722979982,0,1,2.268929,8,53,5.549887,0
,1.736996954,1,1,3.268263,27,42,5.509923,1,1.649247269,1,1,3.067915,7,49,4.693797,0
,2.631835290,1,1,4.405991,10,33,4.850712,1,11.325181177,0,1,4.104817,12,43,4.408289,1
,2.499554846,1,1,5.046849,10,39,4.641669,1,1.753958309,1,0,1.852803,13,55,4.766442,0
,10.798710395,0,1,1.132757,6,50,4.850811,0,5.488611585,0,0,7.169972,10,53,5.352583,0
,0.404717449,1,0,7.609000,3,45,5.174546,1,6.785996019,0,0,1.977703,9,58,4.706487,0
,5.457547924,1,0,2.626157,5,50,4.736275,1,6.090189523,1,1,4.484178,18,26,5.439283,1
,4.043925973,0,1,3.879262,10,57,4.985775,0,5.207490385,0,0,2.532634,3,40,4.701095,1
,6.700785566,0,0,1.046443,9,58,4.706487,0,8.063779312,0,0,3.863380,12,42,4.960784,1
,3.443084267,0,0,4.055330,13,55,6.196016,0,3.742564597,1,1,1.921520,30,43,5.659616,1
,7.347604840,0,1,2.289630,7,53,5.764246,0,4.522371497,0,0,1.743172,10,54,5.717564,0
,3.141453858,0,0,4.195614,12,59,6.130060,0,2.610096567,1,0,4.644975,10,63,5.178184,0
,5.676462759,0,1,3.145090,10,47,4.974681,0,0.823474027,0,1,1.847302,7,43,5.966562,1
,5.064240342,1,0,4.547300,20,47,5.318160,0,3.814753975,1,1,7.130591,13,48,5.404638,1
,1.948065239,1,1,2.895129,10,30,4.561979,1,0.307495726,0,1,1.999431,10,25,4.796997,1
,1.087432697,1,0,6.069058,8,39,4.956558,1,11.623339762,0,0,3.126116,16,47,4.983549,1
,5.289909514,0,1,1.343425,3,53,4.921255,0,11.200088930,0,0,1.754649,6,65,5.038911,0
,3.818835120,1,1,7.430773,18,34,6.383217,1,3.581836528,0,1,3.594143,7,57,5.625326,0
,3.978130150,0,0,4.689308,25,34,4.976703,1,3.966332634,0,0,3.086715,23,61,5.420771,0
,5.853514645,0,0,3.617398,7,44,5.080005,1,5.024868105,0,1,4.186934,13,54,4.516129,0
,5.014721549,0,1,6.599893,18,51,5.025885,0,6.123538776,0,1,9.218273,10,60,5.182124,0
,0.250535083,1,0,4.255883,17,60,5.882353,0,5.315638988,0,1,3.552472,9,39,6.211300,1
,8.943173794,0,0,2.027613,6,69,4.621613,0,0.643938181,1,1,2.772486,12,52,4.933737,0
,11.665038561,0,0,2.688076,18,57,5.420771,0,6.386061328,0,1,3.495028,13,38,5.178184,1
,3.754429619,0,1,1.990836,10,64,5.102694,0,1.803470788,1,0,3.253845,8,62,5.516086,0
,4.626745340,1,1,2.398580,13,31,5.486540,1,4.498457366,0,1,3.056049,5,56,5.943168,0
,11.781940349,0,1,4.172124,3,43,6.165568,1,0.636151513,1,0,4.104137,27,64,4.985775,0
,7.071992117,0,0,7.530330,23,31,4.819277,1,2.054706340,1,1,4.429761,12,43,5.359078,1
,4.081329354,0,0,4.672746,13,30,5.474375,1,4.410354002,0,1,1.943130,20,47,4.821142,0
,10.830422327,0,0,1.406592,6,65,5.038911,0,2.021974852,0,1,4.598517,17,30,4.839637,1
,0.984429151,0,1,3.572986,10,35,4.683626,1,2.791930722,1,1,6.938314,13,51,4.778846,0
,1.319126947,1,0,6.658133,13,47,5.487805,1,3.403215755,0,1,2.820119,10,66,4.724556,0
,1.114048386,1,1,7.170733,33,40,5.796012,1,1.726013194,1,1,1.201856,13,52,5.957490,0
,4.236758247,0,0,1.039859,7,50,5.940885,0,2.399982465,1,0,1.548747,17,60,4.980119,0
,10.192507553,0,1,2.264828,6,53,4.899540,0,4.533598947,0,1,6.302760,25,43,5.699880,1
,12.328957704,0,1,2.500607,17,65,6.250000,0,5.429531916,0,1,2.670772,5,30,4.907975,1
,1.984378344,1,0,3.649803,3,44,4.984073,1,1.674410282,1,1,8.669409,8,56,5.178184,0
,7.865770458,0,1,1.528876,5,34,4.921255,1,1.729541664,1,1,3.158925,7,55,4.967597,0
,6.968346506,0,0,3.422479,9,44,5.773003,1,8.564097952,0,1,3.386506,10,54,5.510658,0
,1.338175242,1,1,6.946359,17,27,5.624977,1,1.230913660,1,0,3.834779,12,59,4.672253,0
,10.367612136,0,0,3.435516,5,61,5.352583,0,1.067652497,1,1,2.511027,14,67,5.642155,0
,0.620071897,0,0,1.027446,13,35,5.772393,1,0.911272438,1,1,5.543977,8,40,4.902511,1
,1.303407766,1,1,4.257949,8,49,4.652018,1,5.699432411,1,1,4.889994,10,52,4.686909,0
,1.547658188,1,1,1.968839,8,55,4.930935,0,7.358874620,0,1,4.041300,20,55,5.481173,0
,5.571559299,0,1,4.007010,10,53,5.009940,0,3.676279729,0,1,1.946200,5,68,4.550068,0
,6.684024578,1,1,3.979688,17,56,4.615620,0,1.379151664,1,0,3.149651,15,61,5.455447,0
,0.629674483,1,0,2.702666,5,48,4.980119,1,4.742057996,0,0,4.425371,3,52,5.294117,0
,2.015700490,1,1,1.872647,8,48,5.590170,1,3.043686427,1,0,1.679111,7,65,5.007613,0
,6.327232707,0,0,1.423567,8,65,4.607373,0,2.035903059,1,1,3.731202,16,47,5.404638,1
,9.197425396,0,1,3.015058,7,53,4.498833,0,7.587141539,0,0,1.617228,11,34,5.381357,0
,10.418610171,0,1,2.941383,13,44,4.606335,1,7.188468185,0,0,2.067625,3,59,4.619330,0
,5.295264044,1,0,1.067770,7,43,4.733485,1,2.325087890,1,0,3.215898,13,51,5.070603,0
,6.965401833,0,0,1.519877,13,52,4.881905,0,1.138303286,1,1,1.082815,5,63,4.786756,0
,3.677413744,1,1,5.675761,17,64,4.606335,0,5.864792939,1,1,3.875541,8,42,4.423004,1
,0.760840306,1,1,3.999757,7,53,4.808812,0,3.656370708,1,0,5.468261,8,46,4.427997,1
,7.977481081,0,1,2.167356,10,53,5.555122,0,3.666981993,1,1,6.772681,12,49,4.419417,0
,8.704343821,0,0,4.100223,15,40,5.115846,1,5.110168507,0,0,1.458944,4,43,4.302066,1
,1.647651386,0,1,6.557593,5,56,4.465782,0,1.090246173,1,0,1.680488,6,42,4.899540,1
,2.257441262,1,1,3.015033,3,42,4.843404,1,1.953757578,1,1,3.024511,12,45,4.983549,1
,5.885632485,0,1,2.956574,13,43,5.479188,1,8.194390338,0,1,7.579026,13,35,5.359112,1
,10.215228869,0,0,3.524750,17,41,4.841229,1,4.927100198,0,1,2.800356,7,33,4.802921,1
,2.083129417,1,0,2.583179,7,41,4.230605,1,8.860463429,0,1,2.872429,9,57,4.811160,0
,1.122232053,1,0,3.870709,33,30,5.132883,1,5.192371923,1,1,1.290323,5,33,5.063291,1
,3.076037616,1,0,5.099107,12,28,4.741448,1,2.968698078,1,0,3.007386,12,63,5.313040,0
,6.855113339,1,1,2.081363,4,40,4.946170,1,9.073883685,0,1,1.796952,5,62,5.796012,0
,4.529704834,0,0,2.538278,9,38,4.974681,1,8.954766884,1,1,1.717818,12,30,5.174546,1
,4.204872587,1,0,5.798053,30,71,4.536092,0,0.143972740,1,0,4.512659,17,41,5.661268,1
,4.089147902,0,1,3.590890,10,54,6.084539,0,0.914337276,0,0,4.581733,10,35,4.682535,1
,1.338854033,1,1,1.164560,7,51,5.146990,0,1.504048270,1,1,4.605996,4,25,5.055576,1
,4.054685644,0,0,1.846526,8,41,5.266344,1,8.813238157,0,1,3.038083,8,59,4.680553,0
,4.712040939,0,1,2.912451,12,69,4.759858,0,1.047072230,1,0,3.203211,7,45,4.753973,1
,3.748215334,0,1,4.663353,17,46,5.115846,0,3.691430329,0,1,2.328387,13,67,4.444445,0
,1.622962498,1,1,3.417064,17,49,6.172840,0,1.703839499,1,1,3.225836,7,38,5.021689,1
,3.692680160,0,1,6.947422,8,52,5.418258,1,8.236347756,0,1,3.231294,6,35,5.514311,1
,5.493944034,0,0,1.195176,10,48,5.391265,1,1.740439738,1,0,7.561386,12,68,5.422877,0
,9.467880878,0,0,4.429035,13,47,4.811160,0,5.534484262,1,1,4.177471,8,25,4.901409,1
,3.263101583,0,1,3.750295,8,56,4.689338,1,11.031644791,0,1,3.568039,4,45,6.558120,1
,3.631417685,0,1,1.254824,12,62,5.587602,1,0.664138551,1,1,4.826736,17,22,4.923234,1
,0.551372078,1,0,2.199468,13,47,4.798963,1,1.280112444,1,1,6.409967,5,46,4.819277,1
,3.409553875,1,1,3.845603,7,38,5.454546,1,6.523484412,0,1,6.081367,7,60,4.466325,0
,12.469131390,0,1,3.269804,5,51,4.028379,0,5.796893424,1,0,2.936419,12,57,4.910347,0
,4.892863138,0,0,3.538422,13,31,4.908459,1,9.756522156,0,0,4.437509,17,61,5.994789,0
,5.956365864,0,1,3.241240,12,48,5.077524,0,4.273672079,0,1,1.138052,13,35,6.052149,1
,1.283302727,0,0,1.568443,10,75,4.771733,0,2.593753193,1,1,4.576124,12,56,4.561979,0
,11.022361688,0,1,2.562302,7,46,4.723693,1,4.341809339,0,1,1.490054,8,34,5.070667,1
,1.740294068,1,0,3.022582,4,48,4.991342,1,1.207288043,1,0,5.583565,2,41,5.327739,0
,6.801260934,0,1,1.536236,17,50,6.431975,0,1.239098875,1,1,4.007751,17,40,6.295086,1
,1.923827177,1,0,4.302293,3,53,4.835737,0,0.242325526,1,0,7.432516,17,39,5.169417,1
,1.962638927,1,0,5.646618,9,51,4.593059,0,4.765866071,1,0,3.637664,10,49,4.848485,1
,7.497251748,1,1,4.630456,20,58,4.856782,0,9.718406642,0,0,3.369178,7,58,5.153882,0
,6.159691103,1,1,4.243608,6,51,5.767761,0,8.714778531,1,1,2.099445,13,61,4.899540,0
,1.347854222,0,1,1.391289,7,35,4.977315,1,5.704179743,0,0,2.328115,5,53,4.861484,0
,6.091363504,0,1,1.109748,6,39,5.784654,1,1.326709378,1,0,3.207921,8,30,5.148021,1
,3.637595079,1,1,5.518120,13,42,5.084070,1,0.493235736,1,1,2.206559,10,33,4.681194,1
,6.604112747,0,1,3.543386,5,30,4.381244,1,0.704454326,1,0,5.358995,8,34,6.344507,1
,6.813956765,1,0,2.938351,5,54,5.247021,0,7.243775093,0,1,1.665633,12,65,4.590991,0
,6.249772609,0,0,1.794299,10,48,5.391265,1,9.608482597,0,1,1.714089,27,53,4.478945,0
,0.500085713,1,0,4.301557,15,65,6.157191,0,5.480820839,0,1,3.189965,8,54,4.850712,0
,12.710331207,0,1,5.678058,18,52,4.991342,0,2.313843835,1,1,1.115404,8,48,5.590170,1
,0.056334345,0,1,3.687239,10,70,5.263158,0,6.024738466,0,1,3.098768,10,47,4.974681,1
,1.745185673,1,0,3.256935,7,38,5.021689,1,3.839304875,0,1,3.408107,10,57,4.985775,0
,4.140271763,1,0,4.057547,13,63,4.419417,0,1.711032155,0,1,5.459517,13,61,4.374088,0
,2.945441471,0,0,1.405297,4,61,4.695976,0,3.901363912,0,1,2.377022,8,60,5.207717,0
,6.368603051,0,0,3.065122,15,63,4.943196,0,11.082273203,0,1,3.236662,8,55,5.015292,0
,6.616702935,0,0,4.697332,15,48,5.120432,0,8.089383595,0,1,2.602792,8,31,5.624463,1
,4.535427619,1,1,2.685595,6,33,5.096031,1,0.341383471,1,1,4.014767,20,45,5.584316,1
,4.216489097,0,0,1.873561,10,64,5.645998,0,10.816620562,0,1,5.293613,3,44,4.869480,1
,2.918432540,1,1,1.946789,13,45,4.847189,1,5.836709913,0,1,2.577243,9,36,5.385101,1
,2.466280534,1,0,3.519277,11,57,4.479032,0,1.710999982,1,1,2.286678,13,47,5.217020,1
,0.454346596,1,1,5.080456,23,24,5.401257,1,1.411485907,1,0,4.108376,10,45,4.921529,1
,5.112801510,0,1,1.299017,4,53,5.359112,0,2.124665475,0,1,1.207315,12,25,5.714959,1
,4.209002507,0,1,2.060803,6,56,4.540842,1,2.920484931,0,1,4.500004,33,56,4.904786,0
,4.774739271,1,0,1.534296,3,40,4.631800,1,2.784954072,1,1,4.319636,23,29,5.115846,0
,2.501403175,0,0,3.549613,10,69,4.637013,0,4.784769985,1,0,4.576474,20,42,6.797196,0
,3.484472023,1,1,2.526027,10,19,5.109458,1,0.938932350,1,0,5.619398,17,40,5.094267,1
,0.159616361,0,0,5.332134,15,38,6.033400,1,9.910881360,0,1,1.597291,7,31,5.029849,1
,3.745673862,1,1,3.402812,8,46,5.206833,0,5.694807992,1,0,1.506099,3,69,5.391265,0
,7.085666266,0,1,3.803988,15,65,5.225269,0,7.315893747,0,1,1.828972,12,65,4.590991,0
,0.359526952,1,0,7.432276,12,37,4.762347,1,3.824774877,1,0,5.918078,13,43,4.988877,1
,0.982772838,0,0,3.579282,8,73,4.666667,0,2.296799624,1,0,3.528680,7,56,5.366974,1
,1.171713557,1,1,3.774755,5,21,4.599488,1,0.403966122,0,1,5.021954,18,68,4.713064,0
,8.476544176,0,1,2.615981,6,33,4.508021,1,9.171684436,0,0,3.688894,10,50,5.132883,0
,4.501429805,0,1,3.868264,5,40,4.997703,1,7.027335286,0,1,1.930436,13,37,5.359112,1
,8.157290614,1,0,4.747985,12,50,5.163978,0,4.795064767,0,1,4.249273,10,44,4.635125,1
,0.150124681,0,0,5.595037,15,38,6.033400,1,1.685459787,1,1,5.868948,13,66,4.800717,0
,5.789678885,0,0,2.262091,5,39,4.834520,1,1.807590302,1,1,4.509097,15,42,4.879078,1
,5.849247382,1,0,4.852739,5,37,5.174506,1,6.505686522,0,0,4.137111,22,49,4.402515,1
,5.965507461,0,1,1.939071,6,39,5.784654,1,4.745028981,1,1,1.053258,18,42,5.096089,1
,7.397726906,0,1,2.246315,12,38,5.659616,1,2.582139688,1,1,6.767324,13,26,5.555122,1
,11.620713389,0,0,2.144898,11,43,5.224291,1,2.354342799,1,0,2.147087,5,38,5.151093,1
,5.919569252,0,1,5.462591,12,56,4.984073,0,2.982479436,1,0,4.404439,13,39,5.504342,1
,2.043911104,1,1,4.789044,12,68,5.286123,0,6.471085244,0,0,3.951506,21,49,5.303301,0
,0.399796909,1,1,5.164287,20,40,4.650769,1,10.989655194,0,0,3.660749,8,56,4.508264,0
,5.618858888,1,0,2.844965,5,50,4.736275,1,8.903061433,0,1,2.530546,15,61,4.723693,0
,5.849367329,1,1,5.745089,15,53,4.800717,0,3.330087391,1,1,5.408377,13,42,5.084070,1
,2.450236703,1,1,4.486764,10,42,5.766097,1,9.622271528,1,1,2.089677,7,40,5.077524,1
,7.346688972,1,1,3.503054,7,65,4.423004,0,2.478651273,1,1,4.371018,10,33,4.850712,1
,3.895491310,0,1,4.428101,7,64,4.550068,0,0.353398250,1,0,2.659858,20,46,5.386311,1
,6.413721112,0,1,5.275857,17,33,5.014839,1,11.180874380,0,1,2.211318,7,29,4.952207,1
,8.942891953,0,0,1.842506,7,53,4.536092,0,3.260151861,1,1,4.456568,8,45,5.620375,1
,0.290494463,1,0,5.976378,12,47,4.677072,0,3.055727764,0,1,6.463273,10,59,4.953681,0
,0.766194332,1,0,7.256844,4,58,4.550068,0,1.521106293,1,0,8.229116,15,50,5.385101,0
,3.305986019,1,0,4.139014,9,44,5.553775,1,6.471743709,0,1,2.100047,8,49,4.603557,1
,3.950150655,0,0,3.332262,6,53,5.326697,0,0.990544543,1,1,3.381248,8,46,4.830680,1
,3.924687900,0,0,3.280216,7,56,4.615620,0,7.015305440,0,0,3.387375,6,62,5.615465,0
,1.417217070,1,1,3.455838,13,75,5.229125,0,10.686381647,0,1,4.642725,5,62,4.231886,0
,7.922455770,0,1,2.804753,8,53,5.549887,0,6.190366155,0,0,7.111154,23,31,4.819277,1
,5.878891469,0,1,2.374204,9,38,5.137896,1,1.770688727,1,1,4.540284,17,39,5.052686,1
,2.177196351,1,0,1.940813,10,64,4.850811,0,3.530529798,0,1,6.152448,10,46,5.410267,1
,5.916002765,1,1,4.911436,10,49,4.766442,1,7.569615889,0,1,3.472873,12,61,5.182124,0
,5.298762161,0,0,1.040475,9,33,5.142595,1,1.420976909,0,1,2.265722,10,31,5.324759,1
,5.471085865,0,0,3.298073,9,52,4.923659,0,8.747444579,0,0,1.321746,3,42,4.322629,1
,1.065527489,1,0,5.262461,15,45,4.988877,1,3.061506056,1,1,2.738867,17,44,5.483719,1
,3.101768774,1,1,4.103467,5,65,4.593059,0,3.487619817,0,0,2.029529,2,41,4.245699,1
,6.933932250,0,1,3.596082,5,45,5.585256,1,1.627501819,1,1,3.103932,17,52,4.820110,0
,0.459899324,1,0,3.436148,11,53,5.064476,0,1.196476922,1,0,2.791036,16,56,4.921529,0
,11.015911395,0,1,3.866773,10,52,4.650930,1,6.493480081,0,1,1.706045,5,47,6.004324,0
,4.323063906,1,1,3.060590,7,32,5.381357,0,5.687045109,0,1,1.883918,6,47,6.994941,1
,8.029955685,0,1,1.344839,7,30,5.055576,1,4.521019032,0,1,3.643624,20,50,4.960819,0
,8.275488496,0,0,4.024176,10,44,5.476925,1,0.978927475,0,1,3.900913,10,35,4.683626,1
,4.191739835,0,0,1.734116,5,70,4.921255,0,4.186580516,1,0,1.027840,18,54,5.588332,1
,3.085856359,0,0,3.336492,10,33,5.154913,1,0.536928066,0,1,2.596819,4,32,5.487283,1
,0.194680564,0,1,6.918101,8,30,5.055576,1,3.133393842,0,1,1.236199,5,35,4.615620,1
,1.455958072,1,1,4.255443,15,72,4.896896,0,5.071150410,0,1,4.203627,10,46,5.047441,1
,5.197252612,0,1,1.523496,12,58,4.759858,0,1.686088255,1,0,6.768106,10,48,4.881406,1
,7.146015864,0,1,3.070468,8,33,5.303301,1,3.140740703,1,1,5.457124,16,31,5.391265,1
,2.995375995,0,1,2.412651,10,66,4.724556,0,5.847248092,0,1,6.272926,8,58,5.552737,0
,7.864718976,1,1,6.764246,2,33,5.962848,1,6.584225602,0,1,5.211285,13,39,4.881905,1
,0.754393587,1,0,5.989793,8,50,5.625326,0,2.387815383,1,1,4.734830,12,56,4.561979,0
,1.461969210,1,0,6.026673,7,38,4.850712,1,2.690545692,1,1,4.465125,12,55,5.000000,0
,9.542404099,0,1,1.977479,7,55,4.273348,0,6.744635861,0,1,4.647560,7,60,4.351666,0
,3.172623868,1,1,1.919412,9,37,5.423261,1,9.720524347,0,1,3.777607,5,32,5.055576,1
,3.337827875,1,0,3.729475,3,56,4.680553,0,8.633958091,1,1,2.907124,8,42,4.718646,1
,4.773444195,0,1,3.721566,13,57,4.762347,0,8.238257041,0,1,3.775267,5,63,4.518320,0
,3.198572387,0,1,4.746129,33,56,4.904786,0,3.267366762,1,1,4.905366,12,47,5.813129,1
,0.931400313,1,0,3.859220,13,32,5.659425,1,8.459602638,1,0,4.539977,12,40,4.469809,1
,10.128541824,0,1,3.368283,8,50,4.632703,0,5.195980381,0,1,1.910538,10,36,5.642155,1
,5.852506169,0,1,3.876700,7,34,5.090253,1,6.510188589,0,1,3.003975,7,57,5.063291,0
,3.872470390,0,1,4.161128,10,51,4.913402,0,7.483627495,0,0,1.357098,2,61,4.908459,0
,4.131692208,1,1,5.209945,5,54,4.561979,0,2.694968577,1,1,1.119132,7,39,5.219121,1
,3.787812692,0,1,1.064101,3,50,6.034860,0,1.253495932,1,1,3.296340,7,24,5.180798,1
,1.165154219,1,1,5.900348,17,41,5.090253,1,4.895568901,0,1,7.015869,7,33,5.352583,1
,2.806630443,0,1,2.610060,5,37,5.094267,1,1.191073802,1,0,7.813818,30,50,5.486694,1
,9.939381623,1,1,3.169415,7,56,4.242424,0,3.225863512,1,1,2.121407,10,41,5.197775,1
,0.344448391,1,1,3.092457,3,33,4.944419,1,3.726565385,0,1,2.963556,7,37,6.389871,1
,0.826608778,0,0,6.412055,13,65,4.921255,0,2.980067010,1,0,2.558614,7,49,4.252083,1
,1.236738356,1,0,7.448885,29,39,4.705882,1,3.880393843,0,1,3.072832,7,32,5.837676,1
,3.729013091,0,0,4.712902,23,66,5.153882,0,4.170337678,0,1,3.894594,7,41,5.379040,1
,9.514441755,0,1,3.584329,8,50,4.632703,0,8.285041461,0,0,4.233252,10,54,5.229125,0
,5.500631989,0,0,2.195442,7,38,4.869480,1,3.914963804,0,1,3.152215,7,64,5.135196,0
,2.088464898,1,0,4.029423,15,37,5.087983,1,4.695915733,1,0,3.674853,8,58,4.550068,0
,0.508263771,0,1,1.563390,13,62,5.781450,0,6.423408436,0,1,2.161117,12,26,4.247670,1
,3.814645745,0,0,3.108560,3,66,4.723693,0,5.121649589,1,1,5.113745,14,23,5.261336,1
,6.925514557,0,0,2.878833,5,59,4.870607,0,11.682487699,0,1,1.528803,10,52,5.357143,0
,1.416111472,1,1,4.456550,5,37,4.937851,1,3.878054673,0,1,2.235307,14,41,5.752427,1
,2.087212526,1,0,4.277159,3,63,4.549815,0,11.490053811,0,1,4.920517,22,55,4.577900,0
,1.256101345,1,1,2.184008,5,46,4.841229,1,0.994297127,0,1,3.027359,10,35,4.683626,1
,5.419649140,0,1,7.539007,15,51,5.045987,0,3.663966622,0,1,6.360583,13,53,4.633481,0
,8.512175783,0,0,4.602498,15,56,6.017664,0,2.665814427,1,1,3.777707,17,68,5.045599,0
,8.666369135,0,1,3.935311,7,50,5.000000,0,9.608806822,0,1,8.645478,10,52,4.821142,0
,12.072595912,0,1,5.177734,3,44,4.869480,1,2.886999926,1,1,5.462131,17,23,4.753973,1
,0.803767962,1,1,2.498672,13,42,4.789794,1,10.416958992,0,1,2.119693,14,66,4.364066,0
,1.026551632,1,0,3.098744,21,45,5.993707,0,7.667335959,0,1,3.372685,15,65,5.225269,0
,7.624852608,0,1,3.687177,14,36,5.345836,1,4.142351924,0,1,1.987480,6,36,4.463393,1
,3.429760437,1,1,3.259736,6,29,5.800728,1,7.490867481,0,1,1.704233,4,21,4.731417,1
,3.869903494,0,1,4.590009,10,51,4.913402,0,4.281969914,0,0,3.971873,17,36,4.761905,1
,3.944002883,0,0,3.128415,7,34,5.095541,1,10.968094598,0,0,1.468032,13,63,4.535342,0
,5.212150599,0,1,3.214282,7,46,5.512261,1,3.622001163,0,1,2.950828,7,37,6.389871,1
,0.385548345,1,1,1.019563,12,42,5.732484,1,1.914268609,1,0,4.715359,3,53,4.835737,0
,0.588104360,1,1,9.685002,10,38,4.781461,1,5.557568658,0,1,2.078265,6,37,4.923234,1
,3.822249042,0,1,1.666467,3,37,4.766442,1,10.318248702,0,1,4.334759,5,62,4.231886,0
,2.823237952,1,0,2.649425,3,49,5.624385,0,2.492673108,0,0,3.377316,20,64,4.650930,0
,1.269897482,1,1,4.920774,13,33,5.045987,1,9.357387723,0,1,1.830246,7,31,5.029849,1
,1.645688168,1,1,3.882723,13,66,4.550068,0,4.903651221,0,1,4.779042,8,49,4.503865,0
,4.510936888,1,1,2.734384,13,33,5.145105,1,0.956381372,1,0,6.548832,7,54,5.517594,1
,1.914477820,1,1,3.041214,18,50,5.257000,0,4.989703569,0,1,4.949574,3,23,4.655240,1
,8.088284030,0,1,1.665561,2,61,4.908459,0,5.214449620,0,1,1.544228,7,59,4.801516,0
,0.561215906,1,1,6.643541,8,23,5.519851,1,7.285200687,0,1,2.316143,12,51,5.257000,0
,1.658712539,1,1,2.431853,10,46,4.535342,0,8.969450859,0,1,7.833928,14,52,4.701095,0
,1.862072924,0,1,3.555554,12,66,4.759858,0,0.624934128,0,1,1.302327,8,26,4.827945,1
,1.355911583,1,1,4.869751,10,23,4.850811,1,8.520290436,0,1,3.430283,14,44,5.085716,1
,9.497624757,0,1,2.607802,8,67,4.815713,0,0.984252049,1,1,4.560475,15,44,5.948074,1
,6.512466904,0,1,3.379775,3,51,4.772126,0,9.250421339,0,1,1.691052,15,33,5.401378,1
,3.061095497,1,0,3.273423,10,64,5.357143,0,8.622636667,0,1,1.536653,7,52,5.201327,0
,6.451470142,0,1,3.627664,12,34,5.249339,1,10.723622760,0,1,2.081499,17,61,4.714770,0
,4.626932916,0,1,3.223785,3,37,5.263158,1,1.174603368,1,0,3.604351,7,39,5.823232,1
,3.277367725,0,0,4.956886,15,40,5.209758,1,4.576979733,1,1,3.188168,14,47,5.550788,1
,1.751395382,1,1,6.483675,18,68,4.672253,0,3.228455601,1,1,4.016816,10,66,4.350764,0
,8.824059396,1,1,1.913446,17,64,4.533199,0,6.269713275,0,0,6.857850,10,67,5.624713,0
,6.376485793,0,0,3.339430,15,40,4.719673,1,3.583907430,1,0,3.372652,20,60,5.611407,0
,1.034547914,1,1,4.694744,30,43,5.518136,1,5.683106629,0,1,1.782150,5,34,5.695211,1
,9.372136899,0,0,4.141078,15,51,5.106757,0,8.749752051,0,1,2.565810,7,52,5.095541,0
,2.875352299,1,0,5.988152,10,46,5.555451,1,1.999597713,1,0,7.011769,5,28,5.025885,1
,3.684788517,1,0,4.480205,13,32,5.455447,1,11.260781370,0,0,1.302029,7,46,4.933737,1
,2.667858936,1,1,7.545183,5,36,4.577911,1,0.165100834,0,1,5.976905,15,38,6.033400,1
,3.677171814,1,1,4.256848,13,66,3.875617,0,7.959013539,0,0,1.518858,10,66,4.686909,0
,0.453267962,0,1,1.311068,13,62,5.781450,0,7.981381743,0,0,3.079583,1,48,4.952456,1
,1.091764040,1,0,6.968320,12,45,5.142595,1,2.530478925,1,0,3.404665,22,68,5.164568,0
,7.857383113,0,1,1.535848,11,55,5.334129,0,6.030819284,1,1,4.804286,8,25,4.901409,1
,3.541964972,0,1,1.563701,11,41,5.412659,1,8.166674885,0,1,1.516168,8,45,5.488113,1
,1.452840143,1,1,3.202186,8,33,4.770898,1,0.627102346,1,0,1.644970,5,62,5.077524,0
,2.017670790,1,1,2.584701,7,51,5.034317,1,3.913311846,0,0,1.630889,15,29,5.161291,1
,12.363519356,0,1,4.134158,8,43,4.575657,1,2.280510671,1,0,2.942188,17,59,4.796997,0
,4.341179063,1,0,5.692550,9,29,4.561979,1,0.543884458,1,0,5.872703,9,37,5.517594,1
,5.384412364,0,0,7.086228,10,53,5.352583,0,0.357926877,0,0,1.233673,10,41,4.631800,1
,7.515937889,0,1,3.709324,7,57,4.841229,1,6.610268430,1,1,2.985425,12,26,4.247670,1
,0.767447666,1,0,3.118447,7,30,5.295317,1,2.430434981,1,0,2.116950,13,40,5.014839,1
,8.792517293,0,0,1.345875,10,52,4.454354,0,3.013777351,1,1,3.569921,13,36,4.532735,1
,2.163265932,1,0,4.411266,7,59,5.229125,0,7.968264872,0,0,3.003258,12,42,4.960784,1
,0.844040503,1,1,8.385769,7,42,6.516221,1,2.132462188,1,1,2.068871,10,30,4.561979,1
,10.138973269,0,0,1.010718,7,52,4.798963,0,3.056468891,1,1,1.038576,7,58,5.389681,0
,2.560767053,1,1,6.056348,5,52,4.976703,0,3.663193185,1,0,3.686986,7,31,6.202187,1
,2.291634648,1,1,4.039808,13,41,4.991342,1,8.643643262,1,0,2.566005,4,49,5.474375,0
,9.164207891,1,1,3.264355,14,29,4.820110,1,4.102095343,0,1,1.825318,3,47,6.041007,1
,7.108724919,0,0,1.307443,8,50,5.270361,0,4.659110887,0,0,2.996149,5,63,4.375697,0
,3.004592075,1,0,2.455574,10,44,4.724556,1,8.739417505,0,1,1.396395,7,23,5.228350,1
,7.454305665,0,1,4.503241,8,31,5.164568,1,4.439855779,0,0,1.313272,8,41,4.635125,0
,7.439941799,0,1,3.301486,17,56,4.615620,0,4.982958377,1,0,5.707454,13,29,5.220239,1
,8.799376889,0,0,1.639104,3,54,4.810457,0,1.478564699,1,1,3.254479,5,46,6.495191,1
,8.231881752,0,0,4.892462,13,59,5.576314,0,1.370228371,1,1,5.351432,17,41,5.090253,1
,7.847733076,0,1,1.340585,11,55,5.334129,0,9.260718383,0,0,2.958745,8,67,5.078968,0
,0.713335555,1,0,7.859836,11,44,4.311743,1,1.908614536,1,1,2.698272,9,48,4.704970,0
,0.947645353,1,1,4.364577,9,37,4.276668,1,7.067428947,0,1,3.923005,12,44,5.497474,1
,8.226063233,0,0,4.746083,15,56,6.017664,0,2.140435100,1,1,6.296405,15,72,4.615620,0
,8.220473873,0,1,3.105577,12,61,5.182124,0,3.949257197,1,1,6.190143,20,33,5.429166,1
,7.047474298,1,1,1.336732,7,48,5.416645,0,7.598319805,0,1,3.714017,7,54,5.481173,0
,10.267383089,0,1,1.222931,15,33,5.401378,1,4.023251093,0,0,1.371423,6,52,4.599488,0
,8.131350075,0,0,2.936749,7,76,4.621450,0,4.243175326,0,1,5.441902,11,40,4.923659,1
,1.052185439,1,1,2.927342,15,37,5.296764,1,11.637009767,0,1,4.994733,22,55,4.577900,0
,1.843750316,1,1,5.343756,12,51,5.421687,1,3.217178057,0,1,2.654884,13,55,4.166667,0
,0.050499688,0,1,3.330976,10,70,5.263158,0,2.150433838,1,1,4.745550,12,43,5.359078,1
,4.119472797,0,1,3.218857,14,41,5.357143,1,5.816827159,0,1,1.338042,17,29,5.386785,1
,1.235188234,0,1,5.678901,13,51,4.965363,0,1.902475299,1,1,8.224553,21,44,6.073310,1
,3.472140973,0,1,4.710552,6,39,5.299210,1,2.225435008,0,1,3.108445,3,67,4.224999,0
,3.769734310,1,1,3.474964,14,64,5.295317,0,10.907262619,0,1,8.869520,17,33,4.933737,1
,9.280621948,0,1,3.066142,8,50,4.632703,0,5.440011473,0,1,3.309324,10,55,4.892449,0
,0.720430845,0,1,1.952454,10,35,5.661270,1,9.275824202,0,1,1.432103,7,46,4.821142,1
,6.171754285,0,1,2.782620,5,30,4.907975,1,0.851761242,0,1,4.480100,9,37,4.276668,1
,3.991479129,0,1,2.788143,8,60,5.207717,0,4.027130766,0,1,2.712152,6,56,4.540842,1
,2.150557449,1,1,4.657102,12,31,5.333333,1,8.876155198,0,1,1.333591,7,55,4.273348,0
,2.709879382,1,1,7.772380,5,36,4.577911,1,3.509888648,1,1,4.911962,8,65,4.850712,0
,1.438177766,0,1,6.566161,3,40,5.164568,1,3.654998218,0,1,1.397364,6,50,4.718646,0
,0.491549121,1,1,7.104059,22,57,4.821142,0,0.643228080,1,0,3.836880,7,64,5.028178,0
,5.935350567,1,1,1.242361,8,26,4.960784,1,3.294179628,1,1,6.106902,3,66,4.398887,1
,3.222286648,1,1,1.672454,8,43,4.899540,1,3.880637595,1,1,1.478945,20,37,5.153882,1
,5.916463291,1,1,1.418330,9,53,5.182124,1,8.895888735,0,0,4.131253,23,36,4.930935,1
,7.569828171,0,0,4.078775,15,56,6.017664,0,0.190022413,0,1,3.969467,5,52,5.318160,0
,6.571915947,1,1,3.711871,20,55,5.167555,0,5.651828974,1,0,2.605541,8,68,4.466325,0
,2.792770935,1,0,3.232202,7,61,5.376453,0,2.112533086,1,1,4.284349,13,36,5.077524,1
,6.417937165,0,0,3.174980,7,36,5.159393,1,10.506763444,0,0,3.720029,7,58,5.153882,0
,6.441015197,1,0,1.387931,7,46,5.153882,0,0.416201286,1,0,7.575703,3,45,5.174546,1
,3.400599978,0,1,4.229721,13,53,4.516129,0,2.257456296,1,0,4.318870,13,52,4.839637,0
,1.980817692,1,0,3.172698,7,46,5.229125,1,0.409135638,0,1,1.124168,9,39,5.767761,1
,9.626953711,0,1,4.196797,5,62,4.231886,0,8.873726753,1,1,4.295763,5,29,5.907148,1
,6.526882130,0,0,1.103429,4,62,5.323971,0,2.033732309,1,1,4.651381,17,42,4.800717,1
,3.246602258,1,1,3.451903,15,40,4.891389,1,3.362305370,0,1,4.253282,7,48,5.911692,0
,5.530887640,0,1,1.786893,3,53,4.921255,0,6.631603905,0,1,7.327827,7,44,5.115846,1
,7.165996683,0,1,4.614974,12,42,6.059600,1,8.063983602,0,1,4.630703,13,25,5.363205,1
,5.704537372,0,1,4.382012,20,54,5.291503,0,4.686644559,0,1,2.653673,5,63,4.781461,0
,6.049204838,0,1,2.788475,16,61,4.830680,0,1.036582874,1,0,6.220311,10,46,4.194352,1
,11.037226142,0,1,1.951964,17,56,4.705882,0,5.994978555,0,1,2.766535,6,37,4.923234,1
,4.072662914,1,1,2.209646,6,33,5.096031,1,1.261901212,1,1,4.632160,17,40,6.295086,1
,6.329080725,0,1,4.307175,12,34,5.590170,1,0.657757057,0,0,2.297450,17,66,5.714959,0
,10.014428743,0,1,4.216938,5,62,4.231886,0,9.308329221,0,1,3.008391,14,44,5.085716,1
,2.229046193,1,1,7.814473,5,28,5.025885,1,11.016614279,0,1,2.086410,3,56,4.902511,0
,12.453829413,0,1,4.405956,22,55,4.577900,0,1.480791978,1,0,4.422947,10,51,4.374088,0
,6.073305339,1,1,3.013016,8,42,4.423004,1,10.721930187,0,1,2.417753,10,59,4.921255,0
,2.466598504,1,0,4.432532,7,41,6.171599,1,3.503739686,0,1,1.275846,8,40,5.068487,1
,11.649848491,0,1,4.518982,7,52,5.062724,1,0.582062522,1,0,4.834595,15,65,6.157191,0
,3.602467462,0,0,3.425201,17,54,5.196646,0,9.249038990,0,0,1.177843,5,22,4.881406,1
,3.963687505,0,1,3.850679,17,40,4.860499,1,4.304977773,1,1,6.796459,20,33,5.429166,0
,6.082609915,1,1,1.356877,9,53,5.182124,1,3.009463902,0,1,1.670912,6,50,4.718646,0
,4.757404046,1,0,2.645656,5,57,4.327874,0,3.951021254,0,1,3.233695,17,40,4.860499,1
,6.717183333,0,1,1.588997,7,30,5.055576,1,5.538491000,0,1,3.475693,8,56,4.231886,0
,0.329599467,1,0,6.217489,8,67,4.841229,0,3.427751670,1,1,3.416616,12,59,4.333398,0
,5.739149143,0,0,2.804315,20,59,5.656854,0,4.977774590,0,0,3.624645,3,62,4.631770,0
,5.445877095,0,1,3.540264,13,31,4.908459,1,10.177229403,0,1,3.953893,8,50,4.632703,0
,2.568706762,1,0,3.587995,20,64,4.650930,0,4.019721608,0,1,3.651225,25,52,6.521562,0
,8.533050785,0,1,1.084226,17,73,4.878049,0,1.025160758,0,1,3.647254,10,35,4.683626,1
,4.303917603,0,1,2.133479,3,31,4.594660,1,12.999218050,0,0,1.566818,5,70,5.043558,0
,9.677741577,0,0,1.778715,9,69,5.790636,0,3.838812548,0,1,4.145922,17,46,5.115846,0
,6.949571727,0,0,1.948355,20,34,5.397807,1,9.783264891,0,1,3.596454,8,41,6.014000,1
,4.379084150,0,1,2.724802,6,56,4.540842,1,3.014898616,1,1,3.189670,5,34,5.421687,1
,1.772994029,1,1,3.198882,7,55,4.967597,0,6.667516183,0,1,5.050620,13,33,4.953681,1
,7.165767668,0,0,2.239096,10,61,4.921529,0,2.833907426,1,0,2.945502,15,61,4.687360,0
,7.983392664,0,1,1.161906,8,45,5.488113,1,6.245410346,0,1,1.657422,8,49,4.977630,0
,6.186331715,0,1,1.869346,5,65,5.045987,0,3.357774321,1,0,1.963942,12,68,4.480820,0
,1.383746243,1,0,1.952852,7,55,4.624277,0,1.362630493,1,1,4.468038,10,27,5.015566,1
,6.927369494,0,1,1.243462,12,56,5.154913,0,2.503036259,1,1,1.514198,6,25,4.556451,1
,7.275876713,0,1,1.252355,7,43,4.562997,1,11.576856547,0,0,1.411803,6,65,5.038911,0
,6.818662056,0,1,3.567689,10,60,5.295317,0,3.423068071,0,1,1.089534,6,50,4.718646,0
,1.953605140,1,1,2.907695,18,50,4.535342,0,5.727678074,1,1,3.232356,7,38,5.138322,1
,12.132741075,0,1,3.095812,13,47,5.957490,1,6.053348652,0,1,1.128622,8,49,4.977630,0
,3.694768409,0,1,4.871262,16,54,4.572111,0,6.253365566,0,1,1.742508,17,39,5.500175,1
,3.705009777,0,1,3.267305,10,66,4.622501,0,2.938119453,1,1,5.360723,17,23,4.753973,1
,2.135541826,1,1,3.526161,8,55,5.370431,0,5.668259481,1,0,2.574615,10,29,5.474375,1
,7.080065233,0,1,1.692322,5,34,4.921255,1,8.544738114,0,0,2.767920,10,39,5.423261,1
,7.045502733,0,1,1.396488,5,39,5.696002,1,0.002559307,0,1,3.064055,6,37,4.976703,1
,6.617804649,0,0,3.760495,20,39,4.493895,1,9.681871105,0,1,4.324311,4,45,4.851086,1
,6.044649776,0,1,3.856579,7,63,4.374999,0,6.420862876,1,1,4.622526,10,49,4.766442,1
,1.874660804,1,1,1.159535,9,44,5.059026,1,1.235709549,1,0,5.039398,12,44,5.904718,1
,3.258418523,0,0,5.037710,8,28,4.893999,1,4.830709502,1,0,4.918688,15,54,5.326697,0
,2.922021965,1,0,1.920790,10,34,5.266344,1,3.906092804,0,0,4.699752,11,48,4.869480,0
,0.697623587,1,0,4.116854,5,59,5.642155,0,9.941569914,0,1,8.289183,10,52,4.821142,0
,10.350027004,0,1,3.366197,8,55,5.015292,1,5.527569011,1,1,4.051165,18,26,5.439283,1
,2.606175378,1,1,6.641173,13,26,5.555122,1,2.026639052,1,1,1.038067,10,54,5.294117,1
,7.008150992,0,0,4.232506,12,42,6.059600,1,2.646393091,1,1,6.246285,18,54,5.661270,1
,0.173419837,0,1,4.740231,15,49,4.967444,0,1.716937293,1,0,3.299435,3,25,5.451704,1
,2.830695509,1,0,2.794498,8,68,5.121871,0,3.569302815,1,1,3.840850,7,54,5.294117,0
,4.225227549,0,0,4.298415,10,41,5.021689,1,10.123232628,1,1,3.056420,21,67,4.610694,0
,4.584152274,1,1,3.654068,7,48,5.770498,0,6.633500108,0,0,6.271458,7,60,4.466325,0
,1.716644522,1,0,2.746618,9,43,5.484352,1,2.479548329,1,1,4.672449,10,33,4.850712,1
,1.574094294,0,1,3.839817,5,29,5.407597,1,4.376846888,0,0,1.434719,8,49,4.948717,1
,2.838176338,1,0,5.180557,10,49,5.497474,0,1.597733354,1,0,4.702143,20,69,4.319955,0
,7.761047735,0,1,3.263322,8,33,5.295317,1,2.303782009,1,1,4.107277,12,53,5.078968,0
,8.779316914,1,0,2.786833,4,49,5.474375,0,0.707945846,1,0,6.989046,15,30,5.187748,1
,2.771433325,1,1,5.418486,17,52,5.141796,0,0.576891829,0,1,3.060024,7,60,4.841229,0
,3.349106529,0,0,1.849925,8,60,5.094267,0,1.159469227,1,1,5.524873,8,43,5.063291,1
,1.322922594,1,0,2.465975,8,41,4.718646,1,0.517570282,1,0,3.848236,7,36,5.447472,1
,8.030447971,0,1,1.702648,12,65,4.590991,0,0.213668209,0,1,4.255354,13,61,6.373774,0
,3.283146825,1,0,2.082126,7,54,5.677647,0,1.271388416,1,1,3.141950,10,36,5.333006,1
,7.194447942,0,0,3.324746,1,48,4.952456,1,8.775137830,1,1,3.302513,7,56,4.242424,0
,11.477803069,0,0,4.508761,6,53,4.952207,0,4.495830518,0,1,4.465289,3,45,5.590170,1
,0.842441077,1,0,5.454550,8,50,5.625326,0,10.128304282,0,0,2.302096,11,43,5.224291,1
,5.414548830,1,1,5.351232,10,23,5.132883,1,1.672620775,1,0,1.674226,3,62,4.997703,0
,1.865278895,1,1,5.616745,6,66,5.111615,0,10.968824176,0,1,1.261960,13,57,4.154942,0
,4.239609137,0,1,1.289914,9,58,5.060192,0,6.409953923,1,1,3.524953,10,67,4.419417,0
,1.703956304,1,1,2.970582,5,53,5.521473,0,6.551169501,0,1,1.707130,4,47,5.552737,1
,10.016568240,1,1,2.170398,8,52,5.007613,0,4.105695818,0,0,1.727860,6,52,4.599488,0
,11.702959931,0,1,1.333864,15,41,5.090253,1,0.082156753,1,1,3.966131,23,62,4.635125,0
,2.228627842,1,0,3.881510,7,56,5.366974,1,0.994236710,1,1,6.727516,7,54,5.517594,1
,1.967989864,1,1,3.760652,3,67,4.224999,0,6.943112649,0,0,5.903703,15,39,5.161291,1
,5.752026284,1,1,3.380272,8,42,4.423004,1,3.066375134,0,1,1.500041,5,55,4.668973,0
,5.137912317,0,1,1.842717,9,33,5.219121,1,0.188009532,0,1,6.516876,8,30,5.055576,1
,8.457299556,0,1,4.308274,5,32,5.625000,1,5.110364498,0,1,6.847764,20,51,4.913402,1
,0.432127170,0,0,1.126517,10,41,4.631800,0,5.879574224,0,1,1.449354,5,43,4.864693,1
,4.086557800,0,1,3.752157,23,61,5.420771,0,6.479730458,0,1,1.541884,6,33,5.474375,1
,2.096060458,1,1,4.820653,13,42,6.016540,1,8.740526201,0,1,2.848737,4,54,4.960819,0
,0.987269106,1,1,1.953991,5,55,4.933303,0,3.784784833,0,1,1.146611,13,48,4.960784,0
,7.895338056,0,0,3.431830,10,39,4.519892,1,0.153081605,0,1,5.174208,10,21,5.420764,1
,3.612458119,1,1,3.034005,6,31,5.128117,1,1.204751877,1,0,4.917909,25,53,4.631800,0
,4.797438903,0,1,3.149489,3,43,5.050762,1,3.653500320,1,0,3.462682,3,43,5.735394,1
,4.572908745,0,1,1.790413,5,33,5.115846,1,8.925976739,0,1,1.753543,8,43,5.010377,1
,9.424089447,0,1,3.067573,7,48,4.808812,0,2.721229521,1,0,3.269509,22,68,5.164568,0
,0.617500623,1,0,2.210704,13,47,4.798963,1,0.733668341,1,1,3.274835,10,30,4.998959,1
,1.777653448,1,0,1.163099,13,55,4.766442,0,8.045365441,0,0,4.279056,10,63,5.621055,0
,0.895810048,1,1,6.392544,23,70,4.983549,0,5.886969043,0,0,3.685820,7,40,4.960784,1
,7.924019439,0,0,5.793592,5,38,5.697535,1,11.543164037,0,1,4.212662,5,45,5.329681,0
,2.092920364,1,0,4.339656,13,52,4.839637,0,5.835997393,0,0,5.470233,5,58,4.302066,1
,7.597981449,0,1,3.128263,12,34,4.966996,1,7.735710616,1,0,2.544029,4,49,5.474375,0
,7.710111583,0,1,2.007118,7,48,4.366659,1,3.082220332,1,1,5.178506,11,69,5.112992,0
,6.819277120,0,0,5.871234,15,44,5.025885,1,6.581073376,1,0,3.642454,13,59,4.923659,0
,1.263674974,1,0,6.885249,5,42,4.535342,1,9.058271228,1,1,1.749957,13,34,4.465782,1
,1.561582213,1,0,3.557060,4,48,4.991342,0,5.342270791,0,1,7.101077,15,51,5.045987,0
,7.713528057,0,1,3.927043,7,55,5.062724,0,5.871852210,1,1,3.520430,10,44,4.800717,1
,9.344645272,0,1,1.791652,13,56,5.993707,0,4.890129012,0,1,2.478405,8,38,5.229125,1
,1.194444507,1,0,3.177897,8,54,5.120809,0,10.277324399,0,0,1.573010,13,59,5.809277,0
,2.300039799,1,1,4.265143,8,45,4.724556,1,6.997014031,0,1,2.162587,7,53,5.764246,0
,3.270359373,0,1,6.658219,13,53,4.633481,0,0.593110196,1,0,4.767141,15,65,6.157191,0
,4.606589056,1,0,2.608047,5,57,4.327874,0,3.920409560,0,0,4.776128,23,66,5.153882,0
,4.651463728,0,0,1.835911,5,70,4.921255,0,3.133126574,1,0,5.083054,12,28,4.741448,1
,0.169656895,1,1,4.315575,8,53,5.201327,1,6.501047844,0,1,3.243885,5,43,4.984073,1
,8.068813673,0,1,3.223044,7,58,4.736275,0,7.840334027,0,0,3.619049,17,41,4.376881,1
,3.098918700,0,1,6.186871,10,59,4.953681,0,1.928365503,1,0,3.344411,8,66,4.916011,1
,0.764759076,1,0,7.455564,11,44,4.311743,1,1.757825223,1,0,4.464654,20,69,4.319955,1
,7.459053005,1,0,3.732241,21,49,5.070667,0,5.689613634,0,0,3.397186,9,52,4.923659,0)
testhare <- matrix(testhare,ncol=8,byrow=TRUE)
################################################################################
