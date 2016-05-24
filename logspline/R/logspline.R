#
# Copyright [1993-2016] [Charles Kooperberg]
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
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
        PACKAGE = "logspline")
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
        PACKAGE = "logspline")
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
    cc2[1] <- 5/0
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
        PACKAGE = "logspline")
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
      PACKAGE = "logspline")
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
      PACKAGE = "logspline")

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
        PACKAGE = "logspline")
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
        PACKAGE = "logspline")
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
