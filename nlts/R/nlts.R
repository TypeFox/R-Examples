##############################################################################################
"contingency.periodogram" <- function(x, maxper = 6){
##############################################################################################
#contingency.periodogram is a funcion to estimate the contingency periodogram
#of Pierre Legedre and Pierre Dutielle to test for periodicity in categorical time series.
# I have coded it so as to provide both the Fisher exact test and the asymptotic chi-square.
#
#REFERENCE
#Legendre et al. (1981) The contingency periodogram: A method for identifying rhytms in
#   series of nonmetric ecological data. Journal of Ecology, 69, 965-979.
#
#REQUIRED ARGUMENTS
#x         vector of length n representing the categorical time series.
#maxper    the maximum lag (period) considered
#
#VALUE
#an object is returned consisting of the following components:
#the fischer exact test at each lag,
#the asymtpotic chi-square value, associated d.f. and asymptotic p.val at each lag
############################################################################################
    t1 <- as.factor(x)
    n <- length(x)
    kji <- matrix(NA, nrow = maxper, ncol = 4)
    dimnames(kji) <- list(c(1:maxper), c("exact.p", "chi-val", "df", "asympt.p"))
    for(i in 2:maxper) {
        t2 <- t(1:i)
        kast <- (as.integer(n/i)) * i
        t3 <- cbind(as.factor(t1[1:kast]), as.factor(t2))
        kji[i, 1] <- as.vector(fisher.test(table(t3[, 1], t3[, 2]))[1])$p.value
        t4 <- chisq.test(table(t3[, 1], t3[, 2]))
        kji[i, 2] <- as.vector(t4$statistic)
        kji[i, 3] <- as.vector(t4$parameter)
        kji[i, 4] <- as.vector(t4$p.value)
    }
    res <- as.matrix(kji)
    class(res) <- "contingency.periodogram"
    res
}

##############################################################################################
"plot.contingency.periodogram"<-function(x, ...){
##############################################################################################
#this is the generic plot function for contingency.periodogram objects
##############################################################################################
    object<-x
    n <- length(object[, 2])
    plot(c(2:n), as.vector(object[2:n, 2]), type = "b", xlab = "period",
        ylab = "chi-val")
    lines(2:n, qchisq(0.95, 1:(n - 1)))
}

##############################################################################################
"lin.order.cls" <- function(x, order = 1:5, n.cond = 5, echo = TRUE){
##############################################################################################
#lin.order.cls is a function to estimate the order of a time series using cross-validation of
#the linear autoregressive model. Coefficients are estimated using conditional least squares.
#
#DETAILS
#The time series is normalized prior to cross-validation.
#Note that if the dynamics is highly nonlinear, the nonparametric order selectors (ll.order or
#   nw.order) may be more appropriate.
#If dynamics is approximately linear lin.order.mle will have greater power (although it is
#   a bit slower). I have coded this function to use for comparison with the nonparametric
#   order selectors, because those uses (nonlinear) conditional least-squares.
#
#BACKGROUND
#I coded this functions to estimate the order of ecological time series.
#Bjornstad et al. (1998; Journal of Animal Ecology (1998) 67:110-126; see
#               http://onb.ent.psu.edu/publ/plodia/plodia.pdf)
#
#REQUIRED ARGUMENTS
#x        is a vector representing the time series. The data is scaled to
#       unit variance prior to estimation
#order    a vector representing the orders to be considered. The default is 1:5
#n.cond   is the number of observations to condition on (must be >= maximum order).
#       The default is 5.
#echo     if TRUE a counter over the data points and the order
#               is produced to monitor  progress.
#
#VALUE
#an object of class lin.order.cls is returned consisting of the following components:
#order     represents the search grid of orders
#CVd       is the cross-validation errors across the grid of orders
############################################################################################
    ans <- as.data.frame(matrix(NA, ncol = 2, nrow = length(order)))
    names(ans) <- c("order", "CVd")
    ans[, 1] <- order
    n <- length(x)
    p <- length(order)
    cvseries <- (x - mean(x[(n.cond + 1):n]))/sqrt(var(x[(n.cond + 1):n]))
    xmat <- matrix(0, (n - n.cond), (p + 1))

    for(i in 1:length(order)){
        xmat[, i] <- cvseries[(n.cond + 1 - order[i]):(n - order[i])]
    }

    xmat[, (p + 1)] <- cvseries[(n.cond + 1):n]

    for(j in 1:p) {
        cv <- 0
        conv <- 0
        for(i in 1:(n - n.cond)) {
            dati <- xmat[ - i,  ]
            coef <- solve(t(dati[, (1:j)]) %*% dati[, (1:j)]) %*% t(dati[, (1:j)]) %*% dati[, (p + 1)]
            tpred <- t(xmat[i, 1:j]) %*% coef
            cv <- cv + (tpred - xmat[i, (p + 1)])^2
            if(echo){ cat(i, " ")}
        }
        if(echo){cat("\n order ",j, " of ",p,"!\n")}
        ans[[2]][j] <- cv/(n - n.cond)
    }
    class(ans) <- "lin.order"
    ans
}


##############################################################################################
"summary.lin.order" <- function(object, ...){
##############################################################################################
#this is the generic summary function for lin.order objects
##############################################################################################
    max <- object$order[order(object$CVd)][1]
    cat(paste("The estimated order is ", max, " with a cross-validation error of ", round(
        object$CVd[order(object$CVd)[1]], 3), "\n\n", sep = ""))
    out <- data.frame(order = object[[1]], CVd = object[[2]])
    out
}

##############################################################################################
"plot.lin.order" <- function(x, ...){
##############################################################################################
#this is the generic plot function for lin.order objects
##############################################################################################
    object<-x
    plot(object$CVd ~ object$order, type = "b")
}

##############################################################################################
"ll.order"<- function(x, order = 1:5, step=1, deg = 2, bandwidth = c(seq(0.3, 1.5, by = 0.1), 2:10), cv=TRUE, echo=TRUE){
##############################################################################################
#ll.order is a function to estimate the order of a time series using the nonparametric order
#selection method of Cheng and Tong (1992, 1994) as modified by Yao & Tong (1994;
#see also Fan, Yao & Tong 1996). The method uses leave-one-out cross-validation of the
#locally linear regression against lagged-abundances.
#
#DETAILS
#The time series is normalized prior to cross-validation
#A Gaussian kernel is used for the local regression.
#The bandwidth is optimized using crossvalidation or GCV. If a single bandwidth
#is provided, no cross validation of bandwidth will be carried out.
#Highly nonlinear data will require more narrow bandwidths.
#
#If NA is returned it may be because the min bandwidth considered is too small
#relative to density of data.
#
#STATISTICAL REFERENCES
#Cheng and Tong (1992) On consistent nonparametric order determination and chaos.
#   J Royal Statist Soc B, 54, 427-449.
#Cheng and Tong (1994) Orthogonal projection, embedding dimension and sample size in
#   chaotic time series from a statistical perspective.
#   Phil Trans R Soc Lond A, 348, 325-341.
#Yao and Tong (1994) Quantifying the influence of initial values on non-linear
#   prediction. J Royal Statist Soc B, 56, 701-725.
#Fan, Yao and Tong (1996) Estimation of conditional densities and sensitivity
#   measures in nonlinear dynamical systems. Biometrika, 83, 189-206.
#
#BACKGROUND
#Wilhelm Falck and I coded these functions to estimate the order of nonlinear
#ecological time series.
#Bjornstad et al. (1998; Journal of Animal Ecology (1998) 67:110-126; see
#               http://onb.ent.psu.edu/publ/plodia/plodia.pdf)
#Stenseth et al. (1998; Proceedings of the National Academy of Science USA
#       (1997) 94:5147-5152; see
#       http://onb.ent.psu.edu/publ/hare/hare.pdf)
#
#I'll try to keep an updated set of code at
#http://onb.ent.psu.edu/
#
#REQUIRED ARGUMENTS
#x         is a vector representing the time series.
#order     a vector representing the orders to be considered. The default is 1:5
#bandwidth  a vector representing the grid of bandwidth to be considered for the Gaussian
#           product kernel.
#
#VALUE
#an object of class ll.order is returned consisting of the following components:
#grid      is the cross-validation errors across the grid of orders and bandwidths
#order     is the raw grid of orders
#CV    is cross-val error
#GCV       is Generalized CV
############################################################################################
    require(locfit)

    res<-as.data.frame(matrix(NA, ncol = 6, nrow = length(order)*length(bandwidth)))
    names(res) <- c("order", "CV", "GCV", "bandwidth", "df", "GCV.df")

    bogrid<-expand.grid(bandwidth, order)
    res$order<-bogrid[,2]
    res$bandwidth<-bogrid[,1]

    T <- length(x)

    cvseries <- (x - mean(x[(max(order) + 1):T]))/sqrt(var(x[(max(order) + 1):T]))
    ldata<-mkx(cvseries, order+step-1)

    n<-dim(ldata)[1]

    for(i in 1:dim(bogrid)[1]){
    tmp<-NULL
    if(cv == TRUE){
        tmp<-locfit.raw(lpx(ldata[,1:(bogrid[i,2])], deg=deg, h=bogrid[i,1]),
                y=ldata[,(length(order)+1)], kern='gauss', ev=dat(cv=TRUE))
        res$CV[i]<--2*tmp$dp["lk"]/n
        if(res$CV[i]==0) res$CV[i]<-NA
        res$df[i]<-tmp$dp["df1"]
        }

        tmp<-locfit.raw(lpx(ldata[,1:(bogrid[i,2])], deg=deg, h=bogrid[i,1]),
                y=ldata[,(length(order)+1)], kern='gauss', ev=dat(cv=FALSE))

    #GCV (p 50 in Loader 1999)
    res$GCV[i]<--2*tmp$dp["lk"]*n/(n-tmp$dp["df1"])^2
    if(res$GCV[i]==0) res$GCV[i]<-NA
    #df
    res$GCV.df[i]<-tmp$dp["df1"]

    if(echo == TRUE){
            cat(i, " / ", dim(bogrid)[1], "\n")
    }
    }

    out<-list(grid=res, order=order, deg=deg, step=step, call=deparse(match.call()), cv=cv, x=x)
    class(out) <- "ll.order"
    out
}

##############################################################################################
lpx<-function (x, nn = 0, h = 0, adpen = 0, deg = 2, acri = "none",
    scale = FALSE, style = "none"){
##############################################################################################
#locfit hack to make ll.order work with locfit >1.5
##############################################################################################
    x <- cbind(x)
    z <- as.list(match.call())
    z[[1]] <- z$nn <- z$h <- z$adpen <- z$deg <- z$acri <- z$scale <- z$style <- NULL
    #dimnames(x) <- list(NULL, z)
    if (missing(nn) & missing(h) & missing(adpen))
        nn <- 0.7
    attr(x, "alpha") <- c(nn, h, adpen)
    attr(x, "deg") <- deg
    attr(x, "acri") <- acri
    attr(x, "style") <- style
    attr(x, "scale") <- scale
    class(x) <- "lp"
    x
}

##############################################################################################
"summary.ll.order" <- function(object, GCV=FALSE, ...){
##############################################################################################
#this is the generic summary function for ll.order objects
##############################################################################################
    ans <- as.data.frame(matrix(NA, ncol = 4, nrow = length(object$order)))
    names(ans) <- c("order", "cv.min", "bandwidth.opt", "df")
    if(object$cv==FALSE) GCV = TRUE

    ans$GCV.min<-NA
    ans$GCV.bandwidth.opt <- NA
    ans$GCV.df <- NA

    ans[, 1] <- object$order
    for(i in 1:length(object$order)) {
       if(object$cv==TRUE){
        ans[i, 2] <- min(object$grid$CV[object$grid$order == i])
        wh<- which(object$grid$CV[object$grid$order == i]==ans[i,2])
        ans[i, 3] <- object$grid$bandwidth[object$grid$order == i][wh]
        ans[i, 4] <- object$grid$df[object$grid$order == i][wh]
        }

        ans$GCV.min[i] <- min(object$grid$GCV[object$grid$order == i], na.rm=TRUE)
        wh <- which(object$grid$GCV[object$grid$order == i]==ans$GCV.min[i])
        ans$GCV.bandwidth.opt[i] <- object$grid$bandwidth[object$grid$order == i][wh]
        ans$GCV.df[i] <- object$grid$GCV.df[object$grid$order == i][wh]
    }

    max <- ans$order[order(ans$cv.min)][1]
    cat(paste("The estimated order is ", max, " with a cross-validation error of ", round(ans$
        cv.min[order(ans$cv.min)[1]], 3), "\nand Gaussian bandwidth ", round(as.numeric(
        ans$bandwidth.opt[order(ans$cv.min)][1]), 3), ". (using a local polynomial with ", object$deg,
         " degrees).\n\n", sep = ""))
    ans
}

##############################################################################################
"plot.ll.order" <- function(x, ...){
##############################################################################################
#this is the generic plot function for ll.order objects
##############################################################################################
    object<-x
    ans <- as.data.frame(matrix(NA, ncol = 2, nrow = length(object$order)))
    names(ans) <- c("order", "cv.min")
    order<-object$order
    ans[, 1] <- order
    for(i in 1:length(order)) {
        ans[i, 2] <- min(object$grid$CV[object$grid$order == i])
    }
    plot(cv.min ~ order, data = ans, type = "b")
}

##############################################################################################
"print.ll.order" <- function(x, verbose = FALSE, ...){
##############################################################################################
#this is the generic print function for ll.order objects
#
#ARGUMENTS
#verbose   if FALSE, summary is used. If TRUE, the raw list is echoed
##############################################################################################
    object<-x
    if(!verbose) {
    out <- summary(object)
    print(out)
    cat("\n\nFor a raw listing use print(object, verbose=TRUE)\n")
    }
    if(verbose) {
        print.default(object)
    }
}

##############################################################################################
"predict.ll.order"<- function(object, ...){
##############################################################################################
    require(locfit)

    object=object
    x=object$x
    ans <- as.data.frame(matrix(NA, ncol = 4, nrow = length(object$order)))
    names(ans) <- c("order", "cv.min", "bandwidth.opt", "df")
    ans[, 1] <- object$order
    for (i in 1:length(object$order)) {
        if (object$cv == TRUE) {
            ans[i, 2] <- min(object$grid$CV[object$grid$order ==
                i])
            wh <- which(object$grid$CV[object$grid$order == i] ==
                ans[i, 2])
            ans[i, 3] <- object$grid$bandwidth[object$grid$order ==
                i][wh]
            ans[i, 4] <- object$grid$df[object$grid$order ==
                i][wh]
        }
        if (object$cv == FALSE) {
        ans[i, 2] <- min(object$grid$GCV[object$grid$order ==
            i], na.rm = TRUE)
        wh <- which(object$grid$GCV[object$grid$order == i] ==
            ans$GCV.min[i])
        ans[i, 3] <- object$grid$bandwidth[object$grid$order ==
            i][wh]
        ans[i, 4] <- object$grid$GCV.df[object$grid$order ==
            i][wh]
    }
    }

    ord = ans$order[order(ans$cv.min)][1]
    bw=as.numeric(ans$bandwidth.opt[order(ans$cv.min)][1])
    deg=object$deg
    step=object$step

    res<-data.frame(obs=x, pred=rep(NA, length(x)))

    bogrid<-expand.grid(bw, ord)

    T <- length(x)

    tmu=mean(x[(max(ord) + 1):T])
    tsd=sqrt(var(x[(max(ord) + 1):T]))
    cvseries <- (x - tmu)/tsd
    ldata<-mkx(cvseries, step:(ord+step-1))

    n<-dim(ldata)[1]

    for(k in 1:(T-ord)){
        tmp<-NULL
        tdata=ldata[-k,]
        tmp<-locfit.raw(lpx(tdata[,1:(bogrid[,2])], deg=deg, h=bogrid[,1]),
                y=tdata[,ord+1], kern='gauss', ev=ldata[k,1:(bogrid[,2])])
        res$pred[k+ord]=tsd*predict(tmp)+tmu
        }
      res
}

##############################################################################################
"prediction.profile.ll"<- function(x, step=1:10,order = 1:5, deg = 2, bandwidth = c(seq(0.3, 1.5, by = 0.1), 2:10)){
##############################################################################################
    res<-as.data.frame(matrix(NA, ncol=5, nrow=length(step)))
    names(res)<- c("step", "CV", "order", "bandwidth", "df")
    res$step<-step
    for(k in 1:length(step)){
          tmp<-ll.order(x, order = order, step=step[k], deg = deg, bandwidth = bandwidth, cv=TRUE, echo=FALSE)
          wh<-which(tmp$grid$CV==min(tmp$grid$CV))
          res[k,"CV"]<-tmp$grid$CV[wh]
          res[k,"bandwidth"]<-tmp$grid$bandwidth[wh]
          res[k,"order"]<-tmp$grid$order[wh]
          res[k,"df"]<-tmp$grid$df[wh]
          cat(k, '\n')
          }
    res2<-list(ppll=res)
    class(res2) <- "ppll"
    return(res2)
}

##############################################################################################
"plot.ppll"<-function(x, ...){
##############################################################################################
    object<-x
    plot(object$ppll$step, 1-object$ppll$CV, ylim=c(min(c(0,1-object$ppll$CV)), 1),
    xlab='prediction interval', ylab='predictability', type='b')
}

##############################################################################################
spec.lomb <- function (y=stop("no data arg"), x=stop("no time arg"), freq=NULL){
##############################################################################################
#spec.lomb is a function to estimate the Lomb periodogram for peforming an spectral analysis
#of unevenly sampled data.
#
#The following implentation is solely based on native Splus functions, it should therefore run
#on  all platforms.
#
#The outer code is in public domain, many of the inner functions are part of S-plus (which
#is not in public domain). Upon modifying the code, please comment it. My name should be
#removed as appropriate. I would, of course, be grateful if you notify me about any
#improvements made.
#
#REFERENCE
#Lomb, N.R. 1976, Astrophysics and Space Science 39: 447-462
#
#REQUIRED ARGUMENTS
#y         vector of length n representing the unevenly sampled time series
#x         the a vector (of length n) representing the times of observation
#freq      the frequencies at which the periodogram is to be calculated. If NULL
#          the canonical frequncies are used
#
#VALUE
#an object of class Lomb is returned consisting of the following components:
#freq      the frequencies as supplied
#spec      the estimated amplitudes at the different frequencies
##############################################################################################
  if(is.null(freq)){
    nyear <- max(x)-min(x)+1
    f <- seq(0,.5,length=nyear/2)
  }
  else{
    f <- freq
  }

  # Check arguments
  if (length(y) != length(x)) stop("y and x different lengths");
  if (min(f) < 0 || max(f) > 1) stop("freq must be between 0 and 1");
  if (min(f) == 0 ) f <- f[f>0];    # Get rid of zeros

  nt <- length(x);          # Number of datapoints
  nf <- length(f);          # Number of frequencies
  ones.t <- rep(1,nt);          # Useful unit vectors
  ones.f <- rep(1,nf);

  ## Convert to angular frequencies
  omega <- 2 * pi * f;

  ## Stats of the time series
  hbar <- mean(y);
  hvar <- var(y);
  hdev <- y - hbar;

  ## Calculate the vector of taus
  two.omega.t <- 2 * omega %*% t(x);
  sum.sin <- sin(two.omega.t) %*% ones.t;
  sum.cos <- cos(two.omega.t) %*% ones.t;
  tau <- atan(sum.sin/sum.cos) / (2*omega);

  ## Calculate the trig functions that go into the main expression
  t.m.tau <- (ones.f %*% t(x)) - (tau %*% t(ones.t));
  omega.ttau <- (omega %*% t(ones.t)) * t.m.tau;
  sin.ott <- sin(omega.ttau);
  cos.ott <- cos(omega.ttau);
  z <- ((cos.ott %*% hdev)^2 / ((cos.ott^2) %*% ones.t) +
    (sin.ott %*% hdev)^2 / ((sin.ott^2) %*% ones.t)) / (2 * hvar);

  max <- z == max(z,na.rm=TRUE)
  max <- max[is.na(max)==FALSE]
  P <- 1 - (1-exp(-z[max]))^(length(x))

  res <- list(spec=z[,1], freq=f, f.max=f[max], per.max=1/f[max], p = P)
  class(res) <- "lomb"
  res
}

##############################################################################################
plot.lomb <- function(x, ...){
##############################################################################################
object<-x
plot(object$freq,object$spec, type="l", xlab="frequency", ylab="amplitude")
}

##############################################################################################
summary.lomb <- function(object, ...){
##############################################################################################
list(period=object$per.max,p.val=object$p)
}


##############################################################################################
"add.test"<-function(x, order, n.cond = FALSE){
##############################################################################################
#Chen et al's (1995) Lagrange multiplier test for additivity.
#
#ARGUMENTS
#x        is a vector representing the time series.
#order    a number representing the order to be considered.
#n.cond   is the number of observations to condition on (must be >= order).
##############################################################################################
        require(acepack)
    resid.ace <- function(aceobj){
    aceobj$ty - apply(aceobj$tx, 1, sum)
    }
    if(!n.cond){
        n.cond <- order
        }
    nx <- length(x)
    tmp.mkx <- matrix(0, (nx - n.cond), (n.cond + 1))
    for(i in 1:n.cond)
        tmp.mkx[, i] <- x[(n.cond + 1 - i):(nx - i)]
    tmp.mkx[, (n.cond + 1)] <- x[(n.cond + 1):nx]

    tmp.ace <- ace(tmp.mkx[, 1:order], tmp.mkx[, (n.cond + 1)], lin = 0)
    tmp.resid1 <- resid.ace(tmp.ace)
    h <- 0
    K <- ((order - 1) * order * (order + 7))/6
    tmp.resid2 <- matrix(NA, ncol = K, nrow = dim(tmp.mkx)[1])
    for(i in 1:order)
        for(j in i:order)
            if(i != j) {
                tmp <- apply(tmp.mkx[, c(i, j)], 1, prod)
                h <- h + 1
                tmp.resid2[, h] <- resid.ace(ace(tmp.mkx[, 1:order], tmp, lin = 0))
            }
    for(i in 1:order)
        for(j in i:order)
            for(k in j:order)
                if(i != j | i != k | j != k) {
                  tmp <- apply(tmp.mkx[, c(i, j, k)], 1, prod)
                  h <- h + 1
                  tmp.resid2[, h] <- resid.ace(ace(tmp.mkx[, 1:order], tmp, lin = 0))
                }
    resid.lm <- lm(tmp.resid1 ~ tmp.resid2)
    unlist(list(chisq = round(dim(tmp.mkx)[1] * summary(resid.lm)$r.squared,4), df=K,
               p.val = round(1 - pchisq(dim(tmp.mkx)[1] * summary(resid.lm)$r.squared, K),4)))
}


##############################################################################################
"lin.test" <- function(x, order){
##############################################################################################
#Tsay's (1986) Tukey one-degree-of-freedom test for linearity.
#
#ARGUMENTS
#x        is a vector representing the time series.
#order    a number representing the order to be considered.
##############################################################################################
    nx <- length(x)
    Y <- matrix(0, (nx - order), (order + 1))
    for(i in 1:order)
        Y[, i] <- x[(order + 1 - i):(nx - i)]
    Y[, (order + 1)] <- x[(order + 1):nx]
    D <- dim(Y)
    ur0 <- residuals(switch(as.character(D[2] - 1),
        "1" = lm(Y[, D[2]] ~ Y[, 1]),
        "2" = lm(Y[, D[2]] ~ Y[, 1] + Y[, 2]),
        "3" = lm(Y[, D[2]] ~ Y[, 1] + Y[, 2] + Y[, 3]),
        "4" = lm(Y[, D[2]] ~ Y[, 1] + Y[, 2] + Y[, 3] + Y[, 4]),
        "5" = lm(Y[, D[2]] ~ Y[, 1] + Y[, 2] + Y[, 3] + Y[, 4] + Y[, 5])
        ))
    ur1 <- residuals(switch(as.character(D[2] - 1),
        "1" = lm(ur0 ~ poly(Y[, 1], degree= 2)),
        "2" = lm(ur0 ~ poly(Y[, 1], Y[, 2], degree= 2)),
        "3" = lm(ur0 ~ poly(Y[, 1], Y[, 2], Y[, 3], degree= 2)),
        "4" = lm(ur0 ~ poly(Y[, 1], Y[, 2], Y[, 3], Y[, 4], degree= 2)),
        "5" = lm(ur0 ~ poly(Y[, 1], Y[, 2], Y[, 3], Y[, 4], Y[, 5], degree= 2))
        ))
    m <- switch(as.character(D[2] - 1),
        "1" = 1,
        "2" = 3,
        "3" = 6,
        "4" = 10,
        "5" = 15)
    Fval <- ((sum(ur0^2) - sum(ur1^2))/m)/(sum(ur1^2)/(D[1] - m - 1))
    pval <- 1 - pf(Fval, m, D[1] - m - 1)
    unlist(list(order = order, F = round(Fval, 4), df1 = m, df2 = D[1] - m - 1, p = round(pval, 4)))
}

##############################################################################################
"portman.Q" <- function(x, K){
##############################################################################################
#Ljung-Box test for whiteness of a time series
#
#ARGUMENTS
#x        is a vector representing the time series.
#K        the number of lags in the ACF to be included.
##############################################################################################
    Q <- 0
    n <- length(x)
    p <- acf(x, plot = FALSE, lag.max = K)$acf[2:(K + 1)]
    for(k in 1:K)
        Q <- Q + p[k]^2/(n - k)
    Q <- n * (n + 2) * Q
    res <- list(chisq = round(Q,4), df = K, p.val = round(1 - pchisq(Q, K),4))
    unlist(res)
}

##############################################################################################
"specar.ci" <- function(x, order, resamp = 500, nfreq = 100, echo = TRUE){
##############################################################################################
#specar.ci is a funcion to estimate a confidence interval for the power spectrum and
#in particular a confidence interval for the dominant period.
#The function uses a parametric ar-model to attain the estimate -- using the native
#S-plus spec.ar function. A confidence interval is obtained by resampling the ar-coefficients
#using the variance-covariance matrix from the arima.mle object.
#
#BACKGROUND
#I coded this function to assess whether a virus altered the natural period of a cyclic host
#species:
#Bjornstad et al. (1998; Journal of Animal Ecology (1998) 67:110-126; see
#               http://onb.ent.psu.edu/publ/plodia/plodia.pdf)
#
#The following implentation is solely based on native Splus functions, it should therefore run
#on  all platforms. The drawback is that it may be slow.
#
#I've commented the code crudely to make everything as transparent as possible. I'd greatly
#appreciate any comment or suggestion (onb1@psu.edu). I know the code can be
#optimized greatly, but i'm not a programmer. I'll try to keep an updated set of code at
#http://onb.ent.psu.edu/software.html.
#
#The outer code is in public domain, many of the inner functions are part of S-plus (which
#is not in public domain). Upon modifying the code, please comment it. My name should be
#removed as appropriate. I would, of course, be grateful if you notify me about any
#improvements made.
#
#REQUIRED ARGUMENTS
#x         vector of length n representing the time series. The data will be scaled to
#       unit variance prior to estimation
#order     the order to be used for the ar model. If "aic" is used, the order will
#             be selected automatically using the AIC criterion. Note that if
#             a zero'th order process is chosen by the aic, a first order process
#         will be used.
#resamp      is the number of resamples of the ar-coefficients from the var-covar matrix
#nfreq     is the number of points at which to save the value for the power spectrum (and
#             confidence envelope).
#echo      if True, a counter for each resamp shows the progress.
#
#VALUE
#an object of class specar.ci is returned consisting of the following components:
#order     is the ar-order
#spectrum  $freq is the the spectral frequencies
#          $spec is the estimated power-spectrum of the data (in decibel)
#          $maxfreq is the frequency of maximum power of the data
#resamp    is the summary across the resamples
#resamp$spectrum   gives the quantile summary for the resampling distribution of
#         of the spectral powers
#resamp$maxfreq    gives the full vector of output for the resampled max.frequencies
############################################################################################
    require(stats)
    if(order == "aic") {
        s.ar <- ar(x, aic = TRUE)
        if(s.ar$order == 0) {
            s.ar <- ar.mle(x, order.max = 1, aic = FALSE)
            order <- s.ar$order
        }
        else {
            order <- s.ar$order
        }
    }
    else {
        s.ar <- ar.mle(x, order.max = order, aic = FALSE)
    }
    real <- spec.ar(s.ar, n.freq = nfreq, plot = FALSE)
    trekk <- matrix(NA, ncol = nfreq, nrow = resamp)
    maxfreq <- 1:resamp

    s.ar2<-s.ar

    for(i in 1:resamp) {
        if(order > 1) {
            vs <- svd(s.ar$asy.var.coef)
            vsqrt <- t(vs$v %*% (t(vs$u) * sqrt(vs$d)))
            ans <- matrix(rnorm(order), nrow = 1) %*% vsqrt
            ans <- sweep(ans, 2, s.ar$ar, "+")
        }
        if(order == 1) {
            ans <- rnorm(1, s.ar$ar, sqrt(s.ar$var.coef))
        }
        s.ar2$ar <- as.vector(ans)
        s.ar.mle3 <- spec.ar(s.ar2, n.freq = nfreq, plot = FALSE)
        trekk[i,  ] <- s.ar.mle3$spec
        maxfreq[i] <- s.ar.mle3$freq[match(max(s.ar.mle3$spec), s.ar.mle3$spec)]
        if(echo) {
            cat(i, "\n")
        }
    }
    trekk <- apply(trekk, 2, quantile, probs = c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9,
        0.95, 0.975, 1))
    dimnames(trekk) <- list(c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1),
        NULL)
    res <- list(spectrum = list(freq = real$freq, spec = real$spec, maxfreq = real$freq[match(
        max(real$spec), real$spec)]), order = order, resamp = list(spectrum = trekk,
        maxfreq = maxfreq))
    class(res) <- "specar.ci"
    res
}


##############################################################################################
"summary.specar.ci" <- function(object, period = TRUE, ...){
##############################################################################################
#this is the generic summary function for specar.ci objects
#
#ARGUMENTS
#period    if T, the summary is in terms of period (=1/freq) rather than frequency
##############################################################################################
    if(period == TRUE) {
        ans <- list(period = unlist(list(period = 1/object$spectrum$maxfreq,
            ci.lower = as.vector(1/quantile(object$resamp$maxfreq, probs = c(
            0.975))), ci.upper = as.vector(1/quantile(object$resamp$maxfreq,
            probs = c(0.025))))), order = as.vector(object$order),
            resamp.summary = summary(1/object$resamp$maxfreq))
    }
    if(period == FALSE) {
        ans <- list(frequency = unlist(list(freq = object$spectrum$maxfreq,
            ci.lower = as.vector(quantile(object$resamp$maxfreq, probs = c(
            0.025))), ci.upper = as.vector(quantile(object$resamp$maxfreq,
            probs = c(0.975))))), order = as.vector(object$order),
            resamp.summary = summary(object$resamp$maxfreq))
    }
    ans
}

##############################################################################################
"plot.specar.ci" <- function(x, period = TRUE, ...){
##############################################################################################
#this is the generic plot function for specar.ci objects
#
#ARGUMENTS
#period    if T, the summary is in terms of period (=1/freq) rather than frequency
##############################################################################################
    object<-x
    if(period == TRUE) {
            n <- length(object$spectrum$freq)
            plot((1/object$spectrum$freq)[2:n], (object$spectrum$spec)[
                2:n], ylim = range(object$resamp$spectrum[c(2, 10), 2:n]),
                xlab = "period", ylab = "amplitude", type = "l")
            lines((1/object$spectrum$freq)[2:n], object$resamp$spectrum[
                "0.025", 2:n])
            lines((1/object$spectrum$freq)[2:n], object$resamp$spectrum[
                "0.975", 2:n])
    }
    if(period == FALSE) {
            n <- length(object$spectrum$freq)
            plot(object$spectrum$freq, object$spectrum$spec, ylim =
                range(object$resamp$spectrum[c(2, 10),  ]), xlab =
                "frequency", ylab = "amplitude", type = "l")
            lines(object$spectrum$freq, object$resamp$spectrum["0.025",
                ])
            lines(object$spectrum$freq, object$resamp$spectrum["0.975",
                ])
    }
}

##############################################################################################
mkx<-function(x, lags){
##############################################################################################
# U. Lall and A. Sharma - Lall, U. & Sharma, A. (1996) A nearest neighbor
#bootstrap for time series resampling. Water Resources Research, 32, 679-693.
#
#function to create matrix of lagged time series.
#x is the univariate time series
#lags is the vector of lags. If lags contains 1 and 4 (say) then
#x1 (output) would consist of xt-1, xt-4, xt.
    nx <- length(x)
    nl <- length(lags)
    ml <- max(lags)
    x1 <- matrix(0, (nx - ml), (nl + 1))
    for(i in 1:nl)
        x1[, i] <- x[(ml + 1 - lags[i]):(nx - lags[
            i])]
    x1[, (nl + 1)] <- x[(ml + 1):nx]
    x1
}



