
"sm.ancova" <- function(x, y, group, h, model = "none",
         h.alpha =  NA, weights = NA, covar = diag(1/weights), ...) {

x.name <- deparse(substitute(x))
y.name <- deparse(substitute(y))

data   <- sm.check.data(x, y, weights, group, ...)
x      <- data$x
y      <- data$y
weights<- data$weights
group  <- data$group
nobs   <- data$nobs
ndim   <- data$ndim
opt    <- data$options

if(missing(h))
   h <- h.select(x, y, weights = weights, group = group, ...)
else 
   {if(length(h)!=ndim) stop("length(h) does not match size of x")}

covar.set <- FALSE
if (!missing(covar)) {
  if (!is.na(opt$nbins) & opt$nbins!=0)
      stop("if covar is set, nbins must be 0 or NA")
  if (!all(weights == as.integer(rep(1,length(y))))) 
      stop("if covar is set, weights must not be set")
  covar.set <- TRUE
  }

if (missing(weights) & missing(covar))
  replace.na(opt, nbins, round((nobs > 500) * 8 * log(nobs) / ndim))
replace.na(opt, band, TRUE)
if (model == "none") opt$band <- FALSE

if (ndim == 1) {
  if (is.na(h.alpha)) h.alpha <- 2 * diff(range(x)) / nobs
  replace.na(opt, display, "line")
  replace.na(opt, ngrid, 50)
  replace.na(opt, xlab, x.name)
  replace.na(opt, ylab, y.name)
  }
else {
  opt$display <- "none"
  }

fact <- factor(group)
if (ndim==1) {
   ord     <- order(fact, x)
   xx      <- x[ord]
   }
else {
   ord <- order(fact)
   xx  <- x[ord,]
   }
yy         <- y[ord]
weights    <- weights[ord]
fact       <- fact[ord]
fac.levels <- levels(fact)
nlev       <- length(fac.levels)

rawdata <- list(x = xx, y = yy, fac = fact, nbins = opt$nbins, nobs = nobs, 
                     ndim = ndim, devs = 0)
if ((!is.na(opt$nbins)) & (opt$nbins>0)) {
  for (i in 1:nlev) {
    ind            <- (fact==fac.levels[i])
    if (ndim==1)  {xx.ind <- xx[ind]}  else {xx.ind <- xx[ind,]}
    bins <- binning(xx.ind, yy[ind], nbins=opt$nbins)
    if (i ==1 ) {
      x            <- matrix(as.vector(bins$x), ncol=ndim)
      y            <- bins$means
      fac          <- rep(fac.levels[1], length(bins$means))
      weights      <- bins$x.freq
      rawdata$devs <- bins$devs
      }
    else {
      x            <- rbind(x, matrix(as.vector(bins$x), ncol=ndim))
      y            <- c(y, bins$means)
      fac          <- c(fac, rep(fac.levels[i], length(bins$means)))
      weights      <- c(weights, bins$x.freq)
      rawdata$devs <- c(rawdata$devs, bins$devs)
      }
    }
  if (ndim == 1) x <- as.vector(x)
  weights <- as.integer(weights)
  fac     <- factor(fac)
  covar   <- diag(1/weights)
  }
else {
  x   <- xx
  y   <- yy
  fac <- fact
  }
n <- table(fac)


#--------------------------Model testing------------------------

B  <- diag(0, sum(n))
Sd <- diag(0, sum(n))
istart <- 1
for (i in 1:nlev) {
  irange <- istart:(istart + n[i] - 1)
  wi     <- weights[irange]
  if (ndim==1) {
     xi <- x[irange]
     Sd[irange, irange] <-  sm.weight(xi, xi, h, weights=wi, options=opt)
     }
  else {
     xi  <- x[irange,]
     Sd[irange, irange] <- sm.weight2(xi, xi, h, weights=wi, options=opt)
     }
  B[irange, irange]  <- sm.sigweight(xi, weights=wi)
  istart <- istart + n[i]
  }
if (ndim==1) {
   Ss <- sm.weight(x, x, h, weights=weights, options=opt)
   }
else {
   Ss <- sm.weight2(x, x, h, weights=weights, options=opt)
   }
sigma <- sqrt((y %*% B %*% y)[1, 1] + sum(rawdata$devs))

if (model == "equal") {
  Q <- Sd - Ss
  Q <- t(Q) %*% diag(weights) %*% Q
  obs <- ((y %*% Q %*% y) / sigma^2)[1,1]
  }

if (model == "parallel") {
  D <- matrix(0, ncol = nlev - 1, nrow = sum(n))
  istart <- n[1] + 1
  for (i in 2:nlev) {
        D[istart:(istart + n[i] - 1),i - 1] <- 1
    }
  if (ndim==1) {
     Q <- diag(sum(n)) - sm.weight(x, x, h.alpha, weights=weights, options=opt)
     }
  else {
     Q <- diag(sum(n)) - sm.weight2(x, x, h, weights=weights, options=opt)
     }
  Q <- solve(t(D) %*% t(Q) %*% diag(weights) %*% Q %*% D) %*%
                t(D) %*% t(Q) %*% diag(weights) %*% Q
  alpha <- as.vector(Q %*% y)
  ghat  <- as.vector(Ss %*% (diag(sum(n)) - D %*% Q) %*% y)
  ghati <- as.vector(Sd %*% y)
  obs   <- sum(weights*(as.vector(D %*% alpha) + ghat - ghati)^2) / sigma^2
  Q     <- D %*% Q + Ss %*% (diag(sum(n)) - D %*% Q) - Sd
  Q     <- t(Q) %*% diag(weights) %*% Q
  }

p <- NULL
if (!(model == "none")) {
  if (!covar.set) {
      p   <- p.quad.moment(Q - B * obs, covar, obs,
                         sum(weights)-length(weights))
      }
  else {
     p   <- p.quad.moment.old(Q, covar, obs * sigma^2)
     }
  if (model == "equal")    model.name <- "equality"
  if (model == "parallel") model.name <- "parallelism"
  if (opt$verbose > 0) cat("Test of", model.name, ":  h = ",
             signif(h), "   p-value = ", round(p, 4), "\n")
  }

if (ndim == 1)
  sigma <- sigma / sqrt(nobs - 2 * nlev)
else
  sigma <- sigma / sqrt(nobs)

#--------------------------Graphical display------------------------

if (!(opt$display %in% "none")) {

  replace.na(opt, xlim, range(rawdata$x))
  replace.na(opt, ylim, range(rawdata$y))
  if (length(opt$lty) < nlev) opt$lty <- 1:nlev
  if (length(opt$col) < nlev) opt$col <- 2:(nlev + 1)

  plot(rawdata$x, rawdata$y, type = "n",
       xlab = opt$xlab, ylab = opt$ylab, xlim = opt$xlim, ylim = opt$ylim)
  for (i in 1:nlev)
     text(rawdata$x[fac == fac.levels[i]], rawdata$y[fac == fac.levels[i]],
          as.character(fac.levels[i]), col = opt$col[i])

  if (opt$band & nlev > 2) {
     if (opt$verbose > 0) cat("Band available only to compare two groups.\n")
     opt$band <- FALSE
     }
  if (opt$band & covar.set) {
     if (opt$verbose > 0) cat("Band not available when covariance is set.\n")
     opt$band <- FALSE
     }       
      
  if (!opt$band) {
    for (i in 1:nlev) {
      ind <- (fac == fac.levels[i])
      sm.regression(x[ind], y[ind], h = h, weights = weights[ind],
          ngrid = opt$ngrid, add = TRUE, lty = opt$lty[i], col = opt$col[i])
      }
    }
  else {
    eval.points <- opt$eval.points
    if (any(is.na(eval.points))) {
      start.eval <- max(tapply(x, fac, min))
      stop.eval  <- min(tapply(x, fac, max))
      eval.points <- seq(start.eval, stop.eval, length = opt$ngrid)
      }

    ind <- (fac == fac.levels[1])
    model1 <- sm.regression(x[ind], y[ind], h = h,
           eval.points = eval.points, weights = weights[ind],
           options = opt, display = "none", ngrid = opt$ngrid, 
           add = TRUE, lty = 1)
    ind <- fac == fac.levels[2]
    model2 <- sm.regression(x[ind], y[ind], h = h,
           eval.points = eval.points, weights = weights[ind],
           options = opt, display = "none", ngrid = opt$ngrid, 
           add = TRUE, lty = 2)
    model.y <- (model1$estimate + model2$estimate) / 2
    if (model == "parallel")
         model.y <- cbind(model.y - alpha/2, model.y + alpha/2)
    se <- sqrt((model1$se/model1$sigma)^2 + (model2$se/model2$sigma)^2)
    se <- se * sigma
    upper <- model.y + se
    lower <- model.y - se
    if (model == "equal") {
      upper <- pmin(pmax(upper, par()$usr[3]), par()$usr[4])
      lower <- pmin(pmax(lower, par()$usr[3]), par()$usr[4])
      polygon(c(eval.points, rev(eval.points)), c(lower, rev(upper)),
              border = FALSE, col = 5)
      }
    else if (model == "parallel") {
      upper[,1] <- pmin(pmax(upper[,1], par()$usr[3]), par()$usr[4])
      lower[,1] <- pmin(pmax(lower[,1], par()$usr[3]), par()$usr[4])
      upper[,2] <- pmin(pmax(upper[,2], par()$usr[3]), par()$usr[4])
      lower[,2] <- pmin(pmax(lower[,2], par()$usr[3]), par()$usr[4])
      polygon(c(eval.points, rev(eval.points)),
            c(lower[,1],   rev(upper[,1])),
            density = 20, angle = 90, border = FALSE, col = 5)
      polygon(c(eval.points, rev(eval.points)),
            c(lower[,2],   rev(upper[,2])),
            density = 20, angle =  0, border = FALSE, col = 6)
      }
    for (i in 1:nlev)
      text(rawdata$x[fac == fac.levels[i]], rawdata$y[fac == fac.levels[i]],
           as.character(fac.levels[i]), col = opt$col[i])
    lines(eval.points, model1$estimate, lty = opt$lty[1], col = opt$col[1])
    lines(eval.points, model2$estimate, lty = opt$lty[2], col = opt$col[2])
    }
  }

#-------------------------------Output-----------------------------

r <- list(p = p, model = model, sigma = sigma)
if (model == "parallel")
   r <- list(p = p, model = model, sigma = sigma, alphahat = alpha)
if (!(opt$display == "none") & opt$band) {
   r$upper <- upper
   r$lower <- lower
   r$eval.points <- eval.points
   }
r$data <- list(x=x, y=y, group=fac, nbins=rawdata$nbins, devs=rawdata$devs, 
          weights=weights)
r$call <- match.call()
invisible(r)

}
