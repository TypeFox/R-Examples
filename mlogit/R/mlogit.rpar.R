make.beta <- function(mua, siga, rpar, random.nb, correlation){
  nr <- names(rpar)
  censored <-    nr[rpar == "cn"]
  lognormal <-   nr[rpar == "ln"]
  truncated <-   nr[rpar == "tn"]
  normal  <-     nr[rpar ==  "n"]
  uniform  <-    nr[rpar ==  "u"]
  triangular  <- nr[rpar ==  "t"]
  Ka <- ncol(random.nb)
  R <- nrow(random.nb) 
  
  betaa <- matrix(NA, R, Ka)
  betaa.mu <- betaa.sigma <- betaa

  colnames(betaa) <- colnames(betaa.mu) <- colnames(betaa.sigma) <-
    colnames(random.nb) <- names(mua)
  if (correlation){
    colnames(random.nb) <- NULL
    CC <- makeC(siga)
    sigeta <- tcrossprod(random.nb, CC)
    colnames(sigeta) <- nr
#    mymua <- mua
#    mysiga <- sigeta
#    betaa <- t(mymua + t(mysiga))
    betaa <- t(mua + t(sigeta))
    betaa.mu <- matrix(1, R, Ka)
    betaa.sigma <- random.nb[, rep(1:Ka, Ka:1)]

    for (i in 1:Ka){
      sigi <- i + cumsum(c(0, (Ka-1):1))[1:i]
      if (rpar[i] == "cn"){
        betaa[, i] <- pmax(betaa[, i], 0)
#        betaa.mu[, i] <- (betaa[, i] > 0) * 1 + 0
#        betaa.sigma[, sigi] <- ( (betaa[, i] > 0) * 1 + 0) * betaa.sigma[, sigi]
        betaa.mu[, i] <- as.numeric(betaa[, i] > 0)
        betaa.sigma[, sigi] <- as.numeric(betaa[, i] > 0) * betaa.sigma[, sigi]
      }
      if (rpar[i] == "ln"){
        betaa[, i] <- exp(betaa[, i])
        betaa.mu[, i] <- betaa[, i]
        betaa.sigma[, sigi] <- betaa[, i] * betaa.sigma[, sigi]
      }
    }
  }
  else{
    if(length(censored) > 0){
      sel <- censored
      betaa[,sel] <- pmax(t(mua[sel]+siga[sel]*t(random.nb[,sel,drop=F])),0)
      betaa.mu[,sel] <- as.numeric(betaa[,sel]>0)
      betaa.sigma[,sel] <- betaa.mu[,sel]*random.nb[,sel]
    }
    if(length(lognormal) > 0){
      sel <- lognormal
      betaa[,sel] <- exp(t(mua[sel]+siga[sel]*t(random.nb[,sel,drop=F])))
      betaa.mu[,sel] <- betaa[,sel,drop=F]
      betaa.sigma[,sel] <- betaa.mu[,sel]*random.nb[,sel]
    }
    if(length(normal) > 0){
      sel <- normal
      betaa[,sel] <- t(mua[sel]+siga[sel]*t(random.nb[,sel,drop=F]))
      betaa.mu[,sel] <- 1
      betaa.sigma[,sel] <- random.nb[,sel,drop=F]
    }

    if(length(uniform) > 0){
      sel <- uniform
      etauni <- pnorm(random.nb[,sel,drop=F])
      betaa[,sel] <- t(mua[sel]-siga[sel]+2*t(etauni)*siga[sel])
      betaa.mu[,sel] <- 1
      betaa.sigma[,sel] <- 2*etauni-1
    }

    if(length(triangular) > 0){
      sel <- triangular
      eta05 <- random.nb[,sel,drop=F]<0.5
      betaa.mu[,sel] <- 1
      betaa.sigma[,sel] <- eta05*(sqrt(2*pnorm(random.nb[,sel,drop=F]))-1)+
        !eta05*(1-sqrt(2*(1-pnorm(random.nb[,sel,drop=F]))))
      
      betaa[,sel] <- t(mua[sel]+siga[sel]*t(betaa.sigma[,sel]))
    }
  }
  list(betaa = betaa, betaa.mu = betaa.mu, betaa.sigma = betaa.sigma)
}

gnrpoints <- function(low, up, n = 100){
  low+(up-low)*(0:n)/n
}

halton <- function(prime = 3, length = 100, drop = 10){
  halt <- 0
  t <- 0
  while(length(halt) < length + drop){
    t <- t + 1
    halt <- c(halt, rep(halt, prime - 1) + rep(seq(1, prime - 1, 1) / prime ^ t, each = length(halt)))
  }
  halt[(drop + 1):(length + drop)]
}

make.random.nb <- function(R, Ka, halton){
# Create the matrix of random numbers
  if (!is.null(halton)){
    length.halton <- rep(R,Ka)
    prime <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
               47, 53, 59, 61, 71, 73, 79, 83, 89, 97, 101, 103,
               107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
               173, 179, 181, 191, 193, 197, 199)
    drop.halton <- rep(100,Ka)
    if (!is.na(halton) && !is.null(halton$prime)){
      if (length(halton$prime) != Ka){
        stop("wrong number of prime numbers indicated")
      }
      else{
        prime <- halton$prime
      }
      if (!is.na(halton) && !is.null(halton$drop)){
        if (!length(halton$drop) %in% c(1,Ka)) stop("wrong number of drop indicated")
        if (length(halton$drop) == 1){
          drop.halton <- rep(halton$drop,Ka)
        }
        else{
          drop.halton <- halton$drop
        }
      }
    }
    random.nb <- numeric(0)
    i <- 0
    for (i in 1:Ka){
      random.nb <- cbind(random.nb,qnorm(halton(prime[i],R,drop.halton[i])))
    }
  }
  else{
    random.nb <- matrix(rnorm(R*Ka), ncol = Ka, nrow = R)
  }
  random.nb
}

makeC <- function(x){
  # create the lower triangular C matrix
  K <- (-1+sqrt(1 + 8*length(x)))/2
  mat <- matrix(0, K, K)
  mat[lower.tri(mat, diag = TRUE)] <- x
  mat
}

make.rpar <- function(rpar, correlation, estimate, norm){
  K <- length(rpar)
  nr <- names(rpar)
  rpar <- lapply(rpar, function(x) list(dist = x))
  if (correlation){
    Ktot <- length(estimate)
    v <- estimate[(Ktot-0.5*K*(K+1)+1):Ktot]
    v <- tcrossprod(makeC(v))
    colnames(v) <- rownames(v) <- nr
    sv <- sqrt(diag(v))
    names(sv) <- nr
  }      
  for (i in (1:K)){
    m <- estimate[nr[i]]
    if (!correlation){
      s <- estimate[paste("sd.", nr[i], sep="")]
    }
    else{
      s <- sv[i]
    }
    names(m) <- names(s) <- NULL
    rpar[[i]]$mean <- m
    rpar[[i]]$sigma <- s
    rpar[[i]]$name <- nr[[i]]
    if (!is.null(norm)){
      vn <- estimate[norm]
      names(vn) <- NULL
      rpar[[i]]$norm <- vn
    }
  }
  z <- lapply(rpar,function(x){attr(x,"class")="rpar" ;x})
  if (correlation) attr(z, 'covariance') <- v
  z
}

plot.rpar <- function(x, norm = NULL, type = c("density", "probability"), ...){
  type <- match.arg(type)
  if (type == "density") f <- drpar
  if (type == "probability") f <- prpar
  marg <- .05
  law <- x$dist
  rg <- rg(x)
  low <- rg[1]
  np <- x$name
  neg.values <- ifelse(low < 0,TRUE,FALSE)
  up <- rg[2]
  if (!is.finite(low)) low <- qrpar(x, norm = norm)(0.005)
  if (!is.finite(up)) up <- qrpar(x, norm = norm)(0.995)
  ptstot <- gnrpoints(low, up, 1000)
  ytot <- do.call(f, list(x = x, norm = norm))(ptstot)
  ymax <- max(ytot)*(1 + marg)
  plot(ptstot, ytot, type = "n", ann = F, xaxs = "i", yaxs = "i",
       las = 1, ylim = c(0, ymax),
       xlim = c(low - marg *(up - low), up + marg * (up - low)))
  ma <- paste("Distribution of", np)
  if (neg.values){
    pourc0 <- prpar(x)(0)
    ma <- paste(ma, ":", round(pourc0 * 100, 0),"% of 0")
    if (type == "density"){
      if (low < 0){
        ptsneg <- gnrpoints(low, 0, 10)
        yneg <- do.call(f, list(x = x, norm = norm))(ptsneg)
        polygon(c(low, ptsneg, 0),c(0, yneg, 0), col = "lightblue", border = NA)
      }
    }
    else{
      segments(low - marg * (up - low), pourc0, 0, pourc0, lty = "dotted")
      segments(0, 0, 0, pourc0, lty = "dotted")
    }
  }
  lines(ptstot, ytot)
  if (law == "u" && type == "density"){
    segments(up, 0, up, drpar(x)(up))
    segments(low, 0, low, drpar(x)(low))
  }
  title(main = ma)
}

plot.mlogit <- function(x, par = NULL, norm = NULL, type = c("density", "probability"), ...){
  if (is.null(x$rpar)) stop("the plot method is only relevant for random parameters")
  if (is.null(par)) par <- names(x$rpar)
  if (!is.null(norm)) norm = abs(coef(x)[norm])
  rpar <- x$rpar[par]
  K <- length(rpar)
  if (K > 1){
    nrow <- 1 + (K > 2) + (K > 6)
    ncol <- 1 + (K > 1) + (K > 4)
    opar <- par(mfrow = c(nrow, ncol))
    for (i in names(rpar)){
      plot(rpar(x, i), norm = norm, type = type, ...)
    }
    par(opar)
  }
  else{
    plot(rpar(x, 1), norm = norm, type = type, ...)
  }
}

#rpar extract one or several random parameters from an mlogit object
rpar <- function(x, par = NULL, norm = NULL, ...){
  if (is.null(par)) par <- names(x$rpar)
  if (length(par) == 1){
    result <- x$rpar[[par]]
    if (!is.null(norm)) result$norm <- abs(coef(x)[norm])
  }    
  else{
    result <- x$rpar[par]
    if (!is.null(norm))
      lapply(result, function(x){x$norm <- abs(coef(x)[norm]); x})
  }
  result
}

print.rpar <- function(x, digits = max(3, getOption("digits") - 2),
                       width = getOption("width"), ...){
  dist <- switch(x$dist,
                 "n"  = "normal",
                 "ln" = "log-normal",
                 "cn" = "censored normal",
                 "t"  = "triangular",
                 "u"  = "uniform"
                 )
  npar1 <- switch(x$dist,
                  "n"  = "mean",
                  "ln" = "meanlog",
                  "cn" = "mean",
                  "t"  = "center",
                  "u"  = "center"
                 )

  npar2 <- switch(x$dist,
                  "n"  = "sd",
                  "ln" = "sdlog",
                  "cn" = "sd",
                  "t"  = "span",
                  "u"  = "span"
                 )
  par1 <- x$mean
  par2 <- x$sigma
  cat(paste(dist, " distribution with parameters ",round(par1,3),
            " (", npar1, ")", " and ", round(par2,3),
            " (",npar2,")", "\n",sep=""))
}

summary.rpar <- function(object, ...){
  norm <- object$norm
  rg <- rg.rpar(object, norm)
  Q1 <- qrpar(object, norm)(0.25)
  M <-  qrpar(object, norm)(0.5)
  Q3 <- qrpar(object, norm)(0.75)
  m <- mean(object, norm)
  c('Min.' = as.numeric(rg[1]), '1st Qu.' = as.numeric(Q1), 'Median' = as.numeric(M),
    'Mean' = as.numeric(m), '3rd Qu.' = as.numeric(Q3), 'Max.' = as.numeric(rg[2]))
}

# [cor, cov].rpar extract the covariance or the correlation matrix of
# the random effects
cor.mlogit <- function(x){
  if (is.null(x$rpar) || is.null(attr(x$rpar, 'covariance')))
    stop('cor.mlogit only relevant for random models with correlation')
  cor.mlogit <- cov.mlogit(x)
  K <- nrow(cor.mlogit)
  sd.mlogit <- sqrt(diag(cor.mlogit))
  for (i in 1:K){
    for (j in 1:K){
      cor.mlogit[i,j] <- cor.mlogit[i,j]/sd.mlogit[i]/sd.mlogit[j]
    }
  }
  cor.mlogit
}
cov.mlogit <- function(x){
  if (is.null(x$rpar) || is.null(attr(x$rpar, 'covariance')))
    stop('cov.mlogit only relevant for random models with correlation')
  attr(x$rpar, 'covariance')
}

# if a normalization coefficient is used, m2norm and s2norm transform
# the two parameters of the distribution
m2norm <- function(m, dist, norm){
  switch(dist,
         "n"  = m / norm,
         "ln" = m - log(norm),
         "t"  = m / norm,
         "cn" = m / norm,
         "u"  = m / norm
         )
}
s2norm <- function(s, dist, norm){
  switch(dist,
         "n"  = s / norm,
         "ln" = s,
         "t"  = s / norm,
         "cn" = s / norm,
         "u"  = s / norm
         )
}

# mean, med, rg and stdev methods for rpar and mlogit objects (sd and
# range are not generic, so create a stdev and a rg generic ; median
# is generic, but without ..., so create a med generic
stdev <- function(x, ...){
  UseMethod("stdev")
}
rg <- function(x, ...){
  UseMethod("rg")
}
med <- function(x, ...){
  UseMethod("med")
}
mean.rpar <- function(x, norm = NULL, ...){
  dist <- x$dist
  m <- x$mean
  s <- abs(x$sigma)
  if (!is.null(norm)){
    s <- s2norm(s, dist, norm)
    m <- m2norm(m, dist, norm)
  }
  switch(dist,
         "n"  = m,
         "ln" = exp(m + 0.5 * s^2),
         "u"  = m,
         "t"  = m,
         "cn" = s * dnorm(- m / s) + m * (1 - pnorm(- m / s))
         )
}
med.rpar <- function(x, norm = NULL, ...){
  dist <- x$dist
  m <- x$mean
  s <- abs(x$sigma)
  if (!is.null(norm)){
    s <- s2norm(s, dist, norm)
    m <- m2norm(m, dist, norm) 
  }
  switch(dist,
         "n"  = m,
         "ln" = exp(m),
         "u"  = m,
         "t"  = m,
         "cn" = max(0, m)
         )
}
stdev.rpar <- function(x, norm = NULL, ...){
  dist <- x$dist
  m <- x$mean
  s <- abs(x$sigma)
  if (!is.null(norm)){
    s <- s2norm(s, dist, norm)
    m <- m2norm(m, dist, norm)
  }
  switch(dist,
         "n"  = s,
         "ln" = sqrt(exp(s^2)-1)*exp(m+0.5*s^2),
         "u"  = s^2/3,
         "t"  = s,
         "cn" = sqrt( s^2*(1-pnorm(-m/s))+m*(s*dnorm(-m/s)+m*(1-pnorm(-m/s)))-
           (s*dnorm(-m/s)+m*(1-pnorm(-m/s)))^2)
         )
}

rg.rpar <- function(x, norm = NULL, ...){
  dist <- x$dist
  m <- x$mean
  s <- abs(x$sigma)
  if (!is.null(norm)){
    s <- s2norm(s, dist, norm)
    m <- m2norm(m, dist, norm)
  }
  result <- switch(dist,
                   "n"  = c(-Inf, +Inf),
                   "ln" = c(0, +Inf),
                   "u"  = c(m - s, m + s),
                   "t"  = c(m - s, m + s),
                   "cn" = c(0, +Inf)
                   )
  names(result) <- c('Min.', 'Max.')
  result
}
mean.mlogit <- function(x, par = NULL, norm = NULL, ...){
  if (!is.null(norm)) norm <- abs(coef(x)[norm])
  if (is.null(par)) par <- names(x$rpar)
  sapply(x$rpar[par], function(x) mean(x, norm = norm, ...))
}
med.mlogit <- function(x, par = NULL, norm = NULL, ...){
  if (!is.null(norm)) norm <- abs(coef(x)[norm])
  if (is.null(par)) par <- names(x$rpar)
  sapply(x$rpar[par], function(x) med(x, norm = norm, ...))
}
stdev.mlogit <- function(x, par = NULL, norm = NULL, ...){
  if (!is.null(norm)) norm <- abs(coef(x)[norm])
  if (is.null(par)) par <- names(x$rpar)
  sapply(x$rpar[par], function(x) stdev(x, norm = norm, ...))
}
rg.mlogit <- function(x, par = NULL, norm = NULL, ...){
  if (!is.null(norm)) norm <- abs(coef(x)[norm])
  if (is.null(par)) par <- names(x$rpar)
  if (length(par) > 1)
    sapply(x$rpar[par], function(x) rg(x, norm = norm, ...))
  else rg(x$rpar[[par]])
}

# [qpd]rpar methods for rpar and rlogit objects

qrpar <- function(x, ...){
  UseMethod("qrpar")
}
prpar <- function(x, ...){
  UseMethod("prpar")
}
drpar <- function(x, ...){
  UseMethod("drpar")
}
qrpar.rpar <- function(x, norm = NULL, ...){
  dist <- x$dist
  m <- x$mean
  s <- abs(x$sigma)
  if (!is.null(norm)){
    s <- s2norm(s, dist, norm)
    m <- m2norm(m, dist, norm)
  }
  switch(dist,
         "n"  = function(x=(1:9)/10) qnorm(x, m, s),
         "ln" = function(x=(1:9)/10) qlnorm(x, m, s),
         "u"  = function(x=(1:9)/10) qunif(x, m - s, m + s),
         "t"  = function(x=(1:9)/10) (m - s + sqrt(2 * s^2 * x)) *
         (x <= 0.5)+(m + s - sqrt(2 * s^2 *(1 - x))) *(x > 0.5),
         "cn" = function(x=(1:9)/10) max(0, qnorm(x, m, s))
         )
}
prpar.rpar <- function(x, norm = NULL, ...){
  dist <- x$dist
  m <- x$mean
  s <- abs(x$sigma)
  if (!is.null(norm)){
    s <- s2norm(s, dist, norm)
    m <- m2norm(m, dist, norm)
  }
  switch(dist,
         "n"  = function(x) pnorm(x, m, s),
         "ln" = function(x) plnorm(x, m, s),
         "u"  = function(x) punif(x, m - s, m + s),
         "t"  = function(x) (x >= (m - s) & x < m)*(x - m + s)^2 / (2 * s^2) +
         (x >= m & x <= (m + s)) * (1 - (m + s - x)^2 / (2 * s^2)) + (x > (m + s)) * 1 + 0,
         "cn" = function(x) pnorm(x, m, s)
         )
}
drpar.rpar <- function(x, norm = NULL, ...){
  dist <- x$dist
  m <- x$mean
  s <- abs(x$sigma)
  if (!is.null(norm)){
    s <- s2norm(s, dist, norm)
    m <- m2norm(m, dist, norm)
  }
  switch(dist,
         "n"=function(x) dnorm(x, m, s),
         "ln"=function(x) dlnorm(x, m, s),
         "u"=function(x) (1/s+x*0)*(x >= m-s & x <= m+s) + 0,
         "t"=function(x) (x >= (m-s) & x < m)*(x-m+s)/s^2+
         (x>=m & x <= (m+s))*(s+m-x)/s^2+0,
         "cn"=function(x) dnorm(x,m,s)
         )
}


qrpar.mlogit <- function(x, par = 1, y = NULL, norm = NULL, ...){
  if (is.null(rpar))
    stop("qrpar function only relevant for random parameters models")
  if (length(par) > 1)
    stop("only one parameter should be selected")
  x <- x$rpar[[par]]
  if (is.null(norm)) norm <- abs(coef(x)[norm])
  if (is.null(y)){
    qrpar(x, norm = norm, ...)
  }
  else{
    qrpar(x, norm = norm, ...)(y)
  }
}


prpar.mlogit <- function(x, par = 1, y = NULL, norm = NULL, ...){
  if (is.null(rpar))
    stop("prpar function only relevant for random parameters models")
  if (length(par) > 1)
    stop("only one parameter should be selected")
  x <- x$rpar[[par]]
  if (!is.null(norm)) norm <- coef(x)[norm]
  if (is.null(y)){
    prpar(x, norm = norm, ...)
  }
  else{
    prpar(x, norm = norm, ...)(y)
  }
}

drpar.mlogit <- function(x, par = 1, y = NULL, norm = NULL, ...){
  if (is.null(rpar))
    stop("drpar function only relevant for random parameters models")
  if (length(par) > 1)
    stop("only one parameter should be selected")
  x <- x$rpar[[par]]
  if (!is.null(norm)) norm <- coef(x)[norm]
  if (is.null(y)){
    drpar(x, norm = norm, ...)
  }
  else{
    drpar(x, norm = norm, ...)(y)
  }
}

