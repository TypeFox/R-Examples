efp <- function(formula, data = list(),
                type = c("Rec-CUSUM", "OLS-CUSUM", "Rec-MOSUM", "OLS-MOSUM",
                "RE", "ME", "Score-CUSUM", "Score-MOSUM", "fluctuation"),
                h = 0.15, dynamic = FALSE, rescale = TRUE)
{
    if(!inherits(formula, "formula")) {
      X <- if(is.matrix(formula$x))
             formula$x
           else model.matrix(terms(formula), model.frame(formula))
      y <- if(is.vector(formula$y))
             formula$y
           else model.response(model.frame(formula))
    } else {
      mf <- model.frame(formula, data = data)
      y <- model.response(mf)
      X <- model.matrix(formula, data = data)
    }  
   
    n <- nrow(X)
    if(dynamic)
    {
      Xnames <- colnames(X)
      X <- cbind(y[1:(n-1)],X[2:n,])
      colnames(X) <- c("lag", Xnames)
      y <- y[-1]
      n <- n-1
    }
    k <- ncol(X)
    type <- match.arg(type)
    if(type == "fluctuation") type <- "RE"

    retval <- list(process = NULL,
                   type = type,
                   nreg = k,
                   nobs = n,
                   call = match.call(),
                   formula = formula,
                   par = NULL,
                   type.name = NULL,
                   lim.process = NULL,
                   coef = NULL,
                   Q12 = NULL,
                   datatsp = NULL,
		   rescale = rescale)

    orig.y <- NULL

    switch(type,

           ## empirical process of Recursive CUSUM model

           "Rec-CUSUM" = {
               w <- recresid(X, y)
               sigma <- sqrt(var(w))
               process <- cumsum(c(0,w))/(sigma*sqrt(n-k))
               if(is.ts(data)) {
                   if(NROW(data) == n) process <- ts(process, end = end(data), frequency = frequency(data))
               } else {
	         env <- environment(formula)
                 if(missing(data)) data <- env
                 orig.y <- eval(attr(terms(formula), "variables")[[2]], data, env)
                 if(is.ts(orig.y) && (NROW(orig.y) == n))
                   process <- ts(process, end = end(orig.y),
                                 frequency = frequency(orig.y))
               }
               retval$type.name <- "Recursive CUSUM test"
               retval$lim.process <- "Brownian motion"
           },

           ## empirical process of OLS-based CUSUM model

           "OLS-CUSUM" = {
               fm <- lm.fit(X,y)
               e <- fm$residuals
               sigma <- sqrt(sum(e^2)/fm$df.residual)
               process <- cumsum(c(0,e))/(sigma*sqrt(n))
               if(is.ts(data)) {
                   if(NROW(data) == n) process <- ts(process, end = end(data), frequency = frequency(data))
               } else {
	         env <- environment(formula)
                 if(missing(data)) data <- env
                 orig.y <- eval(attr(terms(formula), "variables")[[2]], data, env)
                 if(is.ts(orig.y) && (NROW(orig.y) == n))
                   process <- ts(process, end = end(orig.y),
                                 frequency = frequency(orig.y))
               }
               retval$type.name <- "OLS-based CUSUM test"
               retval$lim.process <- "Brownian bridge"
           },

           ## empirical process of Recursive MOSUM model

           "Rec-MOSUM" = {
               w <- recresid(X, y)
               nw <- n - k
               nh <- floor(nw*h)
               process <- rep(0, (nw-nh))
               for(i in 0:(nw-nh))
               {
                   process[i+1] <- sum(w[(i+1):(i+nh)])
               }
               sigma <- sqrt(var(w)*(nw-1)/(nw-k))
               process <- process/(sigma*sqrt(nw))
               if(is.ts(data)) {
                   if(NROW(data) == n) process <- ts(process, end = time(data)[(n-floor(0.5 + nh/2))], frequency = frequency(data))
               } else {
	         env <- environment(formula)
                 if(missing(data)) data <- env
                 orig.y <- eval(attr(terms(formula), "variables")[[2]], data, env)
                 if(is.ts(orig.y) && (NROW(orig.y) == n)) {
                   process <- ts(process, end = time(orig.y)[(n-floor(0.5 + nh/2))],
                                 frequency = frequency(orig.y))
                 } else {
                   process <- ts(process, end = (n-floor(0.5 + nh/2))/n,
                                 frequency = n)
                 }
	       }
               retval$par <- h
               retval$type.name <- "Recursive MOSUM test"
               retval$lim.process <- "Brownian motion increments"
           },

           ## empirical process of OLS-based MOSUM model

           "OLS-MOSUM" = {
               fm <- lm.fit(X,y)
               e <- fm$residuals
               sigma <- sqrt(sum(e^2)/fm$df.residual)
               nh <- floor(n*h)
               process <- cumsum(c(0,e))
               process <- process[-(1:nh)] - process[1:(n-nh+1)]
	       process <- process/(sigma*sqrt(n))
               if(is.ts(data)) {
                   if(NROW(data) == n) process <- ts(process, end = time(data)[(n-floor(0.5 + nh/2))], frequency = frequency(data))
               } else {
	         env <- environment(formula)
                 if(missing(data)) data <- env
                 orig.y <- eval(attr(terms(formula), "variables")[[2]], data, env)
                 if(is.ts(orig.y) && (NROW(orig.y) == n)) {
                   process <- ts(process, end = time(orig.y)[(n-floor(0.5 + nh/2))],
                                 frequency = frequency(orig.y))
                 } else {
                   process <- ts(process, end = (n-floor(0.5 + nh/2))/n,
                                 frequency = n)
                 }
	       }
               retval$par <- h
               retval$type.name <- "OLS-based MOSUM test"
               retval$lim.process <- "Brownian bridge increments"
           },

           ## empirical process of recursive estimates fluctuation

           "RE" = {
               m.fit <- lm.fit(X,y)
               beta.hat <- m.fit$coefficients
               sigma <- sqrt(sum(m.fit$residual^2)/m.fit$df.residual)
               process <- matrix(rep(0,k*(n-k+1)), nrow=k)
               Q12 <- root.matrix(crossprod(X))/sqrt(n)
               if(rescale)
               {
                   for(i in k:(n-1))
                   {
                       Qi12 <- root.matrix(crossprod(X[1:i,]))/sqrt(i)
                       process[,(i-k+1)] <- Qi12 %*%
                               (lm.fit(as.matrix(X[1:i,]), y[1:i])$coefficients - beta.hat)
                   }
               }
               else
               {
                   for(i in k:(n-1))
                   {
                       process[,(i-k+1)] <- Q12 %*% (lm.fit(as.matrix(X[1:i,]),
                                                           y[1:i])$coefficients - beta.hat)
                   }
               }
               process <- t(cbind(0, process))*matrix(rep((k-1):n,k),
                                                      ncol=k)/(sigma*sqrt(n))
               colnames(process) <- colnames(X)
               if(is.ts(data)) {
                   if(NROW(data) == n) process <- ts(process, end = end(data), frequency = frequency(data))
               } else {
	         env <- environment(formula)
                 if(missing(data)) data <- env
                 orig.y <- eval(attr(terms(formula), "variables")[[2]], data, env)
                 if(is.ts(orig.y) && (NROW(orig.y) == n)) {
                   process <- ts(process, end = end(orig.y),
                                 frequency = frequency(orig.y))
                 }
	       }
               retval$Q12 <- Q12
               retval$type.name <- "RE test (recursive estimates test)"
               retval$lim.process <- "Brownian bridge"
           },

           ## empirical process of moving estimates fluctuation

           "ME" = {
               m.fit <- lm.fit(X,y)
               beta.hat <- m.fit$coefficients
               sigma <- sqrt(sum(m.fit$residual^2)/m.fit$df.residual)
               nh <- floor(n*h)
               process <- matrix(rep(0,k*(n-nh+1)), nrow=k)
               Q12 <- root.matrix(crossprod(X))/sqrt(n)
               if(rescale)
               {
                   for(i in 0:(n-nh))
                   {
                       Qnh12 <- root.matrix(crossprod(X[(i+1):(i+nh),]))/sqrt(nh)
                       process[, i+1] <-  Qnh12 %*% (lm.fit(
                                             as.matrix(X[(i+1):(i+nh),]), y[(i+1):(i+nh)])$coefficients - beta.hat)
                   }
               }
               else
               {
                   for(i in 0:(n-nh))
                   {
                       process[, i+1] <- Q12 %*% (lm.fit(as.matrix(X[(i+1):(i+nh),]),
		         y[(i+1):(i+nh)])$coefficients - beta.hat)
                   }
               }
               process <- nh*t(process)/(sqrt(n)*sigma)
               colnames(process) <- colnames(X)
               if(is.ts(data)) {
                   if(NROW(data) == n) process <- ts(process, end = time(data)[(n-floor(0.5 + nh/2))], frequency = frequency(data))
               } else {
	         env <- environment(formula)
                 if(missing(data)) data <- env
                 orig.y <- eval(attr(terms(formula), "variables")[[2]], data, env)
                 if(is.ts(orig.y) && (NROW(orig.y) == n)) {
                   process <- ts(process, end = time(orig.y)[(n-floor(0.5 + nh/2))],
                                 frequency = frequency(orig.y))
                 } else {
                   process <- ts(process, end = (n-floor(0.5 + nh/2))/n,
                                 frequency = n)
                 }
	       }
               retval$par <- h
               retval$Q12 <- Q12
               retval$type.name <- "ME test (moving estimates test)"
               retval$lim.process <- "Brownian bridge increments"
           },

           "Score-CUSUM" = {
               fm <- lm.fit(X,y)
               e <- as.vector(fm$residuals)
               sigma2 <- sum(e^2)/n
	       k <- k + 1
               ## Q12 <- sqrt(sigma2) * root.matrix(crossprod(X))/sqrt(n)
               ## Q12 <- rbind(cbind(Q12, 0), 0)
               ## Q12[k,k] <- sqrt(2)*sigma2

	       process <- cbind(X * e, e^2 - sigma2)/sqrt(n)
	       Q12 <- root.matrix(crossprod(process))
	       process <- rbind(0, process)
               process <- apply(process, 2, cumsum)
               process <- t(chol2inv(chol(Q12)) %*% t(process))

	       colnames(process) <- c(names(coef(fm)), "(Variance)")
               if(is.ts(data)) {
                   if(NROW(data) == n) process <- ts(process, end = end(data), frequency = frequency(data))
               } else {
	         env <- environment(formula)
                 if(missing(data)) data <- env
                 orig.y <- eval(attr(terms(formula), "variables")[[2]], data, env)
                 if(is.ts(orig.y) && (NROW(orig.y) == n)) {
                   process <- ts(process, end = end(orig.y),
                                 frequency = frequency(orig.y))
                 }
	       }
               retval$type.name <- "Score-based CUSUM test"
               retval$lim.process <- "Brownian bridge"
               retval$Q12 <- Q12
           },

           "Score-MOSUM" = {
               fm <- lm.fit(X,y)
               e <- as.vector(fm$residuals)
               sigma2 <- sum(e^2)/n
	       k <- k + 1
	       nh <- floor(n*h)
               ## Q12 <- sqrt(sigma2) * root.matrix(crossprod(X))/sqrt(n)
               ## Q12 <- rbind(cbind(Q12, 0), 0)
               ## Q12[k,k] <- sqrt(2)*sigma2

	       process <- cbind(X * e, e^2 - sigma2)/sqrt(n)
	       Q12 <- root.matrix(crossprod(process))
	       process <- rbind(0, process)
               process <- apply(process, 2, cumsum)
               process <- process[-(1:nh),] - process[1:(n-nh+1),]
	       process <- t(chol2inv(chol(Q12)) %*% t(process))

	       colnames(process) <- c(names(coef(fm)), "(Variance)")
               if(is.ts(data)) {
                   if(NROW(data) == n) process <- ts(process, end = time(data)[(n-floor(0.5 + nh/2))], frequency = frequency(data))
               } else {
	         env <- environment(formula)
                 if(missing(data)) data <- env
                 orig.y <- eval(attr(terms(formula), "variables")[[2]], data, env)
                 if(is.ts(orig.y) && (NROW(orig.y) == n)) {
                   process <- ts(process, end = time(orig.y)[(n-floor(0.5 + nh/2))],
                                 frequency = frequency(orig.y))
                 } else {
                   process <- ts(process, end = (n-floor(0.5 + nh/2))/n,
                                 frequency = n)
                 }
	       }
               retval$par <- h
               retval$type.name <- "Score-based MOSUM test"
               retval$lim.process <- "Brownian bridge increments"
               retval$Q12 <- Q12
           })


    if(!is.ts(process))
        process <- ts(process, start = 0, frequency = (NROW(process)-1))

    retval$process <- process

    if(is.ts(data) & NROW(data) == n)
        retval$datatsp <- tsp(data)
    else if(!is.null(orig.y) && is.ts(orig.y) & NROW(orig.y) == n)
        retval$datatsp <- tsp(orig.y)
    else
        retval$datatsp <- c(0, 1, n)

    m.fit <- lm.fit(X,y)
    retval$coefficients <- coefficients(m.fit)
    retval$sigma <-  sqrt(sum(m.fit$residual^2)/m.fit$df.residual)
    class(retval) <- c("efp")
    return(retval)
}


plot.efp <- function(x, alpha = 0.05, alt.boundary = FALSE, boundary = TRUE,
                     functional = "max", main = NULL,  ylim = NULL,
                     ylab = "Empirical fluctuation process", ...)
{
    if(is.null(functional)) fun <- "max"
      else fun <- match.arg(functional, c("max", "range", "maxL2", "meanL2"))
    bound <- boundary(x, alpha = alpha, alt.boundary = alt.boundary, functional = fun)
    pos <- FALSE
    ave <- FALSE

    if(is.null(main)){
        if(alt.boundary & fun == "max" & (x$lim.process %in% c("Brownian motion", "Brownian bridge"))){
            main <- paste(x$type.name, "with alternative boundaries")
        }
        else {
            if(alt.boundary) warning("no alternative boundaries available")
            if(fun == "meanL2")
              main <- paste(x$type.name, "with mean L2 norm")
	    else if(fun == "maxL2")
	      main <- paste(x$type.name, "with max L2 norm")
	    else
	      main <- x$type.name
        }
    }

    if(!is.null(functional) && x$lim.process %in% c("Brownian bridge", "Brownian bridge increments")) {
      z <- as.matrix(x$process)
      k <- ncol(z)

      switch(functional,
        "max" = {
          if(k > 1) {
            z <- apply(abs(z), 1, max)
            pos <- TRUE
          }
        },
        "range" = { stop("no plot available for range functional") },
        "maxL2" = {
	  if(x$lim.process == "Brownian bridge") {
            z <- rowSums(z^2)
	    pos <- TRUE
	  } else {
	    stop("no test/plot available for mean L2 functional")
	  }
        },

        "meanL2" = {
	  if(x$lim.process == "Brownian bridge") {
            z <- rowSums(z^2)
	    ave <- TRUE
	    pos <- TRUE
	  } else {
	    stop("no test/plot available for mean L2 functional")
	  }
        })
      z <- ts(as.vector(z), start = start(x$process), frequency = frequency(x$process))
    } else {
      z <- x$process
    }

    if(is.null(ylim)) {
      ymax <- max(c(z, bound))
      if(pos) ymin <- 0
      else ymin <- min(c(z, -bound))
      ylim <- c(ymin, ymax)
    }

    if(boundary)
        panel <- function(y, ...) {
            lines(y, ...)
            lines(bound, col=2)
            lines(-bound, col=2)
            abline(0,0)
        }
    else
        panel <- function(y, ...) {
            lines(y, ...)
            abline(0,0)
        }

    if(any(attr(z, "class") == "mts"))
        plot(z, main = main, panel = panel, ...)
    else {
        plot(z, main = main, ylab = ylab, ylim = ylim, ...)
        if(boundary) {
            lines(bound, col=2)
            if(!pos) lines(-bound, col=2)
            if(ave) {
              avez <- ts(rep(mean(z), length(bound)), start = start(bound), frequency = frequency(bound))
              lines(avez, lty = 2)
            }
        }
        abline(0,0)
    }
}

pvalue.efp <- function(x, lim.process, alt.boundary, functional = "max", h = NULL, k = NULL)
{
  lim.process <- match.arg(lim.process,
    c("Brownian motion", "Brownian bridge", "Brownian motion increments", "Brownian bridge increments"))
  functional <- match.arg(functional, c("max", "range", "meanL2", "maxL2"))

  switch(lim.process,

  "Brownian motion" = {
    if(functional == "max") {
      if(alt.boundary)
      {
        pval <- c(1, 0.997, 0.99, 0.975, 0.949, 0.912, 0.864, 0.806, 0.739, 0.666,
             0.589, 0.512, 0.437, 0.368, 0.307, 0.253, 0.205, 0.163, 0.129, 0.100,
             0.077, 0.058, 0.043, 0.032, 0.024, 0.018, 0.012, 0.009, 0.006, 0.004,
             0.003, 0.002, 0.001, 0.001, 0.001)
        critval <- (10:44)/10
        p <- approx(critval, pval, x, rule=2)$y
      } else {
        p <- ifelse(x < 0.3, 1 - 0.1465*x,
        2*(1-pnorm(3*x) + exp(-4*(x^2))*(pnorm(x)+pnorm(5*x)-1)-exp(-16*(x^2))*(1-pnorm(x))))
      }
    } else {
      stop("only max functional implemented for Brownian motion")
    }
    p <- 1 - (1-p)^k
  },

  "Brownian bridge" = {
    switch(functional,

    "max" = {
      if(alt.boundary)
      {
        pval <- c(1, 1, 0.997, 0.99, 0.977, 0.954, 0.919, 0.871, 0.812, 0.743,
          0.666, 0.585, 0.504, 0.426, 0.353, 0.288, 0.230, 0.182, 0.142, 0.109,
          0.082, 0.062, 0.046, 0.034, 0.025, 0.017, 0.011, 0.008, 0.005, 0.004,
          0.003, 0.002, 0.001, 0.001, 0.0001)
        critval <- (12:46)/10
        p <- approx(critval, pval, x, rule=2)$y
        p <- 1 - (1-p)^k
      } else {
        p <- ifelse(x<0.1, 1,
        {
          summand <- function(a,b)
          {
            exp(-2*(a^2)*(b^2))*(-1)^a
          }
          p <- 0
          for(i in 1:100) p <- p + summand(i,x)
          1-(1+2*p)^k
        })
      }
    },

    "range" = {
      p <- ifelse(x<0.4,1,
      {
        p <- 0
        for(i in 1:10) p <- p + (4*i^2*x^2 - 1) * exp(-2*i^2*x^2)
        1-(1-2*p)^k
      })
    },

    "maxL2" = {
      if(k > 25) {
        warning("number of regressors > 25, critical values for 25 regressors used")
        k <- 25
      }
      critval <- get("sc.maxL2")[as.character(k), ]
      p <- approx(c(0, critval), c(1, 1-as.numeric(names(critval))), x, rule=2)$y
    },

    "meanL2" = {
      if(k > 25) {
        warning("number of regressors > 25, critical values for 25 regressors used")
        k <- 25
      }
      critval <- get("sc.meanL2")[as.character(k), ]
      p <- approx(c(0, critval), c(1, 1-as.numeric(names(critval))), x, rule=2)$y
    })
  },

  "Brownian motion increments" = {
    if(alt.boundary) stop("alternative boundaries for Brownian motion increments not available")
    if(functional != "max") stop("only max functional for Brownian motion increments available")

    crit.table <- cbind(c(3.2165, 2.9795, 2.8289, 2.7099, 2.6061, 2.5111, 2.4283, 2.3464, 2.2686, 2.2255),
                        c(3.3185, 3.0894, 2.9479, 2.8303, 2.7325, 2.6418, 2.5609, 2.4840, 2.4083, 2.3668),
                        c(3.4554, 3.2368, 3.1028, 2.9874, 2.8985, 2.8134, 2.7327, 2.6605, 2.5899, 2.5505),
                        c(3.6622, 3.4681, 3.3382, 3.2351, 3.1531, 3.0728, 3.0043, 2.9333, 2.8743, 2.8334),
                        c(3.8632, 3.6707, 3.5598, 3.4604, 3.3845, 3.3102, 3.2461, 3.1823, 3.1229, 3.0737),
                        c(4.1009, 3.9397, 3.8143, 3.7337, 3.6626, 3.5907, 3.5333, 3.4895, 3.4123, 3.3912))
    tablen <- dim(crit.table)[2]
    tableh <- (1:10)*0.05
    tablep <- c(0.2, 0.15, 0.1, 0.05, 0.025, 0.01)

    ## crit.table gives Table 1 of Chu, Hornik, Kuan (1995)
    ## but the corresponding test statistic is scaled differently
    ## by a factor of sqrt(h).
    crit.table <- crit.table * sqrt(tableh)

    tableipl <- numeric(tablen)
    for(i in (1:tablen)) tableipl[i] <- approx(tableh, crit.table[,i], h, rule = 2)$y
    p <- approx(c(0,tableipl), c(1,tablep), x, rule = 2)$y
    p <- 1 - (1-p)^k
  },

  "Brownian bridge increments" = {
    if(k > 6) k <- 6
    switch(functional,
    "max" = {
      crit.table <- get("sc.me")[(((k-1)*10+1):(k*10)), ]
      tablen <- dim(crit.table)[2]
      tableh <- (1:10)*0.05
      tablep <- c(0.1, 0.05, 0.025, 0.01)
      tableipl <- numeric(tablen)
      for(i in (1:tablen)) tableipl[i] <- approx(tableh, crit.table[,i], h, rule = 2)$y
      p <- approx(c(0,tableipl), c(1,tablep), x, rule = 2)$y
    },
    "range" = {
      p <- ifelse(x<0.53,1,
      {
        p <- 0
        for(i in 1:10) p <- p + (-1)^(i-1) * i * pnorm(-i*x)
        1-(1-8*p)^k
      })
    })
  })

  return(p)
}


print.efp <- function(x, ...)
{
    cat("\nEmpirical Fluctuation Process:", x$type.name, "\n\n")
    cat("Call: ")
    print(x$call)
    cat("\n")
}

sctest <- function(x, ...)
{
    UseMethod("sctest")
}

sctest.formula <- function(formula, type = c("Rec-CUSUM", "OLS-CUSUM",
  "Rec-MOSUM", "OLS-MOSUM", "RE", "ME", "fluctuation", "Score-CUSUM", "Nyblom-Hansen",
  "Chow", "supF", "aveF", "expF"), h = 0.15, alt.boundary = FALSE, functional = c("max", "range",
  "maxL2", "meanL2"), from = 0.15, to = NULL, point = 0.5, asymptotic = FALSE, data = list(), ...)
{
  type <- match.arg(type)
  if(type == "fluctuation") type <- "RE"
  functional <- match.arg(functional)
  dname <- paste(deparse(substitute(formula)))

  if(type == "Nyblom-Hansen") {
    type <- "Score-CUSUM"
    functional <- "meanL2"
  }

  if(type == "Chow")
  {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    modelterms <- terms(formula, data = data)
    X <- model.matrix(modelterms, data = data)
    METHOD <- "Chow test"
    k <- ncol(X)
    n <- length(y)

    ytsp <- NULL
    if(is.ts(data)){
        if(NROW(data) == n) {
          ytime <- time(data)
          ytsp <- tsp(data)
	}
    } else {
        env <- environment(formula)
        if(missing(data)) data <- env
        orig.y <- eval(attr(terms(formula), "variables")[[2]], data, env)
        if(is.ts(orig.y) & NROW(orig.y) == n){
            ytime <- time(orig.y)
            ytsp <- tsp(orig.y)
        }
    }

    ts.eps <- getOption("ts.eps")

    if(length(point) > 1) {
      if(!is.null(ytsp) && point[2] <= ytsp[3])
        point <- which(abs(ytime-(point[1]+(point[2]-1)/ytsp[3])) < ts.eps)
      else
        stop(paste(sQuote("point"), "does not specify a valid potential change point"))
    }
    else {
      if(point < 1)
      {
        point <- floor(point*n)
      }
    }

    if(!((point>k) & (point<(n-k)))) stop("inadmissable change point")
    e <- lm.fit(X,y)$residuals
    u <- c(lm.fit(as.matrix(X[(1:point),]),y[1:point])$residuals,
      lm.fit(as.matrix(X[((point+1):n),]),y[((point+1):n)])$residuals)
    STATISTIC <- ((sum(e^2)-sum(u^2))/k)/((sum(u^2))/(n-2*k))
    names(STATISTIC) <- "F"
    if(asymptotic)
    {
      STATISTIC <- STATISTIC * k
      PVAL <- 1 - pchisq(STATISTIC, k)
    }
    else
      PVAL <- 1 - pf(STATISTIC, k, (n-2*k))
    RVAL <- list(statistic = STATISTIC, p.value = PVAL, method = METHOD, data.name = "formula")
    class(RVAL) <- "htest"
  }
  else if(type %in% c("Rec-CUSUM", "OLS-CUSUM", "Rec-MOSUM", "OLS-MOSUM", "RE", "ME", "Score-CUSUM"))
  {
    process <- efp(formula, type = type, h = h, data = data, ...)
    RVAL <- sctest(process, alt.boundary = alt.boundary, functional = functional)
  }
  else if(type %in% c("supF", "aveF", "expF"))
  {
    fs <- Fstats(formula, from = from, to = to, data=data, ...)
    RVAL <- sctest(fs, type = type)
  }

  RVAL$data.name <- dname
  return(RVAL)
}

sctest.efp <- function(x, alt.boundary = FALSE, functional = c("max", "range", "maxL2", "meanL2"), ...)
{
    h <- x$par
    type <- x$type
    lim.process <- x$lim.process
    functional <- match.arg(functional)
    dname <- paste(deparse(substitute(x)))
    METHOD <- x$type.name
    x <- as.matrix(x$process)
    if(lim.process %in% c("Brownian motion", "Brownian bridge"))
      x <- x[-1, , drop = FALSE]
    n <- nrow(x)
    k <- ncol(x)

    if(alt.boundary) {
      if(functional == "max" & lim.process %in% c("Brownian motion", "Brownian bridge")) {
        trim <- max(round(n * 0.001, digits = 0), 1)
        j <- 1:n/n
        if(lim.process == "Brownian bridge") {
          x <- x * 1/sqrt(j * (1-j))
	  x <- x[(trim+1):(n-trim), , drop = FALSE]
        } else {
          x <- x * 1/sqrt(j)
          x <- x[trim:n, , drop = FALSE]
        }
        METHOD <- paste(METHOD, "with alternative boundaries")
      } else {
        if(lim.process == "Brownian motion") {
          j <- 1:n/n
          x <- x * 1/(1 + 2*j)
        }
        warning(paste("no alternative boundaries for", lim.process, "with", functional, "functional"))
      }
    } else {
      if(lim.process == "Brownian motion") {
        j <- 1:n/n
        x <- x * 1/(1 + 2*j)
      }
    }


    switch(functional,
      "max" = { STAT <- max(abs(x)) },
      "range" = {
        myrange <- function(y) diff(range(y))
        STAT <- max(apply(x, 2, myrange))
        METHOD <- paste(METHOD, "with range norm")
      },
      "maxL2" = {
        STAT <- max(rowSums(x^2))
        METHOD <- paste(METHOD, "with max L2 norm")
      },
      "meanL2" = {
        STAT <- mean(rowSums(x^2))
        METHOD <- paste(METHOD, "with mean L2 norm")
      })

    switch(type,
      "Rec-CUSUM" = { names(STAT) <- "S" },
      "OLS-CUSUM" = { names(STAT) <- "S0" },
      "Rec-MOSUM" = { names(STAT) <- "M" },
      "OLS-MOSUM" = { names(STAT) <- "M0" },
      "RE" = { names(STAT) <- "RE" },
      "ME" = { names(STAT) <- "ME" },
      names(STAT) <- "f(efp)")

    PVAL <- pvalue.efp(STAT, lim.process, alt.boundary, functional = functional, h = h, k = k)
    RVAL <- list(statistic = STAT, p.value = PVAL, method = METHOD, data.name = dname)
    class(RVAL) <- "htest"
    return(RVAL)
}

boundary <- function(x, ...)
{
    UseMethod("boundary")
}

boundary.efp <- function(x, alpha = 0.05, alt.boundary = FALSE, functional = "max", ...)
{
    h <- x$par
    type <- x$type
    lim.process <- x$lim.process
    functional <- match.arg(functional, c("max", "range", "maxL2", "meanL2"))
    proc <- as.matrix(x$process)
    n <- nrow(proc)
    k <- ncol(proc)

    bound <- uniroot(function(y) {pvalue.efp(y, lim.process, alt.boundary = alt.boundary,
      functional = functional, h = h, k = k) - alpha}, c(0,20))$root
    switch(lim.process,
           "Brownian motion" = {
	     if(functional == "max") {
               if(alt.boundary)
                   bound <- sqrt((0:(n-1))/(n-1))*bound
               else
                   bound <- bound + (2*bound*(0:(n-1))/(n-1))
             } else {
	       stop("only boundaries for Brownian motion with max functional available")
	     }
           },
           "Brownian bridge" = {
	     if(!(functional == "range")) {
               if(alt.boundary & functional == "max")
               {
                   j <- (0:(n-1))/(n-1)
                   bound <- sqrt(j*(1-j))*bound
               }
               else
                   bound <- rep(bound,n)
             } else {
	       stop("no boundaries Brownian bridge with range functional available")
	     }
           },
           "Brownian motion increments" = {
	     if(functional == "max") {
	       bound <- rep(bound, n)
	     } else {
               stop("only boundaries for Brownian motion increments with max functional available")
	     }
	     if(alt.boundary) warning("no alternative boundaries available for Brownian motion increments")
	   },
           "Brownian bridge increments" = {
	     if(functional == "max") {
	       bound <- rep(bound, n)
	     } else {
               stop("only boundaries for Brownian bridge increments with max functional available")
	     }
	     if(alt.boundary) warning("no alternative boundaries available for Brownian bridge increments")
	   }
    )

    bound <- ts(bound, end = end(x$process), frequency = frequency(x$process))
    return(bound)
}

lines.efp <- function(x, functional = "max", ...)
{
    if(is.null(functional)) stop("lines() cannot be added for `functional = NULL'")
    else functional <- match.arg(functional, c("max", "range", "maxL2", "meanL2"))

    if(x$lim.process %in% c("Brownian bridge", "Brownian bridge increments")) {
      z <- as.matrix(x$process)
      k <- ncol(z)

      switch(functional,
        "max" = {
          if(k > 1) {
            z <- apply(abs(z), 1, max)
            pos <- TRUE
          }
        },
        "range" = { stop("no plot available for range functional") },
        "maxL2" = {
	  if(x$lim.process == "Brownian bridge") {
            z <- rowSums(z^2)
	    pos <- TRUE
	  } else {
	    stop("no test/plot available for max L2 functional")
	  }
        },
        "meanL2" = {
	  if(x$lim.process == "Brownian bridge") {
            z <- rowSums(z^2)
	    ave <- TRUE
	    pos <- TRUE
	  } else {
	    stop("no test/plot available for mean L2 functional")
	  }
        })
      z <- ts(as.vector(z), start = start(x$process), frequency = frequency(x$process))
    } else {
      z <- x$process
    }

    lines(z, ...)
}

