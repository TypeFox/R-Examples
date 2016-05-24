densityMclust <- function(data, ...) 
{
  mc <- match.call()
  obj <- Mclust(data, ...)
  obj$call <- mc
  d <- dens(modelName = obj$modelName, data = data, 
            parameters = obj$parameters, logarithm = FALSE)
  obj$density <- d
  class(obj) <- c("densityMclust", "Mclust")
  return(obj)
}

predict.densityMclust <- function(object, newdata, what = c("dens", "cdens"), ...)
{
  if(!inherits(object, "densityMclust")) 
    stop("object not of class \"densityMclust\"")
  if(missing(newdata))
    { newdata <- object$data }
  what <- match.arg(what)
  
  if(what == "dens")
    { d <- dens(modelName = object$modelName, 
                data = newdata, 
                parameters = object$parameters) }
  else
    { d <- cdens(modelName = object$modelName, 
                 data = newdata, 
                 parameters = object$parameters)
      dim <- dim(d)
      attributes(d) <- NULL 
      d <- array(d, dim)
    }
  return(d)
}

plot.densityMclust <- function(x, data = NULL, what = c("BIC", "density", "diagnostic"), ...) 
{
  object <- x # Argh.  Really want to use object anyway

  what <- match.arg(what, several.ok = TRUE)
  if(object$d > 1) 
    what <- setdiff(what, "diagnostic")
  oldpar <- par(no.readonly = TRUE)
  # on.exit(par(oldpar))
  
  plot.densityMclust.density <- function(...)
  { 
    if(object$d == 1)      plotDensityMclust1(object, data = data, ...)
    else if(object$d == 2) plotDensityMclust2(object, data = data, ...)
    else                   plotDensityMclustd(object, data = data, ...)
  }
  
  plot.densityMclust.bic <- function(...)
  { 
    # this add right axis for bic diff
    # oldpar <- par(no.readonly = TRUE)
    # on.exit(par(oldpar))
    # mar <- oldpar$mar
    # mar[4] <- max(mar[4],3)
    # par(mar = mar)
    # plot.mclustBIC(object$BIC, ...)
    # yaxp <- par("yaxp")
    # bicdiff <- seq(0, yaxp[1] - object$bic, length = 100)
    # bicdiff <- pretty(bicdiff, yaxp[3]+1)
    # axis(4, at = object$bic+bicdiff, labels = signif(bicdiff,2))
    plot.mclustBIC(object$BIC, ...)
  }
  
  plot.densityMclust.diagnostic <- function(...)
  { densityMclust.diagnostic(object, ...) }
  
  if(interactive() & length(what) > 1)
    { title <- "Model-based density estimation plots:"
      # present menu waiting user choice
      choice <- menu(what, graphics = FALSE, title = title)
      while(choice != 0)
           { if(what[choice] == "BIC")        plot.densityMclust.bic(...)
             if(what[choice] == "density")    plot.densityMclust.density(...)
             if(what[choice] == "diagnostic") plot.densityMclust.diagnostic(...)
             # re-present menu waiting user choice
             choice <- menu(what, graphics = FALSE, title = title)
           }
  } 
  else 
    { if(any(what == "BIC"))        plot.densityMclust.bic(...)
      if(any(what == "density"))    plot.densityMclust.density(...)
      if(any(what == "diagnostic")) plot.densityMclust.diagnostic(...)
  }
 
  invisible()
}


plotDensityMclust1 <- function(x, data = NULL, hist.col = "lightgrey", hist.border = "white", breaks = "Sturges", ...) 
{
  object <- x # Argh.  Really want to use object anyway
  mc <- match.call(expand.dots = TRUE)
  mc$x <- mc$data <- mc$hist.col <- mc$hist.border <- mc$breaks <- NULL
  xlab <- mc$xlab
  if(is.null(xlab)) 
    xlab <- deparse(object$call$data)
  ylab <- mc$ylab
  if(is.null(ylab)) 
    ylab <- "Density"
  #
  xrange <- extendrange(object$data, f = 0.1)
  xlim <- eval(mc$xlim, parent.frame())
  if(!is.null(xlim)) 
    xrange <- range(xlim)
  ylim <- eval(mc$ylim, parent.frame())
  #
  eval.points <- seq(from = xrange[1], to = xrange[2], length = 1000)
  d <- predict.densityMclust(object, eval.points)
  #
  if(!is.null(data)) 
  { h <- hist(data, breaks = breaks, plot = FALSE)
    plot(h, freq = FALSE, col = hist.col, border = hist.border, main = "",
         xlim = range(h$breaks, xrange), 
         # ylim = range(0, ylim, h$density, max(d)+diff(range(d))*0.1),
         ylim =  if(!is.null(ylim)) range(ylim) else range(0, h$density, d),
         xlab = xlab, ylab = ylab)
    box()
    mc[[1]] <- as.name("lines")
    mc$x <- eval.points
    mc$y <- d
    mc$type <- "l"
    eval(mc, parent.frame())
  }
  else
  { mc[[1]] <- as.name("plot")
    mc$x <- eval.points
    mc$y <- d
    mc$type <- "l"
    mc$xlim <- xlim
    mc$ylim <- if(!is.null(ylim)) range(ylim) else range(0, d)
    mc$ylab <- ylab
    mc$xlab <- xlab
    eval(mc, parent.frame())
  }
  invisible()
}

plotDensityMclust2 <- function(x, data = NULL, nlevels = 11, levels = NULL, col = grey(0.6), 
                               points.pch = 1, points.col = 1, points.cex = 0.8, ...) 
{
  # This function call surfacePlot() with a suitable modification of arguments
  object <- x # Argh.  Really want to use object anyway
  mc <- match.call(expand.dots = TRUE)
  mc$x <- mc$points.pch <- mc$points.col <- mc$points.cex <- NULL
  mc$nlevels <- nlevels; mc$levels <- levels
  mc$col <- col
  
  if(is.null(data)) 
    { addPoints <- FALSE
      mc$data <- object$data } 
  else
    { addPoints <- TRUE }
  
  # set mixture parameters
  par <- object$parameters
  # these parameters should be missing 
  par$variance$cholSigma <- par$Sigma <- par$Vinv <- NULL
  if(is.null(par$pro)) par$pro <- 1  # LS: bug?
  #par$variance$d <- 2 # LS: bug?
  par$variance$cholsigma <- par$variance$sigma
  for(k in seq(par$variance$G))
     { par$variance$cholsigma[,,k] <- chol(par$variance$sigma[,,k]) }
  mc$parameters <- par
  # now surfacePlot() is called
  mc[[1]] <- as.name("surfacePlot")
  out <- eval(mc, parent.frame())
  if(addPoints)
    points(data, pch = points.pch, col = points.col, cex = points.cex)
  #
  invisible(out)
}

plotDensityMclustd <- function(x, data = NULL, nlevels = 11, levels = NULL, col = grey(0.6), 
                               points.pch = 1, points.col = 1, points.cex = 0.8, gap = 0.2, ...) 
{
  # This function call surfacePlot() with a suitable modification of arguments
  
  object <- x # Argh.  Really want to use object anyway
  mc <- match.call(expand.dots = TRUE)
  mc$x <- mc$points.pch <- mc$points.col <- mc$points.cex <- mc$gap <- NULL
  mc$nlevels <- nlevels; mc$levels <- levels
  mc$col <- col
  
  if(is.null(data)) 
    { data <- mc$data <- object$data
      addPoints <- FALSE }
  else
    { data <- as.matrix(data)
      addPoints <- TRUE  }
  
  nc <- object$d
  oldpar <- par(mfrow = c(nc, nc), 
                mar = rep(c(gap,gap/2),each=2), 
                oma = c(4, 4, 4, 4),
                no.readonly = TRUE)
  on.exit(par(oldpar))

  for(i in seq(nc))
     { for(j in seq(nc)) 
          { if(i == j) 
              { plot(0,0,type="n",xlab="",ylab="",axes=FALSE)
                text(0,0, colnames(data)[i], cex=1.5, adj=0.5)
                box()
            } 
            else 
              { # set mixture parameters
                par <- object$parameters
                if(is.null(par$pro)) par$pro <- 1
                par$mean <- par$mean[c(j,i),,drop=FALSE]
                par$Vinv <- NULL
                par$variance$d <- 2
                sigma <- array(dim = c(2, 2, par$variance$G))
                for(g in seq(par$variance$G))
                  sigma[,,g] <- par$variance$sigma[c(j,i),c(j,i),g]
                par$variance$sigma <- sigma
                par$variance$Sigma <- NULL
                par$variance$cholSigma <- NULL
                par$variance$cholsigma <- NULL
                mc$parameters <- par
                mc$data <- object$data[,c(j,i)]
                mc$axes <- FALSE
                mc[[1]] <- as.name("surfacePlot")
                out <- eval(mc, parent.frame())
                box()
                if(addPoints)
                  points(data[,c(j,i)], pch = points.pch, 
                         col = points.col, cex = points.cex)
              }
              if(i == 1 && (!(j%%2))) axis(3)
              if(i == nc && (j%%2))   axis(1)
              if(j == 1 && (!(i%%2))) axis(2)
              if(j == nc && (i%%2))   axis(4)
          }
  }
  #
  invisible(out)
}

dens <- function(modelName, data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  aux <- list(...)
  cden <- cdens(modelName = modelName, data = data,
                logarithm = TRUE, parameters = parameters, warn = warn)
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  G <- if(oneD) { length(parameters$mean) }
       else     { ncol(as.matrix(parameters$mean)) }
  pro <- parameters$pro
  if(is.null(pro))
    stop("mixing proportions must be supplied")
  noise <- (!is.null(parameters$Vinv))
  if(G > 1)
    { if(noise) 
        { # proN <- pro[length(pro)]
          pro <- pro[-length(pro)]
          # pro <- pro/sum(pro)
      }
      if(any(proz <- pro == 0)) 
        { pro <- pro[!proz]
          cden <- cden[, !proz, drop = FALSE]
      }
      cden <- sweep(cden, 2, FUN = "+", STATS = log(pro))
  }
  # logsumexp
  maxlog <- apply(cden, 1, max)
  cden <- sweep(cden, 1, FUN = "-", STATS = maxlog)
  den <- logb(apply(exp(cden), 1, sum)) + maxlog
  if(noise) 
    den <- den + parameters$pro[G+1]*parameters$Vinv
  if(!logarithm) den <- exp(den)
  den
}

cdens <- function(modelName, data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  modelName <- switch(EXPR = modelName,
                      X = "E",
                      XII = "EII",
                      XXI = "EEI",
                      XXX = "EEE",
                      modelName)
  checkModelName(modelName)
  funcName <- paste("cdens", modelName, sep = "")
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name(funcName)
  mc$modelName <- NULL
  eval(mc, parent.frame())
}

densityMclust.diagnostic <- function(object, type = c("cdf", "qq"), 
                                     col = c("black", "green4"), 
                                     lwd = c(2,2), lty = c(1,2),
                                     legend = TRUE, grid = TRUE, main = TRUE, ...)
{
# Diagnostic plots for density estimation 
# (only available for the one-dimensional case)
# 
# Arguments:
# object = a 'densityMclust' object
# data = the data vector
# type = type of diagnostic plot:
# cdf = the fitted distribution function vs the empirical distribution function;
# qq = the fitted distribution function evaluated over the observed points vs 
#      the quantile from a uniform distribution.
#
# Reference: 
# Loader C. (1999), Local Regression and Likelihood. New York, Springer, 
#   pp. 87-90)

  if(!any(class(object) == "densityMclust"))
    { stop("first argument must be an object of class 'densityMclust'") }
  if(object$d > 1)
    { warning("only available for one-dimensional data") 
      return() }  
  type <- match.arg(type, c("cdf", "qq"), several.ok = TRUE)
  main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)

  data <- as.numeric(object$data)
  n <- length(data)
  cdf <- cdfMclust(object, data = data, ngrid = min(n*10,1000), ...)
  
  oldpar <- par(no.readonly = TRUE)
  if(interactive() & length(type) > 1) 
    { par(ask = TRUE)
      on.exit(par(oldpar)) }
  
  if(any(type == "cdf"))
  { # Fitted CDF vs Emprical CDF    
    empcdf <- ecdf(data)
    plot(empcdf, do.points = FALSE, 
         col = col[2], lwd = lwd[2], lty = lty[2],
         xlab = deparse(object$call$data), 
         ylab = "Cumulative Distribution Function",
         panel.first = if(grid) grid(equilogs=FALSE) else NULL,
         main = NULL, ...)
    if(main) title(main = "CDF plot", cex.main = 1.1)
    lines(cdf, col = col[1], lwd = lwd[1], lty = lty[1])
    rug(data)
    if(legend)
      { legend("bottomright", legend = c("Est.CDF", "Emp.CDF"), 
               ncol = 1, inset = 0.05, cex = 0.8,
               col = col, lwd = lwd, lty = lty) }
  }
  
  if(any(type == "qq"))
  { # Q-Q plot
    q <- quantileMclust(object, p = ppoints(n))
    plot(q, sort(data),
         xlab = "Quantiles from estimated density", 
         ylab = "Sample Quantiles", 
         panel.first = if(grid) grid(equilogs=FALSE) else NULL,
         main = NULL, ...)
    if(main) title(main = "Q-Q plot", cex.main = 1.1)
    with(list(y = sort(data), x = q),
         { i <- (y > quantile(y, 0.25) & y < quantile(y, 0.75))
           abline(lm(y ~ x, subset = i), lty = 2) 
         })
# P-P plot
# cdf <- cdfMclust(object, data, ...)
# plot(seq(1,n)/(n+1), cdf$y, xlab = "Uniform quantiles", 
#    ylab = "Cumulative Distribution Function",
#      main = "Diagnostic: P-P plot")
# abline(0, 1, lty = 2)
  }

  invisible()
} 

cdfMclust <- function(object, data, ngrid = 100, ...)
{
# Cumulative Density Function
# (only available for the one-dimensional case)
#
# Returns the estimated CDF evaluated at points given by the optional
# argument data. If not provided, a regular grid of ngrid points is used. 
#
# Arguments:
# object = a 'densityMclust' object
# data = the data vector
# ngrid = the length of rectangular grid 
  
  if(!any(class(object) == "densityMclust"))
    { stop("first argument must be an object of class 'densityMclust'") }
  
  if(missing(data))
    { eval.points <- extendrange(object$data, f = 0.1)
      eval.points <- seq(eval.points[1], eval.points[2], length.out = ngrid) }
  else
    { eval.points <- sort(as.vector(data))
      ngrid <- length(eval.points) }
  
  G <- object$G
  pro <- object$parameters$pro
  mean <- object$parameters$mean
  var <- object$parameters$variance$sigmasq
  if(length(var) < G) var <- rep(var, G)
  noise <- (!is.null(object$parameters$Vinv))

  cdf <- rep(0, ngrid)
  for(k in seq(G))
     { cdf <- cdf + pro[k]*pnorm(eval.points, mean[k], sqrt(var[k])) }
  if(noise) 
    cdf <- cdf/sum(pro[seq(G)])
  
  out <- list(x = eval.points, y = cdf)    
  return(out)
}

quantileMclust <- function(object, p, ...)
{
# Calculate the quantile of a univariate mixture corresponding to cdf equal to p 
#
# Arguments:
# object = a 'densityMclust' object
# p = vector of probabilities (0 <= p <= 1)

  if(!any(class(object) == "densityMclust"))
    { stop("first argument must be an object of class 'densityMclust'") }
  
  eval.points <- extendrange(object$data, f = 1)
  eval.points <- seq(eval.points[1], eval.points[2], length.out = 10000) 
  cdf <- cdfMclust(object, data = eval.points)
  q <- spline(cdf$y, cdf$x, method = "fmm", xmin = 0, xmax = 1, xout = p)$y
  q[ p < 0 | p > 1] <- NaN
  q[ p == 0 ] <- -Inf
  q[ p == 1 ] <- Inf
  return(q)  
}