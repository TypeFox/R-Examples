nlogLik <- function(object, ...) -logLik(object, ...)

gbreakpoints <- function(formula, data, fit = lm, objfun = nlogLik,
                         h = 0.15, breaks = NULL, order.by = NULL,
			 ic = c("LWZ", "BIC"), hpc = c("none", "foreach"), ...)
{
  n <- NROW(data)
  k <- attr(objfun(fit(formula, data = data, ...)), "df")    
  if(h < 1) h <- floor(n*h)
  h <- round(h, digits = 0)
  if(h > floor(n/2))
    stop("minimum segment size must be smaller than half the number of observations")
  if(is.null(breaks)) breaks <- ceiling(n/h) - 2

  ## use order.by
  if(is.null(order.by)) {
    if(is.ts(data)) order.by <- time(data)
      else if(is.zoo(data)) order.by <- index(data)
      else order.by <- 1:n  
  }
  data <- data[ORDER(order.by), , drop = FALSE]
  order.by <- order.by[ORDER(order.by)]

  hpc <- match.arg(hpc)
  if(hpc == "foreach" && !require("foreach")) {
    warning("High perfomance computing (hpc) support with 'foreach' package is not available, foreach is not installed.")
    hpc <- "none"
  }

  ## compute objfun() for all valid segmentes i, ..., j.
  
  #Z# for the default fit & objfun, this can be speeded up
  #Z# by re-using breakpoints(), otherwise this has to be done
  #Z# in a double for-loop and becomes *really* slow
  if(missing(fit) & missing(objfun)) {    
    #Z# obtain RSS-based breakpoints
    lm_bp <- breakpoints(formula, h = h, data = as.data.frame(data), hpc = hpc)
    
    RSS2obj.triang <- function(x) {
      rval <- x$RSS.triang
      n <- x$nobs
      df <- x$nreg + 1
      h <- x$nobs - length(rval)

      rss2nloglik <- function(i) {
        rss <- rval[[i]]
        ni <- seq(along = rss)
        loglik <- 0.5 * ni * (log(rss) + 1 - log(ni) + log(2 * pi))
        loglik
      }
      rval <- lapply(seq(along = rval), rss2nloglik)
      rval <- lapply(rval, function(x) x[-(1:h)])

      return(rval)
    }
    
    #Z# convert RSS to negative log-likelihood
    obj.triang <- RSS2obj.triang(lm_bp)

  } else {
    #Z# this becomes *really* slow...
    #Z# (but can be alleviated to some degree by hpc option)
    obj.triang <- list()
    if(hpc == "none") {
      for(i in 1:(n-h+1)) {
        obj_j <- rep(0, (n-i-h+2))
        for(j in (i+h-1):n) obj_j[(j-i-h+2)] <- as.numeric(objfun(fit(formula, data = data[(i:j),,drop = FALSE], ...)))
        obj.triang[[i]] <- obj_j
      }
    } else {
      obj.triang <- foreach(i = 1:(n-h+1)) %dopar% {
        obj_j <- rep(0, (n-i-h+2))
        for(j in (i+h-1):n) obj_j[(j-i-h+2)] <- as.numeric(objfun(fit(formula, data = data[(i:j),,drop = FALSE], ...)))
        obj_j
      }
    }
  }  
  #Z# For other special models (combinations of fit and objfun)
  #Z# it would be great to have external methods that generate
  #Z# something that can be easily transformed in an obj.triang
  #Z# object.

  ## function to extract the objective(i,j) from obj.triang
  objective <- function(i,j) {
    if((j < (i + h - 1)) || i < 1 || j > n) NA
      else obj.triang[[i]][j-i-h+2]
  }

  ## compute optimal previous partner if observation i is the mth break
  ## store results together with objective in obj.table

  ## breaks = 1
  index <- h:(n-h)
  obj.break <- sapply(index, function(i) objective(1,i))
  obj.table <- cbind(index, obj.break)
  rownames(obj.table) <- as.character(index)

  ## breaks >= 2

  extend.obj.table <- function(obj.table, breaks)
  {
    if((breaks*2) > ncol(obj.table)) {
      for(m in (ncol(obj.table)/2 + 1):breaks)
      {
        my.index <- (m*h):(n-h)
        my.obj.table <- obj.table[,c((m-1)*2 - 1, (m-1)*2)]
        my.obj.table <- cbind(my.obj.table, NA, NA)
        for(i in my.index)
        {
          pot.index <- ((m-1)*h):(i - h)
          obj.break <- sapply(pot.index, function(j) my.obj.table[as.character(j), 2] + objective(j+1,i))
          opt <- which.min(obj.break)
          my.obj.table[as.character(i), 3:4] <- c(pot.index[opt], obj.break[opt])
        }
        obj.table <- cbind(obj.table, my.obj.table[,3:4])
      }
      colnames(obj.table) <- as.vector(rbind(paste("break", 1:breaks, sep = ""),
                                             paste("objective", 1:breaks, sep = "")))
    }
    return(obj.table)
  }

  obj.table <- extend.obj.table(obj.table, breaks)

  ## extract optimal breaks

  extract.breaks <- function(obj.table, breaks)
  {
    if((breaks*2) > ncol(obj.table)) stop("compute obj.table with enough breaks before")
    index <- obj.table[, 1, drop = TRUE]
    obj.break <- sapply(index, function(i) obj.table[as.character(i),breaks*2] + objective(i + 1, n))
    opt <- index[which.min(obj.break)]
    if(breaks > 1) {
      for(i in ((breaks:2)*2 - 1))
        opt <- c(obj.table[as.character(opt[1]),i], opt)
    }
    names(opt) <- NULL
    return(opt)
  }

  opt <- extract.breaks(obj.table, breaks)

  RVAL <- list(breakpoints = opt,
               obj.table = obj.table,
	       obj.triang = obj.triang,
	       objective = objective,
	       extract.breaks = extract.breaks,
	       extend.obj.table = extend.obj.table,
	       nobs = n,
	       npar = k,
	       fit = fit,
	       fitname = deparse(substitute(fit)),
	       objfun = objfun,
	       objname = deparse(substitute(objfun)),
	       icname = match.arg(ic),
	       call = match.call(),
	       index = order.by)
  class(RVAL) <- c("gbreakpointsfull", "gbreakpoints", "breakpointsfull", "breakpoints")
  RVAL$breakpoints <- breakpoints(RVAL)$breakpoints
  return(RVAL)
}

breakpoints.gbreakpointsfull <- function(obj, breaks = NULL, ...)
{
  if(is.null(breaks)) breaks <- obj$icname
  if(is.character(breaks))
  {
    breaks <- match.arg(breaks, c("BIC", "LWZ"))
    IC <- summary(obj)$objective[breaks,]
    breaks <- if(any(is.na(IC))) length(IC)-1 else which.min(IC)-1
  }
  if(breaks < 1)
  {
    breakpoints <- NA
    objective <- obj$objective(1, obj$nobs)
  } else {
    obj.tab <- obj$extend.obj.table(obj$obj.table, breaks)
    breakpoints <- obj$extract.breaks(obj.tab, breaks)
    bp <- c(0, breakpoints, obj$nobs)
    objective <- sum(apply(cbind(bp[-length(bp)]+1,bp[-1]), 1,
                     function(x) obj$objective(x[1], x[2])))
  }
  RVAL <- list(breakpoints = breakpoints,
               objective = objective,
               nobs = obj$nobs,
	       npar = obj$npar,
	       call = match.call(),
	       fitname = obj$fitname,
	       objname = obj$objname,
               index = obj$index)
  class(RVAL) <- c("gbreakpoints", "breakpoints")
  return(RVAL)
}


print.gbreakpoints <- function(x, format.times = NULL, ...)
{
  if(any(is.na(x$breakpoints))) lbp <- 0
    else lbp <- length(x$breakpoints)
  cat(paste("\n\t Optimal ", lbp + 1, "-segment partition for `", x$fitname,"' fit: \n\n", sep = ""))
  cat("Call:\n")
  print(x$call)
  cat("\nBreakpoints at observation number:\n")
  cat(x$breakpoints,"\n")
  cat("\nCorresponding to breakdates:\n")
  cat(as.character(breakdates(x, format.times = format.times)),"\n")
}

breakdates.gbreakpoints <- function(obj, format.times = NULL, breaks = NULL, ...)
{
  if(is.null(format.times)) format.times <- !identical(obj$index, 1:obj$nobs)
  if("gbreakpointsfull" %in% class(obj)) obj <- breakpoints(obj, breaks = breaks)

  if(is.na(obj$breakpoints)[1])
    NA
  else {
    if(format.times) obj$index[obj$breakpoints] else obj$breakpoints/obj$nobs
  }
}

summary.gbreakpoints <- function(object, ...)
{
  print(object)
  cat(paste("\nObjective function `", object$objname, "': ", format(object$objective),"\n", sep = ""))
}

summary.gbreakpointsfull <- function(object, breaks = NULL,
  sort = TRUE, format.times = NULL, ...)
{
  if(is.null(format.times)) format.times <- !identical(object$index, 1:object$nobs)
  if(is.null(breaks)) breaks <- ncol(object$obj.table)/2
  n <- object$nobs
  objective <- c(object$objective(1, n), rep(NA, breaks))
  BIC <- c(AIC(breakpoints(object, breaks = 0), k = log(n)), rep(NA, breaks))
  LWZ <- c(AIC(breakpoints(object, breaks = 0), k = 0.299 * log(n)^2.1), rep(NA, breaks))
  names(objective) <- as.character(0:breaks)
  bp <- breakpoints(object, breaks = breaks)
  bd <- as.character(breakdates(bp, format.times = format.times))
  objective[breaks + 1] <- bp$objective
  BIC[breaks + 1] <- AIC(bp, k = log(n))
  LWZ[breaks + 1] <- AIC(bp, k = 0.299 * log(n)^2.1)
  bp <- bp$breakpoints
  if(breaks > 1) {
  for(m in (breaks-1):1)
  {
    bp <- rbind(NA, bp)
    bd <- rbind(NA, bd)
    bpm <- breakpoints(object, breaks = m)
    if(sort) {
      pos <- apply(outer(bpm$breakpoints, bp[nrow(bp),],
                   FUN = function(x,y) abs(x - y)), 1, which.min)
      if(length(pos) > unique(length(pos))) {
        warning("sorting not possible", call. = FALSE)
	sort <- FALSE
      }
    }
    if(!sort) pos <- 1:m
    bp[1,pos] <- bpm$breakpoints
    bd[1,pos] <- as.character(breakdates(bpm, format.times = format.times))
    objective[m+1] <- bpm$objective
    BIC[m+1] <- AIC(bpm, k = log(n))
    LWZ[m+1] <- AIC(bpm, k = 0.299 * log(n)^2.1)
  }} else {
    bp <- as.matrix(bp)
    bd <- as.matrix(bd)
  }
  rownames(bp) <- as.character(1:breaks)
  colnames(bp) <- rep("", breaks)
  rownames(bd) <- as.character(1:breaks)
  colnames(bd) <- rep("", breaks)
  objective <- rbind(objective, BIC, LWZ)
  rownames(objective) <- c(object$objname, "BIC", "LWZ")
  RVAL <- list(breakpoints = bp,
               breakdates = bd,
	       objective = objective,
	       objname = object$objname,
	       fitname = object$fitname,
	       icname = object$icname,
	       call = object$call)
  class(RVAL) <- "summary.gbreakpointsfull"
  return(RVAL)
}

print.summary.gbreakpointsfull <- function(x, ...)
{
  bp <- x$breakpoints
  breaks <- ncol(bp)
  bd <- x$breakdates
  objective <- x$objective
  bp[is.na(bp)] <- ""
  bd[is.na(bd)] <- ""
  rownames(bp) <- paste("m = ", rownames(bp), "  ", sep = "")
  rownames(bd) <- paste("m = ", rownames(bd), "  ", sep = "")
  objective <- rbind(0:(ncol(objective) - 1), format(objective))
  rownames(objective) <- c("m", x$objname, "BIC", "LWZ")
  colnames(objective) <- rep("", breaks + 1)

  cat("\n\t Optimal (m+1)-segment partition for `", x$fitname,"' : \n\n", sep = "")
  cat("Call:\n")
  print(x$call)
  cat("\nBreakpoints at observation number:\n")
  print(bp, quote = FALSE)
  cat("\nCorresponding to breakdates:\n")
  print(bd, quote = FALSE)
  cat("\nFit:\n")
  print(objective, quote = FALSE)
  invisible(x)
}

plot.gbreakpointsfull <- function(x, breaks = NULL, ...)
{
  plot(summary(x, breaks = breaks), ...)
}

plot.summary.gbreakpointsfull <- function(x, type = "b", col = c(1, 4), lty = c(1, 1),
  legend = TRUE, main = NULL, xlab = "Number of breakpoints", ylab = NULL, ...)
{
  breaks <- as.numeric(colnames(x$objective))
  obj <- x$objective[x$objname,]
  icname <- x$icname
  IC <- x$objective[icname,]
  
  col <- rep(col, length.out = 2)
  lty = rep(lty, length.out = 2)
  
  if(any(is.na(IC))) {
    if(is.null(main)) main <- ""
    if(is.null(ylab)) ylab <- x$objname
    plot(breaks, obj, ylab = ylab, xlab = xlab, main = main,
         type = type, col = col[1], lty = lty[1], ...)
  } else {
    if(is.null(main)) main <- paste(icname, "and Negative Log-Likelihood")
    if(is.null(ylab)) ylab <- ""
    plot(breaks, IC, ylab = ylab, xlab = xlab, main = main,
	 type = type, col = col[1], lty = lty[1], ...)
    par(new = TRUE)
    plot(breaks, obj, type = type, axes = FALSE, col = col[2], lty = lty[2],
         xlab = "", ylab = "")
    if(legend) legend("topright", legend = c(icname, "neg. Log-Lik."),
                      lty = lty, col = col, bty = "n")
    axis(4)
    par(new = FALSE)
  }
}


logLik.gbreakpoints <- function(object, ...)
{
  if(is.null(object$npar)) return(NA)
  n <- object$nobs
  nb <- length(object$breakpoints[!is.na(object$breakpoints)])
  df <- object$npar * (nb + 1) + nb
  logL <- -object$objective
  attr(logL, "df") <- df
  class(logL) <- "logLik"
  return(logL)
}

logLik.gbreakpointsfull <- function(object, breaks = NULL, ...)
{
  bp <- breakpoints(object, breaks = length(object$breakpoints))
  logL <- logLik(bp)
  return(logL)
}

AIC.gbreakpoints <- function(object, ..., k = 2)
{
  if(is.null(object$npar)) return(NA)
  NextMethod()
}

AIC.gbreakpointsfull <- function(object, breaks = NULL, ..., k = 2)
{
  if(is.null(object$npar)) return(NA)
  if(is.null(breaks)) breaks <- 0:(ncol(object$obj.table)/2)
  RVAL <- NULL
  for(m in breaks)
    RVAL <- c(RVAL, AIC(breakpoints(object, breaks = m), k = k))
  names(RVAL) <- breaks
  return(RVAL)
}

confint.gbreakpoints <- function(object, parm, level = 0.95, ...)
{
  stop(paste("no confidence intervals available for", dQuote("gbreakpoints"), "objects"))
}
