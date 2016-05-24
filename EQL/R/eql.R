.eql <- function(y, mu, phi, v, d) {
  value <- -1 / 2 * log(2 * pi * phi * v) - d / (2 * phi)
  return(value)
}
    

eql <- function(formula, param.space, family=powerVarianceFamily(),
                phi.method=c("pearson", "mean.dev"),
                include.model=TRUE, smooth.grid=10,
                do.smooth=dim(family)==1, verbose=1, ...) {
  if (length(param.space) != dim(family) || !is.list(param.space)) {
    stop("'param.space' must be a list including a vector for each ",
         "parameter of the variance family!")
  }
  if (class(family) != "varianceFamily") {
    stop("'family' must be an object of class 'varianceFamily'!")
  }
  if (do.smooth && dim(family) > 1) {
    warning("Smoothing currently available only for one dimensional ",
            "variance families!")
    do.smooth <- FALSE
  }
  if (smooth.grid <= 0 || !isInteger(smooth.grid)) {
    warning("'smooth.grid' must be an integer greater or equal than 1!")
    smooth.grid <- 10
  }
  if (!do.smooth) {
    smooth.grid <- NULL
  }
  params <- family$params
  param.space <- lapply(param.space, sort)
  paramValues <- expand.grid(param.space)
  if (is.null(names(param.space))) {
    names(paramValues) <- names(params)
  } else if (any(!(names(param.space) %in% names(params)))) {
    stop("Unknown parameters in 'param.space'!\nNeeded arguments:",
         paste("'", names(params), "'", sep="", collapse=", "), "\nPassed arguments:",
         paste("'", names(param.space), "'", sep="", collapse=", "), "!")
  }
  nSpace <- dim(paramValues)[1]
  phi.method <- match.arg(phi.method)
  verbose <- as.numeric(verbose)
  if (verbose >= 1) {
    cat("~~ Extended Quasi Likelihood ~~\n")
    cat("* Dimension of parameter space: ", nSpace, "\n")
    cat("* Starting loop:\n")
  }
  if (verbose == 1) {
    if (nSpace >= 20) {
      cat("\n0%       50%        100%\n")
      cat("+---------+---------+\n")
      cat("|")
    }
  }
  getEQL <- function(params, total) {
    params <- as.list(params)
    nr <- params$nr
    progress <- params$progress
    params$nr <- NULL
    params$progress <- NULL
    paramNames <- names(params)
    fam <- do.call(family$family, params)
    arg <- merge(list(formula=formula, family=fam, y=TRUE), list(...))
    model <- try(do.call(glm, arg), silent=TRUE)
    if (any(class(model) == "try-error")) {
      warning("A model could not be evaluated!\n  ", model, call.=FALSE)
      eql <- NA
    } else {
      mu <- model$fitted
      y <- model$y
      d <- fam$dev.resids(y, mu, model$prior.weights)
      phi <- switch(phi.method,
                    "pearson"  = summary(model)$dispersion,
                    "mean.dev" = sum(d) / model$df.residual)
      v <- fam$variance(y)
      if (any(v == 0)) {
        warning("Zero values in the variance function! EQL cannot be evaluated!",
                call.=FALSE)
      }
      eql <- sum(.eql(y, mu, phi, v, d))
    }
    if (verbose == 1 & !is.null(progress)) {
      if (progress) {
        cat("|")
      }
      if (nr == total) {
        cat("\n")
      }
    }
    if (verbose >= 2 | (verbose == 1 & is.null(progress))) {
      if (nr == 1) {
        cat("\n")
      }
      nLength <- length(paramNames)
      fmtParams <- format(unlist(params), digits=3, nsmall=3)
      firstLine <- paste("\t+ ", paramNames[1], " = ", fmtParams[1],
                         "\tEQL = ", format(eql, digits=4, nsmall=4),
                         " (", nr, "/", total , ")",
                         "\n", sep="")
      cat(firstLine)
      if (nLength >= 2) {
        rest <- paste("\t ", paste(paramNames[2:nLength], "=",
                                   fmtParams[2:nLength],
                                   collapse="\n\t"), "\n")
        cat(rest)
      }
    }
    return(eql)
  }
  paramValues$nr <- 1:nSpace
  if (nSpace >= 20) {
    progress <- rep(FALSE, nSpace)
    progress[round(0.05 * 1:20 * nSpace)] <- TRUE
    paramValues$progress <- progress
  }
  eqls <- unlist(apply(paramValues, 1, getEQL, total=nSpace))
  paramValues$nr <- NULL
  paramValues$progress <- NULL
  if (verbose >= 1) {
    cat("\n* End loop\n")
  }
  
  if (all(is.nan(eqls) | is.infinite(eqls) | is.na(eqls))) {
    stop("No valid eql values! Try a different parameter space!")
  }
  # code snippet from tweedie::tweedie.profile
  if (do.smooth) {
    if (any(is.nan(eqls)) | any(is.infinite(eqls)) | any(is.na(eqls))) {
      retain.these <- !((is.nan(eqls)) | (is.infinite(eqls)) | 
                        (is.na(eqls)))
      paramValues <- subset(paramValues, retain.these)
      eqls <- eqls[retain.these]
      if (verbose >= 1) {
        cat("! Smooth perhaps inaccurate--eql contains Inf or NA.\n")
      }
    }
    if (length(eqls) <= 4) {
      if (verbose >=1) {
        cat("! Less than 5 valid EQL-values--skip smoothing!\n")
      } else {
        warning("Less than 5 valid EQL-values--skip smoothing!\n")
      }
    } else {
      if (verbose >= 1) {
        cat("* Start smoothing: \n")
        cat("\tSmoothing on an equidistant grid between EQL values with", smooth.grid,
            "nodes\n")
      }
      nodes <- unlist(paramValues)
      ss <- splinefun(nodes, eqls)
      nodes <- cbind(nodes[-length(nodes)], nodes[-1])
      paramValuesNew <- unique(c(apply(nodes, 1, function(row)
                                       seq(row[1], row[2], length=smooth.grid))))
      origValues <- paramValuesNew %in% nodes
      eqls <- ss(paramValuesNew)
      paramValues <- as.data.frame(paramValuesNew)
      names(paramValues) <- names(params)
      if (verbose >= 1) {
        cat("* End smoothing\n")
      }
    }
  } else {
    # no smooth
    origValues <- rep(TRUE, nSpace)
  }
  # end snippet
  eqlMax <- max(eqls, na.rm=TRUE)
  paramMax <- unique(subset(paramValues, eqls == eqlMax))
  value <- list(eql=eqls, param=paramValues, eql.max=eqlMax, param.max=paramMax,
                dim=dim(family), smooth=do.smooth, is.smoothed=!origValues,
                smooth.grid=smooth.grid)
  if (include.model) {
    if (dim(paramMax)[1] == 1) {
      family <- do.call(family$family, as.list(paramMax))
      arg <- merge(list(formula=formula, family=family), list(...))
      model <- list(model=do.call(glm, arg))
    } else {
      model <- apply(paramMax, 1, function(p) {
        family <- do.call(family$family, as.list(p))
        return(glm(formula, family, ...))
      })
    }
    value <- c(value, model)
  }
  class(value) <- "eql"
  return(value)
}

dim.eql <- function(x) {
  return(x$dim)
}

length.eql <- function(x) {
  return(length(x$eql))
}

print.eql <- function(x,...) {
  cat("EQL-Max:", x$eql.max)
  fmtParams <- format(unlist(x$param.max), digits=3, nsmall=3)
  cat("\tat:", paste(names(x$param.max), "=", fmtParams, collapse="\t"))
  cat("\nEQL Summary:\n")
  print(summary(x$eql))
  invisible(x)
}

plot.eql <- function(x, do.points=(dim(x) == 1 &&
                                   sum(!x$is.smoothed) <= 20),
                     do.ci=TRUE, alpha=0.95, do.bw=TRUE,
                     show.max=TRUE, ...) {
  extra.args <- list(...)
  if (x$dim == 1) {
    if (do.bw) {
      if (!is.null(extra.args$col)) {
        warning("'do.bw' is set to 'TRUE', skipping argument 'col'!")
      }
      extra.args$col <- col <- 1
    } else {
      col <- 2
    }
    default.args <- list(main="Profile EQL-Plot", xlab=names(x$param), ylab="EQL",
                         col=col, lwd=2, type="l", pch=19, cex.sub=0.75,
                         sub=paste("n = ", length(x$eql), ifelse(x$smooth,
                           " (smoothed)", ""), sep=""))
    args <- c(list(unlist(x$param), x$eql), merge(extra.args, default.args))
    do.call("plot", args)
    if (do.points) {
      if (sum(!x$is.smoothed) > 20) {
        warning("There are more than 20 nodes! Plot may be overloaded!")
      }
      points(unlist(x$param)[!x$is.smoothed], x$eql[!x$is.smoothed],
             col=args$col, pch=args$pch)
    }
    if (show.max) {
      param.max <- unlist(x$param.max)
      eql.min <- par("usr")[3]
      segments(param.max, eql.min, param.max, x$eql.max, lty=3, col=args$col)
      if (length(x$eql) < 50) {
        warning("There are only ", length(x$eql), " EQL values! Maximum may be inaccurate!")
      }
    }
    if (do.ci) {
      civ <- x$eql.max - 1 / 2 * qchisq(alpha, 1)
      
      #code snippet from MASS:::boxcox.default
      plims <- par("usr")
      y0 <- plims[3]
      scal <- (1/10 * (plims[4] - y0))/par("pin")[2]
      scx <- (1/10 * (plims[2] - plims[1]))/par("pin")[1]
      text(min(x$param) + scx, civ + scal, paste(alpha*100, "%", sep=""))
      #end snippet
      ind <- range(which(x$eql >= civ))
      civ.x <- x$param[c(ind, ind + c(-1, 1)),]
      civ.y <- x$eql[c(ind, ind + c(-1, 1))]
      z <- civ.x[3:2] + (civ.y[3:2] - civ) * (civ.x[c(1, 4)] - civ.x[3:2]) /
        (civ.y[3:2] - civ.y[c(1, 4)])
      abline(v=z, lty=3, col=args$col)
      abline(h=civ, lty=3, col=args$col)
     }
  } else if (x$dim == 2) {
    eql.panel <- function(..., at, contour=FALSE, labels=NULL) {
      panel.levelplot(..., at=at, contour=contour, labels=labels)
      if (show.max) {
        lpoints(x$param.max, pch=args$pch, lwd=args$lwd, col=args$col)
      }
      if (do.points) {
        lpoints(x$param, pch=args$pch, col=args$col)
      }
      if (do.ci) {
        # Bug in panel.contourplot, does not draw negative contourlines
        # add NA to circumvent this problem
        panel.contourplot(..., contour=TRUE,
                          labels=list(labels=c(paste(alpha*100, "%", sep=""), NULL),
                            col=c(args$col, NA)),
                          lty=c(3, NA), at = c(x$eql.max - 1 / 2 * qchisq(alpha, 2), NA))
      }
    }
    if (do.bw) {
      if (!is.null(extra.args$col)) {
        warning("'do.bw' is set to 'TRUE', skipping argument 'col'!")
      }
      if (!is.null(extra.args$col.regions)) {
        warning("'do.bw' is set to 'TRUE', skipping argument 'col.regions'!")
      }
      extra.args$col <- col <- 1
      col.regions <- rev(colorRampPalette(c("white","black"))(50))
      extra.args$col.regions <- col.regions
    } else {
      col.regions <- heat.colors(50)
      col <- col.regions[1]
    }
    default.args <- list(main="Contour EQL-Plot", xlab=names(x$param)[1],
                         ylab=names(x$param)[2], col=col, lwd=2, pch=3,
                         panel=eql.panel, par.settings=list(par.sub.text=list(alpha=1,
                                                              cex=0.75, col=1, font=1,
                                                              lineheight=1)),
                         col.regions=col.regions,
                         sub=paste("n = ", length(x$eql), sep=""))
    args <- c(list(x$eql~x$param[,1]+x$param[,2]), merge(extra.args, default.args))
    return(do.call("levelplot", args))
  } else {
    stop("'plot.eql' is currently defined only for variance families with at most ",
         "two parameters!")
  }
}
