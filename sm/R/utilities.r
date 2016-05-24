"sm.options" <- function (...) {
    if (nargs() == 0) return(.sm.Options)
    current <- .sm.Options
    temp <- list(...)
    if (length(temp) == 1 && is.null(names(temp))) {
        arg <- temp[[1]]
        switch(mode(arg),
               list = temp <- arg,
               character = return(.sm.Options[arg]),
               stop("invalid argument: ", sQuote(arg)))
    }
    if (length(temp) == 0) return(current)
    n <- names(temp)
    if (is.null(n)) stop("options must be given by name")
    changed <- current[n]
    current[n] <- temp
    if (sys.parent() == 0) env <- asNamespace("sm") else env <- parent.frame()
    assign(".sm.Options", current, envir = env)
    invisible(current)
}

# Not sure where this version came from.  It doesn't seem to work.
# "sm.options" <- function (...) {
#     if (nargs() == 0) return(.sm.Options)
#     current <- .sm.Options
#     if (is.character(...))
#         temp <- eval(parse(text = paste(c("list(", ..., ")"))))
#     else temp <- list(...)
#     if (length(temp) == 1 && is.null(names(temp))) {
#         arg <- temp[[1]]
#         switch(mode(arg),
#                list = temp <- arg,
#                character = return(.Options[arg]),
#                stop(paste("invalid argument:", arg)))
#         }
#     if (length(temp) == 0) return(current)
#     n <- names(temp)
#     if (is.null(n)) stop("options must be given by name")
#     changed <- current[n]
#     current[n] <- temp
#     if (sys.parent() == 0) env <- pos.to.env( match(".GlobalEnv", search()) )
#        else env <- parent.frame()
#     assign(".sm.Options", current, envir = env)
#     invisible(current)
#     }


"binning" <- function (x, y, breaks, nbins) {
    binning.1d <- function(x, y, breaks, nbins) {
        f <- cut(x, breaks = breaks)
        if (any(is.na(f)))
            stop("breaks do not span the range of x")
        freq <- tabulate(f, length(levels(f)))
        midpoints <- (breaks[-1] + breaks[-(nbins + 1)])/2
        id <- (freq > 0)
        x <- midpoints[id]
        x.freq <- as.vector(freq[id])
        result <- list(x = x, x.freq = x.freq, table.freq = freq,
            breaks = breaks)
        if (!all(is.na(y))) {
            result$means <- as.vector(tapply(y, f, mean))[id]
            result$sums <- as.vector(tapply(y, f, sum))[id]
            result$devs <- as.vector(tapply(y, f, function(x) 
                              sum((x - mean(x))^2)))[id]
            }
        result
        }
    binning.2d <- function(x, y, breaks, nbins) {
        f1 <- cut(x[, 1], breaks = breaks[, 1])
        f2 <- cut(x[, 2], breaks = breaks[, 2])
        freq <- t(table(f1, f2))
        dimnames(freq) <- NULL
        midpoints <- (breaks[-1, ] + breaks[-(nbins + 1), ])/2
        z1 <- midpoints[, 1]
        z2 <- midpoints[, 2]
        X <- cbind(rep(z1, length(z2)), rep(z2, rep(length(z1), length(z2))))
        X.f <- as.vector(t(freq))
        id <- (X.f > 0)
        X <- X[id, ]
        dimnames(X) <- list(NULL, dimnames(x)[[2]])
        X.f <- X.f[id]
        result <- list(x = X, x.freq = X.f, midpoints = midpoints,
            breaks = breaks, table.freq = freq)
        if (!all(is.na(y))) {
            result$means <- as.numeric(tapply(y, list(f1, f2), mean))[id]
            result$devs <- as.numeric(tapply(y, list(f1, f2),
                function(x) sum((x - mean(x))^2)))[id]
            }
        result
        }
    binning.3d <- function(x, y, breaks, nbins) {
        f1   <- cut(x[, 1], breaks = breaks[, 1])
        f2   <- cut(x[, 2], breaks = breaks[, 2])
        f3   <- cut(x[, 3], breaks = breaks[, 3])
        freq <- table(f1, f2, f3)
        dimnames(freq) <- NULL
        midpoints <- (breaks[-1, ] + breaks[-(nbins + 1), ])/2
        z1  <- midpoints[, 1]
        z2  <- midpoints[, 2]
        z3  <- midpoints[, 3]
        X   <- as.matrix(expand.grid(z1, z2, z3))
        X.f <- as.vector(freq)
        id  <- (X.f > 0)
        X   <- X[id, ]
        dimnames(X) <- list(NULL, dimnames(x)[[2]])
        X.f <- X.f[id]
        result <- list(x = X, x.freq = X.f, midpoints = midpoints,
            breaks = breaks, table.freq = freq)
        if (!all(is.na(y))) {
            result$means <- as.numeric(tapply(y, list(f1, f2, f3), mean))[id]
            result$devs  <- as.numeric(tapply(y, list(f1, f2, f3),
                                        function(x) sum((x - mean(x))^2)))[id]
            }
        result
        }
    if (length(dim(x)) > 0) {
        if (!isMatrix(x))
            stop("wrong parameter x for binning")
        ndim <- dim(x)[2]
        if (ndim > 3)
            stop("binning can be carried out only with 1-3 variables")
        if (missing(y))
            y <- rep(NA, nrow(x))
        if (missing(nbins))
            nbins <- round(log(nrow(x)) / log(2) + 1)
        if (missing(breaks)) {
            breaks <- cbind(seq(min(x[, 1]), max(x[, 1]), length = nbins + 1),
                            seq(min(x[, 2]), max(x[, 2]), length = nbins + 1))
            if (ndim == 3) 
               breaks <- cbind(breaks, seq(min(x[, 3]), max(x[, 3]), length = nbins + 1))
            breaks[1, ] <- breaks[1, ] - rep(10^(-5), ncol(breaks))
            }
        else nbins <- nrow(breaks) - 1
        if (max(abs(breaks)) == Inf | is.na(max(abs(breaks))))
            stop("illegal breaks")
        if (ndim == 2) 
           result <- binning.2d(x, y, breaks = breaks, nbins = nbins)
        else
           result <- binning.3d(x, y, breaks = breaks, nbins = nbins)
        }
    else {
        x <- as.vector(x)
        if (missing(y))
            y <- rep(NA, length(x))
        if (missing(nbins))
            nbins <- round(log(length(x))/log(2) + 1)
        if (missing(breaks)) {
            breaks <- seq(min(x), max(x), length = nbins + 1)
            breaks[1] <- breaks[1] - 10^(-5)
            }
        else nbins <- length(breaks) - 1
        if (max(abs(breaks)) == Inf | is.na(max(abs(breaks))))
            stop("illegal breaks")
        result <- binning.1d(x, y, breaks = breaks, nbins = nbins)
        }
    result
    }


"replace.na" <- function (List, comp, value) {
    arg <- paste(substitute(List), "$", substitute(comp), sep = "")
    arg.value <- eval(parse(text = arg), parent.frame(1))
    if (any(is.na(arg.value))) {
        change <- paste(arg, "<-", deparse(substitute(value)))
        a <- eval(parse(text = change), parent.frame(1))
        }
    invisible()
    }


# "attach.frame" <- function (data, name, ...) {
    # if (missing(name))
        # name <- deparse(substitute(data))
    # if (is.data.frame(data)) {
        # if (!is.na(pos <- match(name, search()))) {
            # cat(paste(name, "already attached, re-attached in 2nd position\n"))
            # detach(pos = pos)
            # }
        # cat(paste("attaching", name, "\n", sep = " "))
        # attach(what = data, pos = 2, name = name, ...)
        # }
    # else {
        # cat(name)
        # cat(" is not a data.frame\n")
        # }
    # invisible()
    # }


"provide.data" <- function (data, path, options = list()) {
	cat("This function is no longer available in the sm package.\n")
	cat("The data and attach functions should be used instead.\n")
}
	
# "provide.data" <- function (data, path, options = list()) {
    # describe <- sm.options(options)$describe
    # name <- deparse(substitute(data))
    # if (missing(path))
        # path <- system.file("smdata", package="sm")
    # datafile <- file.path(path, paste(name, ".dat", sep = ""))
    # docfile <- file.path(path, paste(name, ".doc", sep = ""))
    # if (!exists(name, where=.GlobalEnv, inherits = FALSE)) {
        # if (file.exists(datafile)) {
            # cat("Data file being loaded\n")
            # env <- .GlobalEnv
            # assign(name, read.table(datafile, header = TRUE), envir = env)
            # attach(what = data, name = name)
            # }
        # else cat("Data file does not exist\n")
        # }
    # else {
        # if (!is.data.frame(data))
            # cat("object exists, not as a data.frame\n")
        # else {
            # cat(paste(name, "already loaded\n"))
            # attach.frame(data, name = name)
            # }
        # }
    # if (describe) {
        # if(file.exists(docfile)) file.show(docfile)
        # else cat("Data description file not found\n")
        # }
    # invisible()
    # }


"sm.check.data" <- function (x, y = NA, weights = NA, group = NA, ...) {
   opt <- sm.options(list(...))

   density <- all(is.na(y))
   if (density) X <- x
      else  X <- cbind(x, y)

   if(all(is.na(weights)) | all(weights == 1))
      X <- cbind(X, 1) 
   else{
      if(!is.na(opt$nbins) & opt$nbins!=0) 
         stop("if weights are set, nbins must be either 0 or NA")
      if(any(weights<0 | is.na(weights))) 
         stop("negative or NA weights are meaningless")
      if(any(weights!=round(weights))) {
         weights <- round(weights/min(weights[weights>0]))
         if(opt$verbose>0) 
            cat("Warning: weights have been rescaled to integer values\n")
         }
      X <- cbind(X, weights)
      }

   ndim <- ncol(X) - 1 - (!density)            # dimensionality of x
   if (!all(is.na(group))) {
      X <- cbind(X, group)
      group.col <- ncol(X)
      }
   if (!all(is.na(opt$h.weights))) {
      X <- cbind(X,opt$h.weights)
      hw.col <- ncol(X)
      } 
   if (any(is.na(X)) & opt$verbose > 0) cat("missing data are removed\n")
   X <- na.omit(data.matrix(X))
   if (ndim > 2 + density) stop("x has too many columns")
   weights <- as.vector(X[, ndim + (!density) + 1])
   if (!density) y <- as.vector(X[, ndim + 1])
   x <- if (ndim == 1) as.vector(X[, 1]) else X[, 1:ndim]
   if (!all(is.na(group))) group <- as.vector(X[, group.col])
   if (!all(is.na(opt$h.weights))) opt$h.weights <- X[, hw.col]
   list(x = x, y = y, weights = weights, group = group, ndim = ndim, 
        nobs = nrow(X), density = density, options = opt)
   }

"britmap" <- function () {
    jump <- c(0, sqrt(diff(britpts$britlat)^2 + diff(britpts$britlong)^2))
    flag <- rep(1, nrow(britpts))
    flag[jump >= 0.6] <- NA
    lines(britpts * flag)
    }

"pause" <- function () {
    if(interactive()) readline("Pause. Press <Enter> to continue...")
    invisible()
    }

"wmean" <- function (x, w)
    sum(x * w)/sum(w)

"wvar" <- function (x, w)
    sum(w * (x - wmean(x, w))^2)/(sum(w) - 1)

if(getRversion() >= "2.15.1")
   utils::globalVariables(c("xyzcheck", "llong", "llat", "X", "Y",
        "britlat", "britlong", "theta", "phi", "h.weights", "nbins",
        "hmult", "long2", "lat2", "invislong", "invislat", "smplot",
        "display", "se", "panel.plot", "method", "h.manual", "se.test",
        "smplot", "band", "xgrid", "xlab", "ylab", "xlim", "ylim", 
        "eval.points", "ndim", "delta", "col.band", "col.points", 
        "yht", "pch", "test", "cex", "zlab", "fmat", "ngrid", "zlim"))
