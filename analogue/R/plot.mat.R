###########################################################################
##                                                                       ##
## plot.mat      - 'plot' method for class 'mat'                         ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## x             - object on which method dispatch applied               ##
## which         - which aspects of 'mat' object to plot if a subset of  ##
##                 the plots is required, specify a subset of the        ##
##                 numbers '1:5'                                         ##
## weighted      - Logical. Should the analysis use weighted mean of env ##
##                 data of analogues as fitted/estimated values?         ##
## k             - number of analogues to use. If missing 'k' is chosen  ##
##                 automatically as the 'k' that achieves lowest RMSE.   ##
## caption       - captions to appear above the plots                    ##
## max.bias      - logical. Should max bias lines be added to residuals  ##
##                 plot?                                                 ##
## n.bias        - number of sections to calculate maximum bias for      ##
## restrict      - logical; restrict comparison of k-closest model to    ##
##                 k <= 'restrict'.                                      ##
## sub.caption   - common title-above figures if there are multiple;     ##
##                 used as 'sub' (s.'title') otherwise.  If 'NULL', as   ##
##                 by default, a possibly shortened version of           ##
##                 'deparse(x$call)' is used.                            ##
## main          - title to each plot-in addition to the above 'caption' ##
## ask           - logical; if 'TRUE', the user is _ask_ed before each   ##
##                 plot, see 'par(ask=.)'.                               ##
## panel         - panel function.  The useful alternative to 'points',  ##
##                 'panel.smooth' can be chosen by 'add.smooth = TRUE'.  ##
## add.smooth    - logical indicating if a smoother should be added to   ##
##                 fitted & residuals plots; see also 'panel' above.     ##
## ...           - arguments passed to other graphics functions          ##
##                                                                       ##
###########################################################################
plot.mat <- function(x,
                     which = c(1:3,5),
                     weighted = FALSE,
                     k,
                     caption = c("Inferred vs Observed", "Residuals vs Fitted",
                       "Leave-one-out errors", "Average bias", "Maximum bias"),
                     max.bias = TRUE,
                     n.bias = 10,
                     restrict = 20,
                     sub.caption = NULL,
                     main = "",
                     ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                     ...,
                     panel = if (add.smooth) panel.smooth else points,
                     add.smooth = getOption("add.smooth")
                     )
  {
    if (!inherits(x, "mat"))
      stop("use only with \"mat\" objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 6))
      stop("'which' must be in 1:6")
    show <- rep(FALSE, 6)
    show[which] <- TRUE
    if(missing(k)) {
      auto <- TRUE
      if(weighted)
        ##k <- which.min(x$weighted$rmse)
        k <- x$weighted$k
      else
        ##k <- which.min(x$standard$rmse)
        k <- x$standard$k
    }
    ## set-up code would go here
    if(any(show[1:2])) {
      Est <- fitted(x, k, weighted = weighted)$estimated
      Obs <- x$orig.y
      Resi <- resid(x, k = k, weighted = weighted)$residuals#[, k]
    }
    if (any(show[3:5])) {
      n.obs <- nrow(x$orig.x)
      if (n.obs > restrict) {
        n.analogs <- restrict
      } else {
        n.analogs <- n.obs
      }
      k.analogs <- deparse(substitute(k <= a, list(a = n.analogs)))
      xlabel <- substitute(paste("No. of analogues (", k <= foo, ")",
          sep = ""), list(foo = n.analogs))
    }
    if (is.null(sub.caption)) {
        cal <- x$call
        if (!is.na(m.f <- match("formula", names(cal)))) {
            cal <- cal[c(1, m.f)]
            names(cal)[2] <- ""
        }
        cc <- deparse(cal, 80)
        nc <- nchar(cc[1])
        abbr <- length(cc) > 1 || nc > 75
        sub.caption <- if (abbr)
            paste(substr(cc[1], 1, min(75, nc)), "...")
        else cc[1]
    }
    ##
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if(show[1]) {
      lims <- range(Est, Obs)
      ylabel <- paste("Fitted (k = ", k, ",weighted = ", weighted, ")", sep = "")
      plot(Obs, Est, type = "n", asp = 1, xlim = lims, ylim = lims,
           ylab = ylabel, xlab = "Observed", ...)
      abline(0, 1, col = "grey", ...)
      panel(Obs, Est, ...)
      if (one.fig)
        title(sub = sub.caption, ...)
      mtext(caption[1], 3, 0.25)
    }
    if(show[2]) {
      ylabel <- paste("Residuals (k = ", k, ", weighted = ", weighted, ")", sep = "")
      plot(Obs, Resi, type = "n", ylab = ylabel, xlab = "Observed", ...)
      abline(h = 0, col = "grey", ...)
      abline(h = mean(Resi), col = "blue", lty = "dashed")
      if(max.bias) {
        groups <- cut(Obs, breaks = n.bias)
        bias <- aggregate(as.vector(Resi), list(group = groups), mean)$x
        ## turn cut intervals into numeric
        interv <- lapply(strsplit(sapply(levels(groups),
                                         function(x) substr(x, 2,
                                                            nchar(x)-1),
                                         USE.NAMES = FALSE), ","),
                         as.numeric)
        ## reformat cut intervals as 2 col matrix for easy plotting
        interv <- matrix(unlist(interv), ncol = 2, byrow = TRUE)
        ## add bias indicators per group
        arrows(interv[,1], bias, interv[,2], bias,
               length = ifelse(one.fig, 0.05, 0.01), angle = 90, code = 3,
               col = "blue")
      }
      panel(Obs, Resi, ...)
      if (one.fig)
        title(sub = sub.caption, ...)
      mtext(caption[2], 3, 0.25)
    }
    if(show[3]) {
      if(weighted) {
        dat <- x$weighted$rmse[1:n.analogs]
      } else {
        dat <- x$standard$rmse[1:n.analogs]
      }
      plot(1:n.analogs, dat, type = "n",
           ylab = paste("RMSEP (weighted = ", weighted, ")", sep = ""),
           xlab = xlabel, ...)
      if(one.fig) {
        lines(1:n.analogs, dat, type = "b", pch = "", ...)
        text(1:n.analogs, dat, labels = as.character(seq(1, n.analogs)),
             cex = 0.7, ...)
      } else {
        lines(1:n.analogs, dat, type = "b", ...)
      }
      if (one.fig)
        title(sub = sub.caption, ...)
      mtext(caption[3], 3, 0.25)
    }
    if(show[4]) {
      if(weighted) {
        dat <- x$weighted$avg.bias[1:n.analogs]
      } else {
        dat <- x$standard$avg.bias[1:n.analogs]
      }
      ## take absolute value of bias, makes plot easier
      dat <- abs(dat)
      plot(1:n.analogs, dat, type = "n",
           ylab = paste("Average bias (weighted = ", weighted, ")", sep = ""),
           xlab = xlabel, ...)
      if(one.fig) {
        lines(1:n.analogs, dat, type = "b", pch = "", ...)
        text(1:n.analogs, dat, labels = as.character(seq(1, n.analogs)),
             cex = 0.7, ...)
      } else {
        lines(1:n.analogs, dat, type = "b", ...)
      }
      if (one.fig)
        title(sub = sub.caption, ...)
      mtext(caption[4], 3, 0.25)
    }
    if(show[5]) {
      if(weighted) {
        dat <- x$weighted$max.bias[1:n.analogs]
      } else {
        dat <- x$standard$max.bias[1:n.analogs]
      }
      ## take absolute value of bias, makes plot easier
      dat <- abs(dat)
      plot(1:n.analogs, dat, type = "n",
           ylab = paste("Maximum bias (weighted = ", weighted, ")", sep = ""),
           xlab = xlabel, ...)
      if(one.fig) {
        lines(1:n.analogs, dat, type = "b", pch = "", ...)
        text(1:n.analogs, dat, labels = as.character(seq(1, n.analogs)),
             cex = 0.7, ...)
      } else {
        lines(1:n.analogs, dat, type = "b", ...)
      }
      if (one.fig)
        title(sub = sub.caption, ...)
      mtext(caption[5], 3, 0.25)
    }
    if (!one.fig && par("oma")[3] >= 1)
        mtext(sub.caption, outer = TRUE, cex = 1.25)
    invisible()
  }
