##****************************************************************************
## Plots for Renouv objects.
## The plotting positions used here are taken from the code of the ismev
## package.
##****************************************************************************

PPplot <- function(x, ...) {
  UseMethod("PPplot")
}

QQplot <- function(x, ...) {
  UseMethod("QQplot")
}

##======================================================================
## This must be improved! Historical data should be plotted using colors
## and sylbolds taken from RLpar() and a legend could be drawn.
##======================================================================
PPplot.Renouv <- function(x,
                          showHist = FALSE,
                          legend = FALSE,
                          par = NULL,
                          posOptions = NULL,
                          ...) {
    
  ## if (x$history.MAX$flag || x$history.MAX$flag)  {
  ##   warning("'x' contains historical information which",
  ##           " is ignored at the time. Try 'SandT'.")
  ## }

    if (is.null(par)) rp <- RLpar()
    else rp <- par
    rp1 <- list()
    
    if (!showHist) {
        if (x$nb.OT == 0) {
            stop("No \"OT\" data in 'x'. Try 'showHist = TRUE'")
        }
        samp <- sort(x$y.OT)
        rp1[["col"]] <- rep(rp$OT$col, length(samp))
        rp1[["pch"]] <- rep(rp$OT$pch, length(samp))
        rp1[["bg"]] <- rep(rp$OT$bg, length(samp))
        F.theo <- x$funs$F.y(parm = x$estimate[-1], x = samp)
        F.emp <- (1:x$nb.OT) / x$nb.OT
    } else {
        ST <- SandT(x)
        ## extract color, pch, ...
        Nm <- strsplit(ST$groupName, split = "\\.")
        exFun <- function(x, item) rp[[x]][item]
        for (item in c("col", "pch", "bg")) {
            rpi <- unlist(lapply(Nm, exFun, item))
            rp1[[item]] <- rpi[ST$group] 
        }
        ## defince coords
        F.theo <- x$funs$F.y(parm = x$estimate[-1], x = ST$x -x$threshold)
        F.emp <- 1 - ST$S
    }
    plot(x = F.emp, y = F.theo,
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = "Empirical",  ylab = "Model",
         type = "n", ...)
    abline(h = c(0, 1), col = "gray")
    abline(v = c(0, 1), col = "gray")
    abline(a = 0, b = 1)
    points(x = F.emp, y = F.theo,
           col = rp1[["col"]], pch = rp1[["pch"]], bg = rp1[["bg"]])
    
}

##======================================================================
## This must be improved! Historical data should be plotted using colors
## and sylbolds taken from RLpar() and a legend could be drawn.
##======================================================================
QQplot.Renouv <- function(x,
                          showHist = FALSE,
                          legend = FALSE,
                          par = NULL,
                          posOptions = NULL,
                          ...) {
    
    ## if (x$history.MAX$flag || x$history.MAX$flag)  {
    ##   warning("'x' contains historical information which",
    ##           " is ignored at the time. Try 'SandT'.")
    ## }
    if (is.null(par)) rp <- RLpar()
    else rp <- par
    rp1 <- list()
    
    if (!showHist) {
        if (x$nb.OT == 0) {
            stop("No \"OT\" data in 'x'. Try 'showHist = TRUE'")
        }
        samp <- x$threshold + sort(x$y.OT)
        rp1[["col"]] <- rep(rp$OT$col, length(samp))
        rp1[["pch"]] <- rep(rp$OT$pch, length(samp))
        rp1[["bg"]] <- rep(rp$OT$bg, length(samp))
        p.emp <-  ppoints(x$nb.OT)
        q.theo <- x$threshold + x$funs$q.y(parm = x$estimate[-1], p = p.emp)
    } else {
        ST <- SandT(x)
        ## extract color, pch, ...
        Nm <- strsplit(ST$groupName, split = "\\.")
        exFun <- function(x, item) rp[[x]][item]
        rp1 <- list()
        for (item in c("col", "pch", "bg")) {
            rpi <- unlist(lapply(Nm, exFun, item))
            rp1[[item]] <- rpi[ST$group] 
        }
        p.emp <- 1 - ST$S
        q.theo <- x$threshold + x$funs$q.y(parm = x$estimate[-1], p = p.emp)
        samp <- ST$x
    }
    
    rr <- range(q.theo, samp, na.rm = TRUE)
    
    plot(x = q.theo, y = samp,
         xlim = rr, ylim = rr,
         xlab = "Model",  ylab = "Empirical",
         type = "n", ...)
    abline(a = 0, b = 1)
    points(x = q.theo, y = samp,
           col = rp1[["col"]], pch = rp1[["pch"]], bg = rp1[["bg"]])
    
}
