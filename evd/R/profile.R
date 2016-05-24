
"profile.evd" <-  function(fitted, which = names(fitted$estimate), conf = 0.999, mesh = fitted$std.err[which]/4, xmin = rep(-Inf, length(which)), xmax = rep(Inf, length(which)), convergence = FALSE, method = "BFGS", control = list(maxit = 500), ...)
{
    if(!inherits(fitted, "evd")) 
        stop("Use only with `evd' objects")
    if(inherits(fitted, "extreme"))
        stop("profiles not implemented for this model")
    if(length(xmin) != length(which))
        stop("`xmin' and `which' must have the same length")
    if(length(xmax) != length(which))
        stop("`xmax' and `which' must have the same length")
    if(length(fitted$estimate) < 2)
        stop("cannot profile one dimensional likelihood")
    if(!is.character(which))
        stop("`which' must be a character vector")
    if(!all(which %in% names(fitted$estimate)))
        stop("`which' contains unrecognized or unestimated parameters")
    if(is.null(fitted$std.err) && missing(mesh))
       stop("fitted model must contain standard errors")
    prof.list <- as.list(numeric(length(which)))
    names(xmin) <- names(xmax) <- names(prof.list) <- which
    if(is.null(names(mesh))) names(mesh) <- which
    mles <- fitted$estimate[which]                   
    for(j in which) {
        print(paste("profiling",j))
        prof1 <- prof2 <- matrix(nrow = 0, ncol = length(fitted$estimate) + 1)
        npmles <- fitted$estimate[!names(fitted$estimate) %in% j]
        start <- as.list(npmles)
        call.args <- c(list(fitted$data, start, 0), as.list(fitted$fixed),
           list(FALSE, FALSE, method, FALSE, control))
        names(call.args) <- c("x", "start", j, names(fitted$fixed),
           "std.err", "corr", "method", "warn.inf", "control")
        dimnames(prof1) <- dimnames(prof2) <- list(NULL, c(j, "deviance",
          names(start)))
        call.fn <- paste("f", class(fitted)[1], sep="")
        if(inherits(fitted, "gev")) {
            call.args$nsloc <- fitted$nsloc
            call.args$prob <- fitted$prob
        }
        if(inherits(fitted, "pot")) {
            call.args$threshold <- fitted$threshold
            call.args$npp <- fitted$npp
            call.args$period <- fitted$period
            call.args$cmax <- fitted$cmax
            call.args$r <- fitted$r
            call.args$ulow <- fitted$ulow
            call.args$rlow <- fitted$rlow
            call.args$mper <- fitted$mper
        }
        if(inherits(fitted, "bvevd")) {
            call.args$model <- fitted$model
            call.args$nsloc1 <- fitted$nsloc1
            call.args$nsloc2 <- fitted$nsloc2
            call.args$sym <- fitted$sym
            call.args$cloc <- fitted$cmar[1]
            call.args$cscale <- fitted$cmar[2]
            call.args$cshape <- fitted$cmar[3] 
        }
        if(inherits(fitted, "bvpot")) {
            call.args$threshold <- fitted$threshold
            call.args$likelihood <- fitted$likelihood
            call.args$model <- fitted$model
            call.args$sym <- fitted$sym
            call.args$cscale <- fitted$cmar[1]
            call.args$cshape <- fitted$cmar[2]          
        }
        lcnt <- TRUE; ppar <- mles[j]
        while(lcnt) {
            ppar <- as.vector(ppar + mesh[j])
            if(ppar >= xmax[j]) ppar <- as.vector(xmax[j])      
            call.args[[j]] <- ppar
            fit.mod <- do.call(call.fn, call.args)
            if(convergence) print(fit.mod$convergence)
            call.args[["start"]] <- as.list(fit.mod$estimate)
            rop <- c(ppar, fit.mod$deviance, fit.mod$estimate)
            prof1 <- rbind(prof1, rop)
            ddf <- fit.mod$deviance - fitted$deviance
            lcnt <- (ddf <= qchisq(conf, 1)) && (ppar != xmax[j])
        }
        call.args[["start"]] <- as.list(npmles)
        lcnt <- TRUE; ppar <- mles[j]
        while(lcnt) {
            ppar <- as.vector(ppar - mesh[j])
            if(ppar <= xmin[j]) ppar <- as.vector(xmin[j])
            call.args[[j]] <- ppar
            fit.mod <- do.call(call.fn, call.args)
            if(convergence) print(fit.mod$convergence)
            call.args[["start"]] <- as.list(fit.mod$estimate)
            rop <- c(ppar, fit.mod$deviance, fit.mod$estimate)
            prof2 <- rbind(prof2, rop)
            ddf <- fit.mod$deviance - fitted$deviance
            lcnt <- (ddf <= qchisq(conf, 1)) && (ppar != xmin[j])
        }
        rop <- c(mles[j], fitted$deviance, npmles)
        prof2 <- prof2[nrow(prof2):1, ,drop = FALSE]
        prof <- rbind(prof2, rop, prof1)
        rownames(prof) <- NULL
        
        rdev <- qchisq(conf, 1) + fitted$deviance
        if(prof[1, "deviance"] == 2e6)  {
          prof <- prof[-1, ,drop = FALSE]
          if(prof[1,"deviance"] <= rdev)
            warning(paste("If", j, "is to satisfy `conf',",
              "`mesh' must be smaller"))
        }
        if(prof[nrow(prof), "deviance"] == 2e6) {
          prof <- prof[-nrow(prof), ,drop = FALSE]
          if(prof[nrow(prof),"deviance"] <= rdev)
            warning(paste("If", j, "is to satisfy `conf',",
              "`mesh' must be smaller"))
        }
        prof.list[[j]] <- prof
    }
    structure(prof.list, deviance = fitted$deviance,
              xmin = xmin, xmax = xmax, class = "profile.evd")
}
    
profile2d <- function (fitted, ...) {
    UseMethod("profile2d")
}

"profile2d.evd" <-  function(fitted, prof, which, pts = 20, convergence = FALSE, method = "Nelder-Mead", control = list(maxit = 5000), ...)
{
    if(!inherits(fitted, "evd")) 
        stop("Use only with `evd' objects")
    if(inherits(fitted, "extreme"))
        stop("profiles not implemented for this model")
    if (!inherits(prof, "profile.evd")) 
        stop("`prof' must be a `profile.evd' object")
    if(length(fitted$estimate) < 3)
        stop("Cannot profile two dimensional likelihood")
    if(missing(which) || !is.character(which) || length(which) != 2)
        stop("`which' must be a character vector of length two")
    if(!all(which %in% names(fitted$estimate)))
        stop("`which' contains unrecognized or unestimated parameters")
    if(!all(which %in% names(prof)))
        stop("`which' contains unprofiled parameters")
    if(is.null(fitted$std.err))
       stop("fitted model must contain standard errors")
    prof.list <- as.list(numeric(3))
    names(prof.list) <- c("trace", which)
    limits1 <- range(prof[[which[1]]][,1])
    limits2 <- range(prof[[which[2]]][,1])
    mles <- fitted$estimate[which]                   
    prof <- matrix(NA, nrow = pts^2, ncol = length(fitted$estimate) + 1)
    parvec1 <- seq(limits1[1], limits1[2], length = pts)
    prof.list[[which[1]]] <- parvec1
    parvec2 <- seq(limits2[1], limits2[2], length = pts)
    prof.list[[which[2]]] <- parvec2
    pars <- expand.grid(parvec1, parvec2)
    start <- as.list(fitted$estimate[!names(fitted$estimate) %in% which])
    # if method unspecified supress optim 1d warnings
    if(missing(method) && length(start) == 1) oldopt <- options(warn = -1)
    call.args <- c(list(fitted$data, start, 0, 0), as.list(fitted$fixed),
       list(FALSE, FALSE, method, FALSE, control))
    names(call.args) <- c("x", "start", which[1], which[2],
       names(fitted$fixed), "std.err", "corr", "method",
       "warn.inf", "control")
    dimnames(prof) <- list(NULL, c(which, "deviance", names(start)))
    call.fn <- paste("f", class(fitted)[1], sep="")
    if(inherits(fitted, "gev")) {
        call.args$nsloc <- fitted$nsloc
        call.args$prob <- fitted$prob
    }
    #if(inherits(fitted, "pot")) {
    #    call.args$threshold <- fitted$threshold
    #    call.args$npp <- fitted$npp
    #    call.args$period <- fitted$period
    #    call.args$cmax <- fitted$cmax
    #    call.args$r <- fitted$r
    #    call.args$ulow <- fitted$ulow
    #    call.args$rlow <- fitted$rlow
    #    call.args$mper <- fitted$mper
    #}
    if(inherits(fitted, "bvevd")) {
        call.args$nsloc1 <- fitted$nsloc1
        call.args$nsloc2 <- fitted$nsloc2
        call.args$model <- fitted$model
        call.args$sym <- fitted$sym
        call.args$cloc <- fitted$cmar[1]
        call.args$cscale <- fitted$cmar[2]
        call.args$cshape <- fitted$cmar[3]
    }
    if(inherits(fitted, "bvpot")) {
        call.args$threshold <- fitted$threshold
        call.args$likelihood <- fitted$likelihood
        call.args$model <- fitted$model
        call.args$sym <- fitted$sym
        call.args$cscale <- fitted$cmar[1]
        call.args$cshape <- fitted$cmar[2]          
    }
    for(i in 1:pts^2) {
        call.args[[which[1]]] <- pars[i,1]
        call.args[[which[2]]] <- pars[i,2]
        fit.mod <- do.call(call.fn, call.args)
        if(convergence) print(fit.mod$convergence)
        prof[i,1] <- pars[i,1]
        prof[i,2] <- pars[i,2]
        prof[i,3] <- fit.mod$deviance
        prof[i,-(1:3)] <- fit.mod$estimate
    }
    prof.list[["trace"]] <- prof
    if(missing(method) && length(start) == 1) oldopt <- options(oldopt)
    if(any(prof[,"deviance"] == 2e6))
        warning("non-convergence present in profile2d object")
    structure(prof.list, deviance = fitted$deviance, class = "profile2d.evd")
}

"plot.profile.evd" <-  function(x, which = names(x), main = NULL,
     ask = nb.fig < length(which) && dev.interactive(), ci = 0.95,
     clty = 2, ...) 
{
    if (!inherits(x, "profile.evd")) 
        stop("Use only with `profile.evd' objects")
    if(!is.character(which))
        stop("`which' must be a character vector")
    if(!all(which %in% names(x)))
        stop("`which' contains unprofiled parameters")
    nb.fig <- prod(par("mfcol"))
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if(is.null(main)) {
        fls <- toupper(substr(which, 1, 1))
        ols <- substr(which, 2, nchar(which))
        cwhich <- paste(fls, ols, sep = "")
        main <- paste("Profile Log-likelihood of", cwhich)
    }
    for(i in which) {
        plot(spline(x[[i]][,1], -x[[i]][,2]/2, n = 75), type = "l",
             xlab = i, ylab = "profile log-likelihood",
             main = main[match(i,which)], ...)
        cdist <- -(attributes(x)$deviance + qchisq(ci, df = 1))/2
        abline(h = cdist, lty = clty)
    }
    invisible(x)
}

"confint.profile.evd" <- function(object, parm, level = 0.95, ...)
{
    if(missing(parm)) 
      parm <- names(object)
    if(!all(parm %in% names(object)))
      stop("`parm' contains unprofiled parameters")
    rdev <- attributes(object)$deviance + qchisq(level, df = 1)
    pct <- c("lower", "upper")
    ci <- array(NA, dim = c(length(parm), 2), dimnames = list(parm, pct))
    # Assumes profile trace is unimodal
    for(i in parm) {
      x <- object[[i]]
      n <- nrow(x)
      th.l <- (x[1, 1] == attributes(object)$xmin[i])
      th.u <- (x[n, 1] == attributes(object)$xmax[i])
      halves <- c(diff(x[,"deviance"]) < 0, FALSE)
      if(x[1,"deviance"] <= rdev && !th.l) {
        warning(paste("cannot calculate lower confidence limit for", i))
        ci[i,1] <- NA
      }
      if(x[1,"deviance"] <= rdev && th.l) ci[i,1] <- x[1, 1]
      if(x[1,"deviance"] > rdev)
        ci[i,1] <- approx(x[halves,2], x[halves,1], xout = rdev)$y
      if(x[n,"deviance"] <= rdev && !th.u) {
        warning(paste("cannot calculate upper confidence limit for", i))
        ci[i,2] <- NA
      }
      if(x[n,"deviance"] <= rdev && th.u) ci[i,2] <- x[n, 1]
      if(x[n,"deviance"] > rdev)
        ci[i,2] <- approx(x[!halves,2], x[!halves,1], xout = rdev)$y 
    }
    ci
}

"plot.profile2d.evd" <-  function(x, main = NULL, ci = c(0.5,0.8,0.9,0.95,0.975, 0.99, 0.995), col = heat.colors(8), intpts = 75, xaxs = "r", yaxs = "r", ...)
{
    if (!inherits(x, "profile2d.evd")) 
        stop("Use only with `profile2d.evd' objects")
    which <- names(x)[2:3]
    if(is.null(main)) {
        fls <- toupper(substr(which, 1, 1))
        ols <- substr(which, 2, nchar(which))
        cwhich <- paste(fls, ols, sep = "")
        main <- paste("Profile Log-likelihood of", cwhich[1], "and", cwhich[2])
    }
    br.pts <- attributes(x)$deviance + qchisq(c(0,ci), df = 2)
    prof <- x$trace[,"deviance"]
    if(any(prof == 2e6))
        warning("non-convergence present in profile2d object")
    prof <- -prof/2
    br.pts <- (-br.pts/2)[length(br.pts):1]
    col <- col[length(col):1]
    
    if(!requireNamespace("akima", quietly = TRUE)) {
        image(x[[which[1]]], x[[which[2]]],
              matrix(prof, nrow = length(x[[which[1]]])),
              col = col, breaks = c(-1e6+1, br.pts),
              main = main, xlab = which[1], ylab = which[2], xaxs = xaxs,
              yaxs = yaxs, ...)
    }
    else {
        lim1 <- range(x[[which[1]]])
        lim2 <- range(x[[which[2]]]) 
        prof.interp <- akima::interp(x$trace[,1], x$trace[,2], prof,
            xo = seq(lim1[1], lim1[2], length = intpts),
            yo = seq(lim2[1], lim2[2], length = intpts))
        image(prof.interp, col = col, breaks = c(min(prof), br.pts),
              main = main, xlab = which[1], ylab = which[2], xaxs = xaxs,
              yaxs = yaxs, ...)
    }
    invisible(x)
}









