# The code is adapted from package 'sem' by J. Fox
# Uses: resempls.R
bootsempls <- function(object, nboot=200, start=c("ones", "old"),
                method=c("ConstructLevelChanges", "IndividualSignChanges", "Standard"),
                verbose=TRUE, strata, ...){
    method <- match.arg(method)
    if(method=="IndividualSignChanges"){
      stop("Not yet implemented.\n",
           "Try 'ConstructLevelChanges' or 'Standard'.")
    }
    refit <- function(){
        data <- data[indices,]
        refitted_model <- resempls(object, data, start, method)
        refitted_model
    }
    if (!require("boot")) stop("package boot not available")
    # the following 2 lines borrowed from boot in package boot
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    warn <- options(warn=-2)
    on.exit(options(warn)) # insure restore even in event of error
    nErrors <- 0
    data <- object$data
    N <- nrow(data)
    if(missing(strata)) strata <- rep(1, N)
    coefficients <- object$coefficients[,2]
    coef_names <- rownames(object$coefficients)
    coefs <- matrix(numeric(0), nrow=nboot, ncol=length(coefficients))
    attr(coefs, "path") <- object$coefficients[,1]
    colnames(coefs) <- coef_names
    outer_weights <- matrix(NA, nrow=nboot, ncol=ncol(data))
    colnames(outer_weights) <- rownames(object$outer_weights)
    clcIndices <- NULL
    tryErrorIndices <- NULL
    bootIndices <- matrix(NA, nrow=nboot, ncol=nrow(data))
    if (verbose) cat("Resample: ")
    for (b in 1:nboot){
      if (verbose){
        treil <- paste(rep(" ", floor(log10(nboot)) - floor(log10(b))), collapse="")
        ndel <- paste(rep("\b", floor(log10(nboot)) + 1), collapse="")
        if(b==1) cat(paste(treil, b, sep=""))
        if(b!=nboot) cat(paste(ndel, treil, b, sep=""))
        else cat(paste(ndel, b, " Done.\n", sep=""))
      }
      for (try in 1:11){
        if (try > 10) stop("more than 10 consecutive convergence failures")
        indices <- sample(N, N, replace=TRUE)
        res <- try(refit(), silent=TRUE)
        if (inherits(res, "try-error")){
          nErrors <- nErrors + 1
          tryErrorIndices <-  rbind(tryErrorIndices, indices)
        }
        else {
          # for construct level changes
          if(method=="ConstructLevelChanges"){
            clcIndices[[b]] <- res$clcIndex
          }
          bootIndices[b,] <- indices
          # Standard
          coefs[b,] <- res$coefficients[,2]
          outer_weights[b,] <- res$outer_weights[res$outer_weights!=0]
          break()
        }
      }
    }
    options(warn)
    if (nErrors > 0) warning("There were ", nErrors,
                             " apparent convergence failures;\n",
                             "  these are discarded from the ",
                              nboot, " bootstrap replications returned.")
    res <- list(t0=coefficients, t=coefs, nboot=nboot, data=data, seed=seed,
                statistic=refit, sim="ordinary", stype="i", call=match.call(),
                tryErrorIndices=tryErrorIndices,  clcIndices= clcIndices,
                bootIndices=bootIndices, outer_weights=outer_weights,
                fitted_model=object, strata=strata)
    class(res) <- c("bootsempls", "boot")
    res
}


print.bootsempls <- function(x, digits = 3, ...){
    t <- x$t
    t0 <- x$t0
    result <- data.frame("Estimate"=t0, "Bias"=colMeans(t, ...) - t0,
        "Std.Error"=apply(t, 2, sd, ...))
    rownames(result) <- attr(t, "path")
    cat("Call: ")
    dput(x$call)
    cat("\n")
    print(result, digits=digits, ...)
    invisible(x)
}


summary.bootsempls <- function(object,
    type=c("perc", "bca", "norm", "basic", "none"), level=0.95, ...){
    if ((!require("boot")) && (type != "none")) stop("boot package unavailable")
    type <- match.arg(type)
    t <- object$t
    t0 <- object$t0
    object$R <- object$nboot
    result <- data.frame("Estimate"=t0, "Bias"=colMeans(t) - t0,
        "Std.Error"=apply(t, 2, sd))
    if (type != "none"){
        p <- length(t0)
        lower <- upper <- rep(0, p)
        low <- if (type == "norm") 2 else 4
        up  <- if (type == "norm") 3 else 5
        noCi <- NULL
        for (i in 1:p){
          ti <- t[!is.na(t[,i]),i] # 14.10.2010
          #if (boot:::const(t[,i], min(1e-08, mean(t[,i], na.rm = TRUE)/1e+06))){
          if (boot:::const(ti, min(1e-08, mean(ti, na.rm = TRUE)/1e+06))){
            lower[i] <- upper[i] <- NA
            noCi <- append(noCi, i)
          }
          else{
            ci <- try(as.vector(boot.ci(object, type=type, index=i,
                      conf=level)[[type, exact=FALSE]]))
            if(inherits(ci, "try-error") && type=="bca"){
                stop("Try to set 'nboot' to the number of observations!\n")
            }
            lower[i] <- ci[low]
            upper[i] <- ci[up]
            }
        }
        result$Lower <- lower
        result$Upper <- upper
        }
    rownames(result) <- colnames(t)
    attr(result, "path") <- attr(t, "path")
    result <- list(table=result, call=object$call, level=level, type=type)
    class(result) <- "summary.bootsempls"
    result
}


print.summary.bootsempls <- function(x, na.print=".", digits = 3, ...){
    cat("Call: ")
    dput(x$call)
    cat("\n")
    if (x$type != "none") {
        cat(paste("Lower and upper limits are for the", 100*x$level,
            "percent", x$type, "confidence interval\n\n"))
        }
    xChar <- format(x$table, digits=digits, ...)
    xChar[is.na(x$table)] <- na.print
    print(xChar, ...)
    invisible(return(x))
}

