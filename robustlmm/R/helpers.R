######################################################
## Functions to generate commonly required objects   ##
#######################################################

## std.b: Return the spherical random effects or
##   "Standardize" the Matrix matrix: \eqn{\Lambda^{-1} matrix / \sigma}{Lambda^-1 matrix / sigma}
##
## @title Standardized values
## @param object rlmerMod object
## @param sigma to use for standardization
## @param matrix matrix to standardize
## @param drop apply drop to result?
## @param t transpose result
## @rdname std
std.b <- function(object, sigma = object@pp$sigma, matrix, drop=TRUE, t=FALSE) 
    object@pp$stdB(sigma, matrix, drop, t)

## std.e: Calculate the standardized residuals or
##   "Standardize" the Matrix sigma: \eqn{R^{-1} matrix / \sigma}{R^-1 matrix / sigma}
##
## @rdname std
std.e <- function(object, sigma = object@pp$sigma, matrix, drop=TRUE) {
    if (missing(matrix)) return(object@resp$wtres / sigma)
    ## for the moment: just divide by sigma
    if (drop) matrix <- drop(matrix)
    matrix/sigma
}

## Calculate scaled squared Mahalanobis distances per group
##
## @title Distances per group
## @param object rlmerMod object
## @param sigma to use for standardization
## @param bs spherical random effects
## @param center whether to return the centered distances.
## @param ... ignored
.dk <- function(object, sigma, center, bs = b.s(object), ...) {
    ua <- uArranged(object, bs/sigma)
    unlist(lapply(seq_along(object@blocks), function(bt) {
        us <- ua[[bt]]
        s <- ncol(us)
        if (s == 1) return(us)
        ## else: square, sum and subtract s
        ret <- rowSums(us*us)
        if (center) ret <- ret - object@pp$kappa_b[bt]*s
        ret
        }))
}
## these helper functions are for the non-centered case, i.e.,
## for estimating the random effects themselves.
## modularize distance function: compute sum(bs^2) - s
.d <- function(bs, s=length(bs)) {
    if (s == 1) return(bs)
    if (is.matrix(bs)) sqrt(rowSums(bs*bs)) else sqrt(sum(bs*bs))
}
## same function, but assume we've already summed 
.d2 <- function(sbs2, s) {
    if (s == 1) stop("s must be larger than 1") ## disable this test?
    sqrt(sbs2)
}
## inner derivative of .d2:
.Dd2 <- function(sbs2, s) {
    ##if (s == 1) stop("s must be larger than 1")
    1/.d2(sbs2,s)
}


## dist.b: Calculate the distance from 0 standardized by sigma.
##   This is just value divided by sigma for uncorrelated
##   observations. For correlated items, this is the Mahalanobis
##   distance from 0. If \code{shifted=TRUE} and correlated items,
##   the squared distances are centered by -kappa_b*s. This is
##   required to compute the weights used for the size of the
##   covariance matrix of the random effects.
##
## @title Calculate distance
## @param object object to use
## @param sigma scale for standardization
## @param center whether to use the centered distances
## @param ... passed on to internal functions.
## @rdname dist
dist.b <- function(object, sigma = object@pp$sigma, center=FALSE, ...) {
    db <- .dk(object, sigma, center, ...)
   ## need to take square root if not centering and dim > 1
    if (!center && any(object@dim > 1)) {
        bidx <- object@ind %in% which(object@dim > 1)
        db[bidx] <- sqrt(db[bidx])
    }
    db[object@k]
}

## dist.e: Calculate dist for residuals
##   always assume they are uncorrelated
##
## @rdname dist
dist.e <- function(object, sigma = object@pp$sigma) {
    std.e(object, sigma) ## just the usual rescaled residuals
}

## wgt.b: Calculate the robustness weights psi(d) / d,
##   standardized by sigma. The robustness weights are calculated
##   with d the Mahalanobis distance. Each group of correlated items
##   then gets a constant weight.
##   The robustness weights for the random effects themselves are
##   different than the ones used for estimating the size of the
##   covariance matrix of the random effects. Those are additionally
##   centered. That way, inlier can also be downweighted.
##   If \code{center=TRUE}, then the centered distances are used to
##   compute the robustness weights and the weight function given
##   by rho.sigma.b is used.
##
## @title Calculate robustness weights
## @param object object to use
## @param sigma scale for standardization
## @param center whether return the centered robustness weights, see Details.
## @rdname wgt
## @export
wgt.b <- function(object, sigma = object@pp$sigma, center = FALSE) {
    db <- dist.b(object, sigma, center)
    rho <- rho.b(object, if (center) "sigma" else "default")
    ret <- numeric()
    for (bt in seq_along(object@blocks)) {
        bind <- as.vector(object@idx[[bt]])
        ret <- c(ret, rho[[bt]]@wgt(db[bind]))
    }
    ret
}

## wgt.e: robustness weights of residuals
##
## @param use.rho.sigma return the weights computed using rho.sigma.
## @rdname wgt
## @export
wgt.e <- function(object, sigma = object@pp$sigma, use.rho.sigma = FALSE)
    if (use.rho.sigma) object@rho.sigma.e@wgt(dist.e(object, sigma)) else
       object@rho.e@wgt(dist.e(object, sigma))

### Calculate robustness weights * squared effect
### Return sensible result in the infinite case
## Assume: x infinite <=> y infinite
.wgtxy2 <- function(rho, x, y)
    .wgtxy(rho, x, y*y)
## wgt(x) * y
.wgtxy <- function(rho, x, y) {
    ret <- rho@wgt(x) * y
    ret[is.infinite(x) & is.infinite(y)] <- if (rho@psi(Inf) == 0) 0 else Inf
    ret
}

## Find blocks of correlated random effects
##
## @title Find blocks in Lambda
## @param obj reTrms-object
## @param Lambdat transpose of matrix Lambda (U_b(theta))
## @param Lind vector of indices, mapping theta to Lambdat
## @return list of blocks
findBlocks <- function(obj, Lambdat=obj$Lambdat(), Lind=obj$Lind) {
    LambdaInd <- Lambdat
    LambdaInd@x[] <- as.double(Lind)
    LambdaInd <- t(LambdaInd)
    LambdaInd <- as(LambdaInd, "matrix") ## to avoid attempt to apply non function error
    bend <- unique(apply(LambdaInd != 0, 2, function(x) max(which(x))))
    nblocks <- length(bend)
    bstart <- c(1, bend[-nblocks]+1)
    bidx <- lapply(1:nblocks, function(i) seq.int(bstart[i],bend[i]))
    blocks <- lapply(bidx, function(idx) LambdaInd[idx,idx])
    bind <- match(blocks, ublocks <- unique(blocks))
    k <- unlist(lapply(1:nblocks, function(i) rep(i, length(bidx[[i]]))))
    bdim <- sapply(ublocks, NCOL)
    bidx <- lapply(1:length(ublocks), function (i) 
                   matrix(unlist(bidx[bind == i]),nrow = bdim[i]))
    q <- sapply(bidx, length)
    list(blocks = ublocks, ind = bind, idx = bidx, dim = bdim, q = q, k = k)
}

lchol <- function(x) {
    r <- try(chol.default(x), silent=TRUE)
    ## if chol fails, return sqrt of diagonal
    if (is(r, "try-error")) {
        Diagonal(x = sqrt(diag(x)))
    } else r
}

#######################################################
## Summary / printing methods                        ##
#######################################################

.summary.cor.max <- 20

## Print method
## along the lines of printMerenv of lme4
##' @S3method print summary.rlmerMod
print.summary.rlmerMod <- function(x, digits = max(3, getOption("digits") - 3),
                           correlation = NULL, symbolic.cor = FALSE,
                           signif.stars = getOption("show.signif.stars"),
                           ranef.comp = c("Variance", "Std.Dev."),
                           show.resids = TRUE, ...) {
    ## check if doFit = FALSE is in call
    if (!is.null(x$call$doFit) && !x$call$doFit) {
        cat("Unfitted rlmerMod object. Use update(object, doFit=TRUE) to fit it.\n")
         return(invisible(x))
    }
    ## title
    cat(x$methTitle, "\n")
    ## these are in fromLme4.R
    .prt.call(x$call); cat("\n")
    if (show.resids)
        ## need residuals.merMod() rather than residuals():
        ##  summary.merMod has no residuals method
        .prt.resids(x$residuals, digits=digits)
    .prt.VC(x$varcor, digits=digits, useScale= x$useScale,
	    comp = ranef.comp, ...)
    .prt.grps(x$ngrps, nobs= x$devcomp$dims[["n"]])
    
    ## fixed effecs
    ## this part is 1:1 from printMerenv
    p <- nrow(x$coefficients)
    if (p > 0) {
	cat("\nFixed effects:\n")
	printCoefmat(x$coefficients, zap.ind = 3, #, tst.ind = 4
		     digits = digits, signif.stars = signif.stars)
	if(is.null(correlation)) { # default
	    correlation <- p <= .summary.cor.max
	    if(!correlation) {
		nam <- deparse(substitute(x))
		if(length(nam) > 1 || nchar(nam) >= 32) nam <- "...."
		message(sprintf(paste(
		    "\nCorrelation matrix not shown by default, as p = %d > %d.",
		    "Use print(%s, correlation=TRUE)  or",
		    "	 vcov(%s)	 if you need it\n", sep="\n"),
				p, .summary.cor.max, nam, nam))
	    }
	}
	else if(!is.logical(correlation)) stop("'correlation' must be NULL or logical")
	if(correlation) {
	    if(is.null(VC <- x$vcov)) VC <- vcov(x, correlation=TRUE)
	    corF <- VC@factors$correlation
	    if (is.null(corF)) {
		message("\nCorrelation of fixed effects could have been required in summary()")
		corF <- cov2cor(VC)
	    } ## else {
	    p <- ncol(corF)
	    if (p > 1) {
		rn <- rownames(x$coefficients)
		rns <- abbreviate(rn, minlength=11)
		cat("\nCorrelation of Fixed Effects:\n")
		if (is.logical(symbolic.cor) && symbolic.cor) {
		    corf <- as(corF, "matrix")
		    dimnames(corf) <- list(rns,
					   abbreviate(rn, minlength=1, strict=TRUE))
		    print(symnum(corf))
		} else {
		    corf <- matrix(format(round(corF@x, 3), nsmall = 3),
				   ncol = p,
				   dimnames = list(rns, abbreviate(rn, minlength=6)))
		    corf[!lower.tri(corf)] <- ""
		    print(corf[-1, -p, drop=FALSE], quote = FALSE)
		} ## !symbolic.cor
	    }  ## if (p > 1)
        } ## if (correlation)
    } ## if (p>0)
    ## robustness weights
    summarizeRobWeights(x$wgt.e, digits=3,
                        header="\nRobustness weights for the residuals:")
    summarizeRobWeights(x$wgt.b, digits=3,
                        header="\nRobustness weights for the random effects:")
    ## rho functions
    cat("\nRho functions used for fitting:\n")
    cat("  Residuals:\n")
    cat("    eff:", .sprintPsiFunc(x$rho.e, short=TRUE), "\n")
    cat("    sig:", .sprintPsiFunc(x$rho.sigma.e, short=TRUE), "\n")
    for (bt in seq_along(x$rho.b)) {
        cat("  Random Effects, variance component ", bt, " (", names(x$rho.b)[bt], "):\n", sep="")
        cat("    eff:", .sprintPsiFunc(x$rho.b[[bt]], short=TRUE), "\n")
        cat("    vcp:", .sprintPsiFunc(x$rho.sigma.b[[bt]], short=TRUE), "\n")
    }
    invisible(x)
}

.methTitle <- function(object)
    sprintf("Robust linear mixed model fit by %s", object@method)

##' @S3method print rlmerMod
print.rlmerMod <- function(x, digits = max(3, getOption("digits") - 3),
                           correlation = NULL, symbolic.cor = FALSE,
                           signif.stars = getOption("show.signif.stars"),
                           ranef.comp = "Std.Dev.", ...) {
    ## check if doFit = FALSE is in call
    if (!is.null(x@call$doFit) && !x@call$doFit) {
        cat("Unfitted rlmerMod object. Use update(object, doFit=TRUE) to fit it.\n")
         return(invisible(x))
    }
    dims <- (devC <- x@devcomp)$dims
    ## title
    cat(.methTitle(x), "\n")
    .prt.call(x@call)
    useScale <- as.logical(dims[["useSc"]])

    varcor <- VarCorr(x)
    .prt.VC(varcor, digits=digits, comp = ranef.comp, ...)
    ngrps <- sapply(x@flist, function(x) length(levels(x)))
    .prt.grps(ngrps, nobs= dims[["n"]])
    if(length(cf <- fixef(x)) >= 0) {
	cat("Fixed Effects:\n")
	print.default(format(cf, digits = digits),
		      print.gap = 2L, quote = FALSE, ...)
    } else cat("No fixed effect coefficients\n")
    invisible(x)
}

##' @S3method summary rlmerMod
## this follows the lines of summary.merMod of lme4
summary.rlmerMod <- function(object, ...) {
    resp <- object@resp
    devC <- object@devcomp
    dd <- devC$dims
    cmp <- devC$cmp
    sig <- sigma(object)
    
    coefs <- cbind("Estimate" = fixef(object),
                   "Std. Error" = sig * sqrt(diag(object@pp$unsc())))
    if (nrow(coefs) > 0) {
        coefs <- cbind(coefs, coefs[,1]/coefs[,2], deparse.level=0)
        colnames(coefs)[3] <- "t value"
    }
    varcor <- VarCorr(object)

    structure(list(methTitle=.methTitle(object), devcomp=devC,
                   ngrps=sapply(object@flist, function(x) length(levels(x))),
                   coefficients=coefs, sigma=sig,
                   vcov=vcov(object, correlation=TRUE, sigm=sig),
                   varcor=varcor, # and use formatVC(.) for printing.
                   call=object@call,
                   wgt.e=wgt.e(object),
                   wgt.b=wgt.b(object),
                   rho.e=rho.e(object),
                   rho.sigma.e=rho.e(object, "sigma"),
                   rho.b=rho.b(object),
                   rho.sigma.b=rho.b(object, "sigma"),
                   residuals=residuals(object, scaled=TRUE)
                   ), class = "summary.rlmerMod")
}

##' Use \code{compare} to quickly compare the estaimated parameters of
##' the fits of multiple lmerMod or rlmerMod objects.
##'
##' @title Create comparison charts for multiple fits
##' @param ... objects to compare, or, for the \code{\link{xtable}}
##' functions: passed to the respective \code{\link{xtable}} function.
##' @param digits number of digits to show in output
##' @param dnames names of objects given as arguments (optional)
##' @param show.rho.functions whether to show rho functions in output.
##' @keywords models utilities
##' @examples
##' \dontrun{
##'   fm1 <- lmer(Yield ~ (1|Batch), Dyestuff)
##'   fm2 <- rlmer(Yield ~ (1|Batch), Dyestuff)
##'   compare(fm1, fm2)
##'   require(xtable)
##'   xtable(compare(fm1, fm2))
##'   str(getInfo(fm1))
##' }
##' @export
compare <- function(..., digits = 3, dnames = NULL,
                    show.rho.functions = TRUE) {
    linfos <- list(...)
    if (!missing(dnames) && !is.null(dnames)) names(linfos) <- dnames
    linfos <- lapply(linfos, getInfo)
    ## check if all methods work at least on the same dataset
    if (length(unique(sapply(linfos, function(x) x$data))) > 1)
        warning("Comparison for objects not fitted to the same dataset")
    ## local helper functions
    .NULLtoNA <- function(lst)
        lapply(lst, function(x) if (is.null(x)) NA else x)
    .getComp <- function(lst, comp) 
        sapply(lst, function(x) .NULLtoNA(x[comp]))
    .getComp2 <- function(lst, comp)
        sapply(lst, function(x) x[comp])
    .dropComp <- function(lst, comp)
        lapply(lst, function(x) { x[comp] <- NULL; x })
    .getNames <- function(lst) 
        unique(unlist(lapply(lst, names)))
    .combineComp <- function(lst) {
        names <- .getNames(lst)
        tmp <- .getComp2(lst, names)
        if (!is.matrix(tmp))
            tmp <- matrix(tmp, ncol = length(tmp))
        rownames(tmp) <- names
        format(tmp, digits = digits)
    }
    lnames <- .getNames(linfos)
    ## header
    call <- match.call()
    call$digits <- NULL
    call$names <- NULL
    call$show.rho.functions <- NULL
    split <- rep("", length(linfos))
    ## prepare coef slot:
    ## add stderr
    linfos <- lapply(linfos, function(linfo) {
        cf <- paste(format(linfo[["coef"]], digits=digits), " (",
                    format(linfo[["stderr"]], digits=digits), ")", sep="")
        names(cf) <- names(linfo[["coef"]])
        linfo[["coef"]] <- cf
        linfo[["stderr"]] <- NULL
        linfo
    })
    ## combine
    ## coefficients
    ret <- rbind(Coef=split, 
                 .combineComp(.getComp(linfos, "coef")))
    ## variance components
    ret <- rbind(ret, NULL=split, VarComp=split, 
                 .combineComp(.getComp(linfos, "varcomp")))
    ## correlations if there are any
    if (any(sapply(linfos, function(linfo) !is.null(linfo$correlations))))
        ret <- rbind(ret, NULL=split, Correlations=split,
                     .combineComp(.getComp(linfos, "correlations")))
    ## sigma
    ret <- rbind(ret, NULL=split, format(.getComp(linfos, "sigma"), digits = digits))
    rownames(ret)[nrow(ret)] <- "sigma"
    ## drop the items already included
    linfos <- .dropComp(linfos, c("data", "coef", "varcomp", "correlations", "sigma"))
    ## drop rho functions if requested
    if (!show.rho.functions) 
        linfos <- .dropComp(linfos, grep("^rho", .getNames(linfos), value=TRUE))
    ## show the rest if there is any
    if (length(.getNames(linfos)) > 0) {
        ret <- rbind(ret, NULL=split)
        for (name in .getNames(linfos)) {
            ret <- rbind(ret, format(.getComp(linfos, name), digits = digits))
            rownames(ret)[nrow(ret)] <- name
        }
    }
    ## clean up and finish
    ret <- gsub("\\s*(NA|NULL)", "", ret)
    colnames(ret) <- if (missing(dnames) && is.null(dnames)) as.character(call)[-1] else dnames
    rownames(ret)[rownames(ret) == "NULL"] <- ""
    class(ret) <- "comparison.table"
    ret
}

##' @S3method print comparison.table
print.comparison.table <- function(x, ...) {
    class(x) <- "matrix"
    print(x, ..., quote=FALSE)
}

##' @rdname compare
##' @method getInfo lmerMod
##' @S3method getInfo lmerMod
getInfo.lmerMod <- function(object, ...) {
    lsum <- summary(object)
    coefs <- lsum$coefficients
    varcor <- lsum$varcor
    isREML <- .isREML(object)
    .namedVector <- function(mat) {
        if (is.vector(mat)) return(mat)
        names <- rownames(mat)
        ret <- drop(mat)
        names(ret) <- names
        ret
    }
    .getVC <- function(varcor) {
        vc <- lapply(varcor, function(grp) attr(grp, "stddev"))
        ret <- unlist(vc, use.names = FALSE)
        names(ret) <-
            unlist(lapply(1:length(vc), function(i)
                          paste(names(vc[[i]]), names(vc)[i], sep=" | ")))
        ##.namedVector(ret)
        ret
    }
    .getCorr <- function(varcor) {
        ret <- lapply(1:length(varcor), function(i) {
            grp <- varcor[[i]]
            corr <- attr(grp, "correlation")
            if (nrow(corr) == 1) return(NULL)
            names <- outer(colnames(corr), paste("x", rownames(corr)), paste)
            ret <- as.vector(corr[upper.tri(corr)])
            names(ret) <- paste(as.vector(names[upper.tri(names)]),
                                names(varcor)[i], sep = " | ")
            ret
        })
        unlist(ret)
    }
    ret <- list(data = object@call$data,
                coef = .namedVector(coefs[,1,drop=FALSE]),
                stderr = .namedVector(coefs[,2,drop=FALSE]),
                varcomp = .getVC(varcor),
                sigma = sigma(object))
    corrs <- .getCorr(varcor)
    if (length(corrs) > 0) ret$correlations <- corrs
    if (!is(object, "rlmerMod")) {
        if (isREML) {
            ret$REML <- lme4::REMLcrit(object)
        } else {
            ret$deviance <- deviance(object, REML=FALSE)
        }
    }
    ret
}

##' @rdname compare
##' @method getInfo rlmerMod
##' @S3method getInfo rlmerMod
getInfo.rlmerMod <- function(object, ...) {
    linfo <- getInfo.lmerMod(object)
    linfo$REML <- linfo$deviance <- NULL
    linfo$rho.e <- .sprintPsiFunc(rho.e(object), TRUE)
    linfo$rho.sigma.e <- .sprintPsiFunc(rho.e(object, "sigma"), TRUE)
    rho.b <- rho.b(object)
    rho.sigma.b <- rho.b(object, "sigma")
    for (bt in seq_along(object@blocks)) {
        linfo[[paste("rho.b",bt,sep="_")]] <- .sprintPsiFunc(rho.b[[bt]], TRUE)
        linfo[[paste("rho.sigma.b",bt,sep="_")]] <- .sprintPsiFunc(rho.sigma.b[[bt]], TRUE)
    }
    linfo
}

##' The functions \code{xtable.comparison.table} and
##' \code{print.xtable.comparison.table} are wrapper functions for the
##' respective \code{\link{xtable}} and \code{\link{print.xtable}}
##' functions.
##' 
##' @rdname compare
##' @param x object of class "comparison.table" or "xtable.comparison.table"
##' @param caption see \code{\link{xtable}}.
##' @param label see \code{\link{xtable}}.
##' @param align see \code{\link{xtable}}.
##' @param display see \code{\link{xtable}}.
##' @seealso \code{\link{xtable}}
##' @examples
##' require(xtable)
##' xtable(compare(fm1, fm2))
##' @importFrom xtable xtable
##' @export
##' @method xtable comparison.table
xtable.comparison.table <- function(x, caption=NULL, label=NULL, align=NULL,
                                    digits=NULL, display=NULL, ...) {
    rn <- sapply(rownames(x), function(n) {
        switch(n,
               Coef="Coefficients (Std. Error)",
               VarComp="Variance components",
               Correlations="Correlations",
               n) })
    tbl <- cbind(rn, x)
    rownames(tbl) <- NULL
    colnames(tbl) <- c(" ", colnames(x))
    if (is.null(align)) align <- c("r", "r", rep.int("l", ncol(x)))
    xtbl <- xtable(tbl, caption=caption, label=label, align=align,
                   digits=digits, display=display, ...)
    class(xtbl) <- c("xtable.comparison.table", class(xtbl))
    xtbl
}

##' @rdname compare
##' @param add.hlines replace empty lines in comparison table by hlines.
##'   Supersedes \code{hline.after} argument of \code{print.xtable}.
##' @param latexify.namescol replace \dQuote{sigma} and \dQuote{x} in
##'   the first column by latex equivalents.
##' @param include.rownames include row numbers (the object returned by
##'  \code{xtable.comparison.table} includes names in the first column)
##' @importFrom xtable print.xtable
##' @seealso \code{\link{print.xtable}}
##' @method print xtable.comparison.table
##' @export
print.xtable.comparison.table <- function(x, add.hlines=TRUE,
                                          latexify.namescol=TRUE,
                                          include.rownames=FALSE, ...) {
    args <- list(...)
    if (add.hlines) {
        rns <- if (is.factor(x[[1]])) levels(x[[1]])[x[[1]]] else x[[1]]
        emptyCol <- sapply(rns, nchar) == 0
        x <- x[!emptyCol,]
        args$hline.after=c(-1, 0, which(emptyCol)-1:sum(emptyCol),nrow(x))
    }
    if (latexify.namescol) {
        if (!is.null(args$sanitize.text.function)) {
            sanitize <- args$sanitize.text.function
        } else {
            ## from print.xtable
            sanitize <- function(str) {
                result <- str
                result <- gsub("\\\\","SANITIZE.BACKSLASH",result)
                result <- gsub("$","\\$",result,fixed=TRUE)
                result <- gsub(">","$>$",result,fixed=TRUE)
                result <- gsub("<","$<$",result,fixed=TRUE)
                result <- gsub("|","$|$",result,fixed=TRUE)
                result <- gsub("{","\\{",result,fixed=TRUE)
                result <- gsub("}","\\}",result,fixed=TRUE)
                result <- gsub("%","\\%",result,fixed=TRUE)
                result <- gsub("&","\\&",result,fixed=TRUE)
                result <- gsub("_","\\_",result,fixed=TRUE)
                result <- gsub("#","\\#",result,fixed=TRUE)
                result <- gsub("^","\\verb|^|",result,fixed=TRUE)
                result <- gsub("~","\\~{}",result,fixed=TRUE)
                result <- gsub("SANITIZE.BACKSLASH","$\\backslash$",result,fixed=TRUE)
                return(result)
                  }
        }
        ## sanitize text
        for (i in 1:ncol(x)) x[[i]] <- sanitize(x[[i]])
        args$sanitize.text.function <- identity
        x[[1]] <- sapply(x[[1]], function(rn) {
            rn <- sub(" x ", " $\\times$ ", rn, fixed=TRUE)
            rn <- sub("\\bsigma\\b", "$\\\\sigma$", rn)
            rn
        })
    }
    args$include.rownames <- include.rownames
    args$x <- x
    do.call(print.xtable, args)
}
 
##' @S3method update rlmerMod
update.rlmerMod <- function(object, formula., ..., evaluate = TRUE) {
    ## update call
    ## set old object as init object
    ## run it
    if (is.null(call <- object@call))
        stop("object should contain a 'call' component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.)) {
        call$formula <- update.formula(formula(object), formula.)
        extras$init <- NULL
    } else {
        ## set init to object, if not explicitly given (and no new data given)
        if (is.null(extras[["data"]])) {
            if (is.null(extras[["init"]])) {
                lcall <- sys.call(sys.parent())
                extras$init <- object
            }
            ## copy pp and resp (to really get a new object)
            extras$init@pp <- object@pp$copy()
            ## reset calledInit... fields to FALSE:
            fields <- grep("calledInit", names(getRefClass(class(extras$init@pp))$fields()), value=TRUE)
            Map(function(field) extras$init@pp$field(field, FALSE), fields)
            extras$init@resp <- object@resp$copy()
        } else {
            extras$init <- NULL
        }
    }
    if (length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate)
        eval(call, parent.frame())
    else call
}

#######################################################
## predict method                                    ##
#######################################################

##' @importFrom stats predict
##' @S3method predict rlmerMod
predict.rlmerMod <- function(object, ...) {
    class(object) <- "lmerMod"
    predict(object, ...)
}
