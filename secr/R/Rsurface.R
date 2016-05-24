Rsurface <- function (mask,  sigma, usecov = NULL, alpha2 = 1, detectfn = "HHN", z = 1,
                      inverse = FALSE, scale = TRUE) {
    if (ms(mask)) stop ("not ready for multisession masks")
    mm <- nrow(mask)
    detectfn <- valid.detectfn(detectfn, c(4,14:18))  ## converts from character
    tmpmask <- cbind(mask, rep(1,mm))
    miscparm <- c(1,0,0,0)
    if (!is.null(usecov)) {
        miscparm[2] <- 1
        tmpmask <- cbind(tmpmask, exp(alpha2 * covariates(mask)[,usecov]))
    }
    miscparm[3] <- scale
    temp <- .C ( "getdenomext",
                as.integer(detectfn),
                as.double (miscparm),
                as.double(unlist(tmpmask)),
                as.integer(mm),
                as.double (sigma),
                as.double (z),
                invdenom = double(mm),
                scale = double(1))
    if (is.null(covariates(mask)))
        covariates(mask) <- data.frame(matrix(nrow = mm, ncol = 0))
    covariates(mask)[,"Resource"] <-
        if (inverse) temp$invdenom
        else 1/temp$invdenom
    OK <- is.finite(covariates(mask)$Resource)
    covariates(mask)$Resource[!OK] <- NA
    class(mask) <- c('Rsurface', 'mask', 'data.frame')  ## need data.frame to guide ms()
    attr(mask, 'scale') <- temp$scale
    mask
}
############################################################################################

Rsurface.as.data.frame <- function (x) {
    covnames <- names(covariates(x))
    OK <- match('Resource', covnames) ## just one col for now
    covnames <- covnames[OK]
    resources <- covariates(x)[,covnames]
    df <- cbind(x, resources)
    names(df) <- c('x','y',covnames)
    df
}
############################################################################################

print.Rsurface <- function (x, ...) {
#    if (ms(x)) {   ## no need yet for ms()
#        out <- vector('list')
#        for (session in names(x)) {
#            cat ('Session ', session, '\n')
#            print(x[[session]], ...)
#            out[[session]] <- x[[session]]
#        }
#        names(out) <- names(x)
#        out
#    }
#    else {
        df <- Rsurface.as.data.frame(x)
        print(df, ...)
#    }
    invisible(df)
}
############################################################################################

plot.Rsurface <- function (x, covariate = 'Resource', plottype = 'shaded',
     scale = 1, ...) {
    if (ms(x)) {
        breaklist <- lapply(x, plot, covariate, plottype, ...)
        invisible(breaklist)
    }
    else {
        if (length(covariate)>1)
            stop ("whoa... just one at a time")
        if (!(covariate %in% names(covariates(x))))
            stop ("covariate ", covariate, " not found")
        covariates(x)[,covariate] <- covariates(x)[,covariate] * scale
        if (plottype %in% c('contour','persp')) {
            xval <- sort(unique(x$x))
            yval <- sort(unique(x$y))
            if (nrow(x) != length(xval)*length(yval)) {
                x <- rectangularMask(x)
                if(nrow(x) != length(xval)*length(yval))
                    stop ("failed to convert irregular mask to rectangle")
            }
            zmat <- matrix(covariates(x)[,covariate], nrow = length(xval))
            if (plottype == 'contour')
                contour(x=xval, y=yval, z=zmat, ...)
            else
                persp(x=xval, y=yval, z=zmat, ...)
        }
        else {
            class(x) <- c('mask','data.frame')
            covlevels <- plot(x, covariate = covariate, dots = (plottype == 'dots'), ...)
            if (!is.null(covlevels)) invisible(covlevels)
        }
    }
}
############################################################################################
