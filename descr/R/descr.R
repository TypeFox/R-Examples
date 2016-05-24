
.onLoad <- function(libname, pkgname) {
    library.dynam("descr", pkgname, libname, local = FALSE);

    if(is.null(getOption("descr.plot")))
        options(descr.plot = TRUE)
    if(is.null(getOption("descr.warn")))
        options(descr.warn = TRUE)
    if(is.null(getOption("descr.na.replacement")))
        options(descr.na.replacement = "NA")

}

.onUnload <- function(libpath) {
    library.dynam.unload("descr", libpath)
}


# R does not have variable labels.
descr <- function (x)
{
    if (class(x)[1] == "data.frame") {
	l <- length(x)
	bnames <- names(x)
	for (i in 1:l) {
	    lb <- attr(x[[i]], "label")
	    if (length(lb) > 0) {
		cat("\n", bnames[i], " - ", lb, "\n", sep = "")
	    }
	    else {
		cat("\n", bnames[i], "\n", sep = "")
	    }
	    print(summary(x[[i]]))
	}
	return(invisible(NULL))
    }
    else {
	lb <- attr(x, "label")
	if (length(lb) > 0) {
	    cat(deparse(substitute(x)), " - ", lb, "\n", sep = "")
	}
	print(summary(x))
	return(invisible(NULL))
    }
}


# The original versions of the functions freq, hist.kdnc, and LogRegR2 were
# written by Dirk Enzmann <dirk.enzmann@jura.uni-hamburg.de> who has given me
# permission to include them in this package. The original code can be found at
# http://www2.jura.uni-hamburg.de/instkrim/kriminologie/Mitarbeiter/Enzmann/Software/Enzmann_Software.html


# Plot histogram of variable with kernel density estimates and normal curve:
# I had to change the name because the "." was causing R to think that the
# function was a method of hist.
histkdnc <- function (v, breaks = 0, include.lowest = TRUE, right = TRUE,
    main = "Histogram with kernel density and normal curve",
    xlab = deparse(substitute(v)), col = grey(0.90),
    col.cur = c("red", "blue"), lty.cur = c(1, 1),
    xlim = NULL, ylim = NULL, ...) 
{
    v2 <- na.omit(v)
    x <- v2
    h <- hist.default(v2, plot = FALSE)
    if (length(breaks) == 1) 
	breaks <- h$breaks
    dens <- density(v2)
    argv <- list(...)
    if(is.null(ylim))
        ylim <- range(0, h$density, dnorm(x = v2, mean = mean(v2), sd = sd(v2)),
                      dens$y)
    if(is.null(xlim))
        xlim <- range(v2, dens$x)
    hist(v2, freq = FALSE, breaks = breaks, include.lowest = include.lowest, 
	right = right, xlim = xlim, ylim = ylim, col = col, 
	xlab = xlab, main = main, ...)
    lines(density(v2), col = col.cur[1], lty = lty.cur[1])
    curve(dnorm(x, mean = mean(v2), sd = sd(v2)), col = col.cur[2],
          add = TRUE, lty = lty.cur[2])
}


# Print multiple R2 analogs
print.LogRegR2 <- function(x, ...)
{
    cat(formatC(gettext("Chi2", domain = "R-descr"), flag = "-", width = 20), x$Chi2, "\n")
    cat(formatC(gettext("Df", domain = "R-descr"), flag = "-", width = 20), x$df, "\n")
    cat(formatC(gettext("Sig.", domain = "R-descr"), flag = "-", width = 20), x$p, "\n")
    cat(formatC(gettext("Cox and Snell Index", domain = "R-descr"), flag = "-", width = 20), x$CoxR2, "\n")
    cat(formatC(gettext("Nagelkerke Index", domain = "R-descr"), flag = "-", width = 20), x$NagelkerkeR2, "\n")
    cat(formatC(gettext("McFadden's R2", domain = "R-descr"), flag = "-", width = 20), x$RL2, "\n")
    return(invisible(NULL))
}

# Calculates multiple R2 analogs (pseudo R2) of logistic regression:
LogRegR2 <- function(model)
{
    if (!(model$family$family == "binomial" && (model$family$link == "logit" || model$family$link == "probit")))
	stop("No logistic regression model, no pseudo R^2 computed.")

    n    <- dim(model$model)[1]
    Chi2 <- model$null - model$dev
    Df   <- model$df.null - model$df.res
    p    <- 1-pchisq(Chi2,Df)

    Cox  <- 1-exp(-Chi2/n)             # Cox & Snell Index
    Nag  <- Cox/(1-exp(-model$null/n)) # Nagelkerke Index
    RL2  <- Chi2/model$null            # also called McFaddens R2

    x <- list('Chi2'=Chi2,'df'=Df,'p'=p,'RL2'=RL2,'CoxR2'=Cox,'NagelkerkeR2'=Nag)
    class(x) <- "LogRegR2"
    x
}

no.drop.levels <- function(x)
{
    if(sum(is.na(x)) > 0){
        nl <- length(levels(x)) + 1
        lv <- c(levels(x), options("descr.na.replacement"))
        io <- is.ordered(x)
        x <- as.numeric(x)
        x[is.na(x)] <- nl
        x <- factor(x, levels = 1:nl, labels = lv, ordered = io)
    }
    x
}

