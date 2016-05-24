##' Function to draw meta-analysis forest plots
##'
##' This function creates forest plots from meta-analysis data.
##'
##' @param estimate Vector of effect estimates
##' @param se Vector of standard errors
##' @param labels Vector of labels for the individual studies
##'           (default: Study 1, Study 2, etc.)
##' @param CI Confidence interval (default: 0.95)
##' @param xexp Whether the effect values are to be depicted on an
##'           exponential scale (default: \code{FALSE})
##' @param xlab Label for the horizontal axis (default: \eqn{\beta})
##' @param ylab Label for the horizontal axis (default: empty)
##' @param ... Arguments passed to the \code{plot} function,
##'           e.g. \code{main="My plot"}
##' @author Yurii Aulchenko, Lennart C. Karssen
##' @keywords hplot
##' @export
##' @examples
##' beta <- c(0.16, 0.091, 0.072, -0.03)
##' se   <- c(0.07, 0.042, 0.048, 0.12)
##' forestplot(beta, se, main="Example plot")
"forestplot" <-
    function(estimate, se,
             labels=paste("Study", c(1:length(estimate))),
             CI=0.95, xexp=FALSE, xlab=expression(beta), ylab="", ...) {
        hoff      <- 3
        del       <- 10
        mea       <- !is.na(estimate)
        estimate  <- estimate[mea]
        se        <- se[mea]
        labels    <- labels[mea]
        w2        <- 1. / (se * se)
        invsumw2  <- 1. / sum(w2)
        mestimate <- sum(estimate * w2) * invsumw2
        mse       <- sqrt(invsumw2)
        npop      <- length(estimate)
        estimate[npop+1] <- mestimate
        se[npop+1]       <- mse
        labels[npop+1]   <- "Pooled"
        chi2 <- round(estimate * estimate / (se * se), 2)
        p <- sprintf("%5.1e", pchisq(estimate * estimate / (se * se),
                                     1, lower.tail=FALSE))
        ##	p[as.numeric(p)<0] <- "<1.e-16"

        if (CI > 1 || CI < 0) {
            stop("CI argument should be between 0 and 1")
        }

        cimultip <- qnorm(1 - (1 - CI) / 2)
        lower    <- estimate - cimultip * se
        upper    <- estimate + cimultip * se

        if (xexp) {
            estimate <- exp(estimate)
            lower    <- exp(lower)
            upper    <- exp(upper)
        }

        cntr <- 0; if (xexp) cntr <- 1;
        lbnd <- (-.1); if (xexp) lbnd <- 0.9
        rbnd <- (.1); if (xexp) rbnd <- 1.1
        minv <- min(lower)
        minv <- minv - abs(minv / 10)
        minv <- min(lbnd, minv)
        maxv <- max(upper)
        maxv <- maxv + abs(maxv / 10)
        maxv <- max(rbnd, maxv)
        hgt  <- (length(estimate) + 1) * del

        if (any(is.na(estimate))) stop("estimate contains NAs")
        if (any(is.na(se))) stop("se contains NAs")

        plot(x=c(cntr, cntr), y=c(0, hgt), xlim=c(minv, maxv),
             ylim=c(0, hgt), type="l", lwd=2, lty=2,
             xlab=xlab, ylab=ylab, yaxt='n', ...)

        ## Draw the bars for the individual studies
        for (i in c(1:(length(estimate)-1))) {
            points(x=c(lower[i], upper[i]), y=c((i) * del, (i) * del),
                   type="l", lwd=2)
            points(x=c(estimate[i]), y=c((i) * del), pch=19, cex=1)

            labeltext <- bquote(
                .(labels[i]) ~ "(" * chi^2 ~ "=" ~ .(chi2[i]) * ","
                ~ italic(P) ~ "=" ~ .(p[i]) * ")"
                )
            text(estimate[i], i * del + 1, labeltext, pos=3, cex=.7)
        }

        ## Draw diamond of the estimate
        for (i in c(length(estimate))) {
            points(x=c(lower[i], estimate[i]),
                   y=c((i) * del, (i) * del + hoff),
                   type="l", lwd=2)
            points(x=c(estimate[i], upper[i]),
                   y=c((i) * del + hoff, (i) * del),
                   type="l", lwd=2)
            points(x=c(upper[i], estimate[i]),
                   y=c((i) * del, (i) * del - hoff),
                   type="l", lwd=2)
            points(x=c(lower[i], estimate[i]),
                   y=c((i) * del, (i) * del - hoff),
                   type="l", lwd=2)

            labeltext <-bquote(
                .(labels[i]) ~ "(" * chi^2 ~ "=" ~ .(chi2[i]) * ","
                ~ italic(P) ~ "=" ~ .(p[i]) * ")"
                )
            text(estimate[i], i * del + 5, labeltext, pos=3, cex=1)
        }
    }
