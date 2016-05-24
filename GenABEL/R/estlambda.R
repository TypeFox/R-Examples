#' Estimate the inflation factor for a distribution of P-values
#'
#' Estimate the inflation factor for a distribution of P-values or 1df
#' chi-square test. The major use of this procedure is the Genomic Control, but
#' can also be used to visualise the distribution of P-values coming from other
#' tests. Methods implemented include 'median' (median(chi2)/0.455...),
#' regression (of observed onto expected) and 'KS' (optimizing the
#' chi2.1df distribution fit by use of Kolmogorov-Smirnov test)
#'
#'
#' @param data A vector of reals. If all are <=1, it is assumed that this is a
#'   vector of P-values, else it is treated as a vector of chi-squares
#' @param plot Whether the plot should be shown or not (default).
#' @param proportion The proportion of lowest P (or
#'   \eqn{\chi^2}{chi^2) values to be used when estimating the inflation
#'   factor \eqn{\lambda}{lambda}.
#' @param method "regression" (default), "median", or "KS": method to
#'   be used for \eqn{\lambda}{lambda} estimation.
#' @param filter if the test statistics with 0-value of
#'   \eqn{\chi^2}{chi^2} should be excluded prior to estimation of
#'   \eqn{\lambda}{lambda}.
#' @param df Number of degrees of freedom.
#' @param ... arguments passed to the \code{\link{plot}} function.
#' @return A list with elements \item{estimate}{Estimate of \eqn{\lambda}{lambda}}
#'   \item{se}{Standard error of the estimate}
#' @author Yurii Aulchenko
#' @seealso \code{\link{ccfast}}, \code{\link{qtscore}}
#' @keywords htest
#' @examples
#'
#' data(srdta)
#' pex <- summary(gtdata(srdta))[,"Pexact"]
#' estlambda(pex, plot=TRUE)
#' estlambda(pex, method="regression", proportion = 0.95)
#' estlambda(pex, method="median")
#' estlambda(pex, method="KS")
#' a <- qtscore(bt,srdta)
#' lambda(a)
#'
"estlambda" <- function(data, plot=FALSE, proportion=1.0,
                        method="regression", filter=TRUE, df=1,... ) {
        data <- data[which(!is.na(data))]
        if (proportion>1.0 || proportion<=0)
                stop("proportion argument should be greater then zero and less than or equal to one")

        ntp <- round( proportion * length(data) )
        if ( ntp<1 ) stop("no valid measurements")
        if ( ntp==1 ) {
                warning(paste("One measurement, lambda = 1 returned"))
                return(list(estimate=1.0, se=999.99))
        }
        if ( ntp<10 ) warning(paste("number of points is too small:", ntp))
        if ( min(data)<0 ) stop("data argument has values <0")
        if ( max(data)<=1 ) {
#		lt16 <- (data < 1.e-16)
#		if (any(lt16)) {
#			warning(paste("Some probabilities < 1.e-16; set to 1.e-16"))
#			data[lt16] <- 1.e-16
#		}
                data <- qchisq(data, 1, lower.tail=FALSE)
        }
        if (filter)
        {
                data[which(abs(data)<1e-8)] <- NA
        }
        data <- sort(data)
        ppoi <- ppoints(data)
        ppoi <- sort(qchisq(ppoi, df=df, lower.tail=FALSE))
        data <- data[1:ntp]
        ppoi <- ppoi[1:ntp]
#	s <- summary(lm(data~offset(ppoi)))$coeff
#       bug fix thanks to Franz Quehenberger

        out <- list()
        if (method=="regression") {
                s <- summary( lm(data~0+ppoi) )$coeff
                out$estimate <- s[1,1]
                out$se <- s[1,2]
        } else if (method=="median") {
                out$estimate <- median(data, na.rm=TRUE)/qchisq(0.5, df)
                out$se <- NA
        } else if (method=="KS") {
                limits <- c(0.5, 100)
                out$estimate <- estLambdaKS(data, limits=limits, df=df)
                if ( abs(out$estimate-limits[1])<1e-4 || abs(out$estimate-limits[2])<1e-4 )
                        warning("using method='KS' lambda too close to limits, use other method")
                out$se <- NA
        } else {
                stop("'method' should be either 'regression' or 'median'!")
        }

        if (plot) {
                lim <- c(0, max(data, ppoi,na.rm=TRUE))
#		plot(ppoi,data,xlim=lim,ylim=lim,xlab="Expected",ylab="Observed", ...)
                oldmargins <- par()$mar
                par(mar=oldmargins + 0.2)
                plot(ppoi, data,
                     xlab=expression("Expected " ~ chi^2),
                     ylab=expression("Observed " ~ chi^2),
                     ...)
                abline(a=0, b=1)
                abline(a=0, b=out$estimate, col="red")
                par(mar=oldmargins)
        }

        out
}
