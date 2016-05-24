"oddsratio" <- function(x, stratum = NULL, log = TRUE)
    loddsratio(x, strata = stratum, log = log)

## "oddsratio" <-
## function (x, stratum = NULL, log = TRUE) {
##   l <- length(dim(x))
##   if (l > 2 && is.null(stratum))
##     stratum <- 3:l
##   if (l - length(stratum) > 2)
##     stop("All but 2 dimensions must be specified as strata.")
##   if (l == 2 && dim(x) != c(2, 2))
##     stop("Not a 2x2 table.")
##   if (!is.null(stratum) && dim(x)[-stratum] != c(2,2))
##     stop("Need strata of 2x2 tables.")

##   lor <- function (y) {
##     if (any(y == 0))
##       y <- y + 0.5
##     y <- log(y)
##     or <- y[1,1] + y[2,2] - y[1,2] - y[2,1]
##     if (log) or else exp(or)
##   }

##   ase <- function(y) {
##     if (any(y == 0))
##         y <- y + 0.5
##     sqrt(sum(1/y))
##   }

##   if(is.null(stratum)) {
##     LOR <- lor(x)
##     ASE <- ase(x)
##   } else {
##     LOR <- apply(x, stratum, lor)
##     ASE <- apply(x, stratum, ase)
##   }

##   structure(LOR,
##             ASE = ASE,
##             log = log,
##             class = "oddsratio"
##             )}

## "print.oddsratio" <-
## function(x, ...) {
##   if (length(dim(x)) > 1)
##     print(cbind(unclass(x)), ...)
##   else
##     print(c(x), ...)
##   invisible(x)
## }

## "summary.oddsratio" <-
## function(object, ...) {
##   if(!is.null(dim(object)))
##     ret <- object
##   else {
##     LOG <- attr(object, "log")
##     ASE <- attr(object, "ASE")
##     Z <- object / ASE

##     ret <- cbind("Estimate"   = object,
##                  "Std. Error" = if (LOG) ASE,
##                  "z value"    = if (LOG) Z,
##                  "Pr(>|z|)"   = if (LOG) 2 * pnorm(-abs(Z))
##                  )
##     colnames(ret)[1] <- if (LOG) "Log Odds Ratio" else "Odds Ratio"
##   }

##   class(ret) <- "summary.oddsratio"
##   ret
## }


## "print.summary.oddsratio" <-
## function(x, ...) {
##   if(!is.null(attr(x, "log"))) {
##     cat("\n")
##     cat(if(attr(x, "log")) "Log Odds Ratio(s):" else "Odds Ratio(s):", "\n\n")
##     print(as.data.frame(unclass(x)), ...)
##     cat("\nAsymptotic Standard Error(s):\n\n")
##     print(attr(x, "ASE"), ...)
##     cat("\n")
##   } else printCoefmat(unclass(x), ...)
##   invisible(x)
## }

## "plot.oddsratio" <-
## function(x,
##          conf_level = 0.95,
##          type = "o",
##          xlab = NULL,
##          ylab = NULL,
##          xlim = NULL,
## 	 ylim = NULL,
##          whiskers = 0.1,
##          baseline = TRUE,
##          transpose = FALSE,
##          ...)
## {
##   if (length(dim(x)) > 1)
##     stop ("Plot function works only on vectors.")

##   LOG <- attr(x, "log")
##   confidence <- !(is.null(conf_level) || conf_level == FALSE)
##   oddsrange <- range(x)

##   if(confidence) {
##     CI  <- confint(x, level = conf_level)
##     lwr <- CI[,1]
##     upr <- CI[,2]
##     oddsrange[1] <- trunc(min(oddsrange[1], min(lwr)))
##     oddsrange[2] <- ceiling(max(oddsrange[2], max(upr)))
##   }

##   if (transpose) {
##     plot(x = unclass(x),
##          y = 1:length(x),
##          ylab = if (!is.null(ylab)) ylab else "Strata",
##          xlab = if (!is.null(xlab)) xlab else if (LOG) "Log Odds Ratio" else "Odds Ratio",
##          type = type,
##          yaxt = "n",
##          xlim = if(is.null(xlim)) oddsrange else xlim,
##          ...)
##     axis (2, at = 1:length(x), names(x))

##     if (baseline)
##       lines(c(1,1) - LOG, c(0,length(x) + 1), lty = 2, col = "red")

##     if (confidence)
##       for (i in 1:length(x)) {
##         lines(c(lwr[i], upr[i]), c(i, i))
##         lines(c(lwr[i], lwr[i]), c(i - whiskers/2, i + whiskers/2))
##         lines(c(upr[i], upr[i]), c(i - whiskers/2, i + whiskers/2))
##       }
##   } else {
##     plot(unclass(x),
##          xlab = if (!is.null(xlab)) xlab else "Strata",
##          ylab = if(!is.null(ylab)) ylab else if(LOG) "Log Odds Ratio" else "Odds Ratio",
##          type = type,
##          xaxt = "n",
##          ylim = if(is.null(ylim)) oddsrange else ylim,
##          ...)
##     axis (1, at = 1:length(x), names(x))

##     if (baseline)
##       lines(c(0,length(x) + 1), c(1,1) - LOG, lty = 2, col = "red")

##     if (confidence)
##       for (i in 1:length(x)) {
##         lines(c(i, i), c(lwr[i], upr[i]))
##         lines(c(i - whiskers/2, i + whiskers/2), c(lwr[i], lwr[i]))
##         lines(c(i - whiskers/2, i + whiskers/2), c(upr[i], upr[i]))
##       }
##   }

## }

## "confint.oddsratio" <-
## function(object, parm, level = 0.95, ...) {
##   ASE <- attr(object, "ASE")
##   LOG <- attr(object, "log")
##   I <- ASE * qnorm((1 + level) / 2)
##   cbind(
##         lwr = if (LOG) object - I else exp(log(object) - I),
##         upr = if (LOG) object + I else exp(log(object) + I)
##         )
## }








