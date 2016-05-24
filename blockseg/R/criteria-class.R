##' Penalized criteria based on estimation of degrees of freedom
##'
##' Produce a plot or send back the values of some penalized criteria
##' accompanied with the vector(s) of parameters selected
##' accordingly. The default behavior plots the BIC and the AIC (with
##' respective factor \eqn{\log(n)}{log(n)} and \eqn{2}{2}) yet the user can specify any
##' penalty.
##'
##' @param object output of a fitting procedure of the \pkg{blockseg}
##' package (e.g. \code{\link{blockSeg}}). Must be of class
##' \code{blockSeg}.
##' @param Y matrix of observations.
##' @param penalty a vector with as many penalties a desired. The
##' default contains the penalty corresponding to the AIC and the BIC
##' (\eqn{2}{2} and \eqn{\log(n)}{log(n)}). Setting the "names"
##' attribute, as done in the default definition, leads to outputs
##' which are easier to read.
##' @param sigma scalar: an estimate of the residual variance. When
##' available, it is plugged-in the criteria, which may be more
##' relevant. If \code{NULL} (the default), it is estimated as usual
##' (see details).
##' @param log.scale logical; indicates if a log-scale should be used
##' when \code{xvar="lambda"}. Default is \code{TRUE}.
##' @param xvar variable to plot on the X-axis: either \code{"df"}
##' (the estimated degrees of freedom), \code{"lambda"}
##' (\eqn{\lambda_1}{lambda1} penalty level) or \code{"fraction"}
##' (\eqn{\ell_1}{l1}-norm of the coefficients). Default is set to
##' \code{"lambda"}.
##' @param plot logical; indicates if the graph should be plotted on
##' call. Default is \code{TRUE}.
##'
##' @return When \code{plot} is set to \code{TRUE}, an invisible
##' \pkg{ggplot2} object is returned, which can be plotted via the
##' \code{print} method. On the other hand, a list with a two data
##' frames containing the criteria and the chosen vector of parameters
##' are returned.
##'
##' @seealso \code{\linkS4class{blockSeg}}.
##'
##' @note When \code{sigma} is provided, the criterion takes the form
##'
##' \if{latex}{\deqn{\left\|\mathbf{y} - \mathbf{X} \hat{\beta} \right\|^2 +
##' \mathrm{penalty} \times \frac{\hat{\mathrm{df}}}{n} \ \sigma^2.}}
##' \if{html}{\out{ <center> RSS + penalty * df / n * sigma<sup>2</sup> </center>}}
##' \if{text}{\deqn{RSS + penalty * df / n * sigma^2}}
##'
##' When it is unknown, it writes
##'
##' \if{latex}{\deqn{\log\left(\left\|\mathbf{y} - \mathbf{X} \hat{\beta} \right\|^2\right) +
##' \mathrm{penalty} \times \hat{\mathrm{df}}.}}
##' \if{html}{\out{ <center> n*log(RSS) + penalty * df </center>}}
##' \if{text}{\deqn{n*log(RSS) + penalty * df}}
##'
##' Estimation of the degrees of freedom (for the elastic-net, the
##' LASSO and also bounded regression) are computed by applying and
##' adpating the results of Tibshirani and Taylor (see references
##' below).
##'
##' @references Ryan Tibshirani and Jonathan Taylor. Degrees of
##' freedom in lasso problems, Annals of Statistics, 40(2) 2012.
##'
##' @rdname criteria
##'
##' @examples 
##' n <- 100
##' K <- 5
##' mu <- suppressWarnings(matrix(rep(c(1,0),ceiling(K**2/2)), K,K))
##' Y <- rblockdata(n,mu,sigma=.5)$Y
##' res <- blockSeg(Y, 50)
##' criteria(res, Y, sigma=.5)
##'
##' @exportMethod criteria
setGeneric("criteria", function(object, Y, penalty=setNames(c(2, log(length(Y))), c("AIC","BIC")), sigma=NULL,
            log.scale=TRUE, xvar = "lambda", plot=TRUE)
           {standardGeneric("criteria")})

##' @rdname criteria
setMethod("criteria", "blockSeg", definition =
   function(object, Y, penalty=setNames(c(2, log(length(Y))), c("AIC","BIC")), sigma=NULL,
            log.scale=TRUE, xvar = "lambda", plot=TRUE) {
       
     Y.hat <- getCompressYhat(object, Y)
     lambda <- object@Lambda
     N <- length(Y)

     ## Compute generalized cross-validation
     GCV <- deviance(object,Y)/(N*(1+getComplexity(object)/N))^2

     ## compute the penalized criteria
     if (is.null(sigma)) {
         crit <- sapply(penalty, function(pen) log(deviance(object,Y)/N) + pen * getComplexity(object)/N)
     } else {
         crit <- sapply(penalty, function(pen) deviance(object,Y)/N + pen * getComplexity(object)/N * sigma^2)
     }

     ## put together all relevant information about those criteria
     criterion <- data.frame(crit, GCV=GCV, df=getComplexity(object), lambda=lambda, row.names=1:nrow(crit))

     ## recover the associated vectors of parameters
     Y.hat.min  <- setNames(Y.hat[apply(crit, 2, which.min)], names(penalty))
     lambda.min <- setNames(lambda[apply(crit, 2, which.min)], names(penalty))
     ## plot the critera, if required
     if (plot) {

       data.plot <- melt(criterion, id=xvar, measure=1:length(penalty), variable.name="criterion", value.name="value")
       rownames(data.plot) <- 1:nrow(data.plot)

       colnames(data.plot)[1] <- "xvar"

       xlab <- switch(xvar,
                      "df" = "Estimated degrees of freedom",
                      ifelse(log.scale,expression(log[10](lambda)),expression(lambda))
                      )

       d <- ggplot(data.plot, aes_string(x="xvar", y="value", colour="criterion", group="criterion")) +
           geom_line(aes_string(x="xvar",y="value")) + geom_point(aes_string(x="xvar",y="value")) +
               labs(x=xlab, y="criterion's value",  title=paste("Information Criteria"))
       
       if (log.scale & xvar=="lambda") {
           d <- d + scale_x_log10()
       }
       print(d)
       return(invisible(d))
   } else {
       return(list(criterion=criterion, Y.hat.min=Y.hat.min, lambda.min=lambda.min))
   }

 })
