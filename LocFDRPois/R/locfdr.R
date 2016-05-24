#' @title Build a table corresponding to the histogram of data
#'
#' @param data The raw data (before histogramming) to estimate a mixture
#'   density for.
#'
#' @return X A data frame whose two columns are integer and number of times that
#'    integer appears in the data.
BuildTable <- function(data) {
    N <- length(data)
    m <- max(data)

    # Build data frame with counts for each 0, ... , m number of mutations
    X0 <- data.frame(data=0:m)
    X <- data.frame(table(data))
    X <- merge(X0, X, all=T)
    X[is.na(X$Freq), "Freq"] <- 0
    X$data <- as.numeric(X$data)
    return (X)
}

#' @title Compute full mixture density for data, assuming exponential family model
#'
#' @param data The raw data (before histogramming) to estimate a mixture density
#'   for.
#' @param df The number of degrees of freedom to use in the natural spline
#'   in fitting the mixture density.
#' @param knots The positions of the knots to use in the natural spline fit.
#'
#' @return f\_hat The estimated mixture density function (can be called on numeric
#'    vectors.)
#'
#' @importFrom splines ns
MixtureDensity <- function(data, df, knots) {
    X <- BuildTable(data)
    model <- glm(Freq ~ splines::ns(data, df, knots), family="poisson", data=X)

    f_hat <- function(x) {
        if (any(x > max(data))) {
            stop("Predicting on data outside of training range")
        }
        newdf <- data.frame(data=x)
        return (as.numeric(predict(model, newdf, type="response")) / length(data))
    }
    return (f_hat)
}

#' @title Return log-likelihood function giving Poisson version of equation (4.6)
#' in "Microarrays, Empirical Bayes, and the Two-Groups model" by Efron.
#'
#' @param r0 The number of points equal to 0
#' @param r1 The number of points equal to 1
#' @param N0 The number of points equal to either 0 or 1
#' @param N The total number of points
#'
#' @return LL A function of lambda0 and pi0, given the current data log likelihood
LLConstructor <- function(r0, r1, N0, N) {
    LL <- function(params) {
        lambda0 <- params[1]
        pi0 <- params[2]
        P0 <- sum(dpois(c(0, 1), lambda0))
        theta <- pi0 * P0

        ll <- N0 * log(theta) + (N - N0) * log(1 - theta) +
            r0 * dpois(0, lambda0, log=T) + r1 * dpois(1, lambda0, log=T) - N0 * log(P0)
        return (ll)
    }
    return (LL)
}

#' @title Optimize lambda0 and pi0, using Poisson version of analytical approach
#' in Efron's "Microarrays, Empirical Bayes, and the Two-Groups Model"
#'
#' In the paper, he derives formulas for arbitrary support. We use support
#' [0, 1] for estimating the null. We initialize the search with lambda=0.05 and
#'  pi=0.95.
#'
#' @param r0 The number of 0's in the data.
#' @param r1 The number of 1's in the data.
#' @param N The total number of points in the data.
AnalyticalOptim <- function(r0, r1, N) {
    N0 <- r0 + r1
    ll <- LLConstructor(r0, r1, N0, N)
    nll <- function(z) { -ll(z) }
    fitted <- suppressWarnings(optim(c(0.05, .95), fn=nll))
    return (list(lambda=fitted$par[1], pi0=fitted$par[2]))
}

#' @title Estimate the density for the null data
#'
#' @param data The raw data (before histogramming) to estimate a mixture density
#'   for.
#'
#' @return f_hat The estimated mixture density function (can be called on numeric
#'    vectors.)
#'
#' @importFrom dplyr filter
NullDensity <- function(data) {
    N <- length(data)
    X <- BuildTable(data)
    r0 <- dplyr::filter(X, data==0)$Freq
    r1 <- dplyr::filter(X, data==1)$Freq

    analytical_optim <- AnalyticalOptim(r0, r1, N)
    lambda0 <- analytical_optim$lambda
    pi0 <- analytical_optim$pi0
    f0 <- function(x) { dpois(x, lambda0) }
    return (list(f0=f0, pi0=pi0, lambda0=lambda0))
}

#' @title Run local fdr functions
#'
#' @param data The raw data (before histogramming) to estimate a mixture density
#'   for.
#' @param df The number of degrees of freedom to use in the natural spline
#'   in fitting the mixture density.
#' @param knots The positions of the knots to use in the natural spline fit.
#'
#' @return fdr The estimated local fdr function (can be called on numeric vectors).
#' @return pi0 The estimated proportion of data that lies in the poisson component
#' @return lambda0 The estimated parameter for the poisson component.
#' @return f0 The estimated null density function (can be called on numeric
#'    vectors.)
#' @return f The estimated mixture density function (can be called on numeric
#'    vectors.)
LocfdrFuns <- function(data, df, knots) {
    f <- MixtureDensity(data, df, knots)
    null_dens <- NullDensity(data)
    f0 <- null_dens$f0
    pi0 <- null_dens$pi0
    fdr <- function(z) {
        pmin(1, pi0 * f0(z) / f(z))
    }
    return (list(fdr=fdr, pi0=pi0, lambda0=null_dens$lambda0, f0=f0, f=f))
}

#' @title Call and summarize output to LocfdrFuns
#'
#' @param x The raw data (before histogramming) to estimate a mixture density for.
#' @param df The number of degrees of freedom to use in the natural spline
#'   in fitting the mixture density.
#' @param knots The positions of the knots to use in the natural spline fit.
#'
#' @return pi0 The estimated proportion of data that lies in the poisson component
#' @return lambda0 The estimated parameter for the poisson component.
#' @return locfdr\_res The fdr, f0 and f functions evaluated on the support of x
#' @return locfdr\_fig A plot of the estimated locfdr fit
#'
#' @importFrom ggplot2 ggplot geom_bar aes geom_line scale_y_sqrt scale_fill_gradient2
#' @export
SummarizeLocfdr <- function(x, df = 5, knots = c()) {
    # Run locfdr
    sim_funs <- LocfdrFuns(x, df, knots)

    # Evaluate estimated functions on support of data
    f0 <- sim_funs$f0(0:max(x))
    f <- sim_funs$f(0:max(x))

    # Calculate local fdr
    N <- length(x)
    fdr <- sim_funs$fdr(0:max(x))

    # Combine and plot
    locfdr_res <- data.frame(BuildTable(x), f0=N * f0, f=N * f, fdr=fdr)
    locfdr_fig <- ggplot2::ggplot(locfdr_res) +
        ggplot2::geom_bar(ggplot2::aes_string(x="data", y="Freq", fill="fdr"), stat="identity") +
        ggplot2::geom_line(ggplot2::aes_string(x="data", y="f0"), col="midnightblue") +
        ggplot2::geom_line(ggplot2::aes_string(x="data", y="f"), col="plum3") +
        ggplot2::scale_y_sqrt() +
        ggplot2::scale_fill_gradient2(high="steelblue", low="indianred", midpoint=0.2)
    return(list(pi0=sim_funs$pi0, lambda0=sim_funs$lambda0, locfdr_res=locfdr_res,
                locfdr_fig=locfdr_fig))
}
