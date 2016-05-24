### Computing risk measures ####################################################

##' @title Value-at-Risk for normal and t distributions
##' @param alpha confidence level
##' @param mu location
##' @param sigma scale
##' @param df degrees of freedom; Inf for the normal distribution
##' @return Value-at-Risk
##' @author Marius Hofert
VaR_t <- function(alpha, mu=0, sigma=1, df=Inf)
{
    stopifnot(0 <= alpha, alpha <= 1, sigma > 0, df > 0)
    mu + sigma * if(identical(df, Inf)) qnorm(alpha) else qt(alpha, df=df)
}

##' @title Expected shortfall for normal and t distributions
##' @param alpha confidence level
##' @param mu location
##' @param sigma scale
##' @param df degrees of freedom; Inf for the normal distribution
##' @return Expected shortfall
##' @author Marius Hofert
ES_t <- function(alpha, mu=0, sigma=1, df=Inf)
{
    stopifnot(0 <= alpha, alpha <= 1, sigma > 0, df > 0)
    mu + (sigma/(1-alpha)) * if(identical(df, Inf)) dnorm(qnorm(alpha)) else
    dt(qt(alpha, df=df), df=df)*(df+qt(alpha, df=df)^2)/(df-1)
}

##' @title Value-at-Risk for the Pareto distribution
##' @param alpha confidence level
##' @param theta Pareto parameter
##' @return Value-at-Risk
##' @author Marius Hofert
##' @note For Pareto defined via F(x)=1-x^{-theta}, F^-(u)=(1-u)^{-1/theta}
VaR_Par <- function(alpha, theta) qPar(alpha, theta=theta)

##' @title Expected shortfall for the Pareto distribution
##' @param alpha confidence level
##' @param theta Pareto parameter
##' @return Expected shortfall
##' @author Marius Hofert
##' @note For Pareto defined via F(x)=1-x^{-theta}, omit the last '-1'
ES_Par <- function(alpha, theta)
{
    stopifnot(0 <= alpha, alpha <= 1, theta > 0)
    (theta/(theta-1))*VaR_Par(alpha, theta=theta) - 1
}


