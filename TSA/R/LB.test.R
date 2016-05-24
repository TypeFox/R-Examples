`LB.test` <-
function (model, lag = 12, type = c("Ljung-Box","Box-Pierce"),no.error=FALSE,
omit.initial=TRUE) 
{
# programmed by Kung-Sik Chan 
# Date: March 28, 2006
#
# Modified from the Box.test function.
#
# input: model=model fit from the arima function
#        lag=number of lags of the autocorrelation of the residuals to be
#            included in the test statistic. (default=12)
#        type=either Ljung-Box or Box-Pierce 
#        omit.initial=if true, (d+Ds) initial residuals are omitted from the test
# Output: a list containing the following elements
#         statistic = test statistic
#         p.value = p-value
#         parameter = d.f. of the Chi-square test
#         lag = no of lags.
#
x=residuals(model)
d1=length(model$mod$Delta)
if(omit.initial) x=window(x,start=time(x)[d1+1])
narma=sum(eval(model$arma)[1:4])
if(is.null(model$call$fixed)) nparm=narma else nparm=sum(is.na(eval(model$call$fixed)[1:narma]))
if(lag<=nparm) { if(no.error) return(list(p.value=NA)) else stop("number of lags is less than the number of parameters: increase the lag")}    
if (NCOL(x) > 1) 
        stop("x is not a vector or univariate time series")
    DNAME <- paste("residuals from ",deparse(substitute(model)))
    type <- match.arg(type)
    cor <- acf(x, lag.max = lag, plot = FALSE, na.action = na.pass,drop.lag.0=TRUE)
    n <- sum(!is.na(x))
    PARAMETER <- lag-nparm
    obs <- cor$acf[1:lag]
    if (type == "Box-Pierce") {
        METHOD <- "Box-Pierce test"
        STATISTIC <- n * sum(obs^2)
        PVAL <- 1 - pchisq(STATISTIC, PARAMETER)
    }
    else {
        METHOD <- "Box-Ljung test"
        STATISTIC <- n * (n + 2) * sum(1/seq(n - 1, n - lag) * 
            obs^2)
        PVAL <- 1 - pchisq(STATISTIC, PARAMETER)
    }
    names(STATISTIC) <- "X-squared"
    names(PARAMETER) <- "df"
    names(lag)<-"lag"
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD, data.name = DNAME,lag=lag), 
        class = "htest")
}

