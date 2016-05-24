GOF <- function(x, 
    lik=function(theta,x) { 
        return(sum ( dnorm(x,mean=theta[1],sd=theta[2],log=T) ) )
    } ,
    pdist=function(x,theta) { 
        return(pnorm(x,mean=theta[1],sd=theta[2]) ) 
    } , 
    p=c(0,1),
    cutpts=quantile(x), 
    paramNames = paste('parameter',1:length(p)),
    pearson=FALSE,
    ...) 
{
    nlmaxResult <- nlmax(lik,p=p,x=x,...)
    mle <- nlmaxResult$estimate
    names(mle) <- paramNames
    prob <- diff( pdist(cutpts, mle))
    n <- length(x)
    o <- table(cut(x,cutpts))
    e <- prob * n


    pearsonStat <- sum( (o-e)^2 / e)
    lrtStat <- 2 * sum( o * log (o/e) )
    df = length(cutpts) - 2 - 2

    if (pearson) {
        pval <- 1- pchisq(pearsonStat,df=df)
        method= "Pearson Goodness of Fit Test"
        stat = pearsonStat
    } else {
        pval <- 1- pchisq(lrtStat,df=df)
        method= "LRT Goodness of Fit Test"
        stat = lrtStat
    }

    names(df) = 'df'
    names(stat) = "X-squared"

    if (nlmaxResult$code > 2) { 
        warning("MLE estimates have not converged sufficiently.") 
    }

    structure(list(
        nlmax = nlmaxResult,
        statistic = stat, 
        estimate = mle,
        parameter = df, 
        p.value = pval, 
        method = method,
        data.name = deparse(substitute(x)), 
        observed = o, 
        expected = e, 
        residuals = (o - e)/sqrt(e),
        table = cbind(o,e,prob)
        ), 
        class = "htest")
}
data <- c(18.0,6.3,7.5,8.1,3.1,0.8,2.4,3.5,9.5,39.7,
          3.4,14.6,5.1,6.8,2.6,8.0,8.5,3.7,21.2,3.1,
          10.2, 8.3,6.4,3.0,5.7,5.6,7.4,3.9,9.1,4.0)
GOF(data, cutpts=c(0,3,5,9,13,Inf),iterlim=1000,p=c(5,5))$table
GOF(data, cutpts=c(0,3,5,9,13,Inf),iterlim=1000,p=c(5,5))
GOF(data, cutpts=c(0,3,5,9,13,Inf),iterlim=1000,p=c(5,5),pearson=T)
