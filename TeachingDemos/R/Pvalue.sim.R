Pvalue.norm.sim <- function(n=50, mu=0, mu0=0, sigma=1, sigma0=sigma,
                            test=c('z','t'),
                            alternative=c('two.sided', 'less', 'greater',
                            '<>','!=','<','>'),
                            alpha=0.05, B=1e4) {
    test <- match.arg(test)
    alternative <- match.arg(alternative)

    x <- matrix(rnorm( n*B, mu, sigma ), nrow=n)
    xbar <- colMeans(x)

    if( is.na(sigma0) ) sigma0 <- apply(x, 2, sd)

    ts <- (xbar - mu0)/sigma0*sqrt(n)

    pdist <- switch(test,
               z=function(x, lower.tail) pnorm(x, lower.tail=lower.tail),
               t=function(x, lower.tail) pt(x, df=n-1, lower.tail=lower.tail)
                    )

    p.vals <- switch(alternative,
                     '!='=,'<>'=,
                     two.sided = 2*pmin( pdist(ts,TRUE), pdist(ts,FALSE) ),
                     '<'=,
                     less      = pdist(ts, TRUE),
                     '>'=,
                     greater   = pdist(ts, FALSE)
                     )

    op <- par(mfrow=c(2,1))
    hist(p.vals, main='', xlab='P-Values')

    if( !is.na(alpha) ) {
        abline(v=alpha, col='red')
        title(sub=paste( round(mean(p.vals <= alpha)*100, 1), '% <= ',
              alpha))
    }

    qqplot( seq(along=p.vals)/(B+1), p.vals,
           xlab='Theoretical quantiles of Uniform',
           ylab='P-values')
    abline(0,1, col='grey')

    par(op)

    invisible(p.vals)
}

Pvalue.binom.sim <- function(n=100, p=0.5, p0=0.5,
                            test=c('exact','approx'),
                            alternative=c('two.sided', 'less', 'greater',
                            '<>','!=','<','>'),
                            alpha=0.05, B=1e3) {
    test <- match.arg(test)
    alternative <- match.arg(alternative)

    x <- rbinom(B,n,p)

    pdist <- switch(test,
        exact=function(x, lower.tail) {
            if(lower.tail) {
                pbinom(x, n, p0)
            } else {
                pbinom(pmax(0,x-1), n, p0, lower.tail=FALSE)
            }
        },
        approx=function(x, lower.tail) {
            xbar <- x/n
            ts <- (xbar - p0)/sqrt( p0*(1-p0)/n )
            pnorm(ts, lower.tail=lower.tail)
        }
                    )

    p.vals <- switch(alternative,
                 '!='=,'<>'=,
                 two.sided = pmin(1,2*pmin( pdist(x,TRUE), pdist(x,FALSE) ) ),
                 '<'=,
                 less      = pdist(x, TRUE),
                 '>'=,
                 greater   = pdist(x, FALSE)
                     )
    op <- par(mfrow=c(2,1))
    hist(p.vals, main='', xlab='P-Values') #, col='grey', prob=TRUE)
#    lines( hist(p.vals, breaks=c(0,pbinom(0:n,n,p0)), plot=FALSE),
#          border='green')

    if( !is.na(alpha) ) {
        abline(v=alpha, col='red')
        title(sub=paste( round(mean(p.vals <= alpha)*100, 1), '% <= ',
              alpha))
    }

    qqplot( seq(along=p.vals)/(B+1), p.vals,
           xlab='Theoretical quantiles of Uniform',
           ylab='P-values')
    abline(0,1, col='grey')

    par(op)

    invisible(p.vals)
}






run.Pvalue.norm.sim <- function() {
    lst <- list( Sim=list( n=list('numentry', init=50),
                   mu=list('numentry', init=0),
                   sigma=list('numentry',init=1),
                   B=list('numentry',init=10000),
                   alpha=list('numentry', init=0.05)
                 ),
                 Test=list( test=list('radiobuttons', values=c('z','t'),
                              init='z'),
                            mu0=list('numentry', init=0),
                            sigma0=list('numentry', init=1),
                            alternative=list('radiobuttons',
                              values=c('!=','<','>'), init='!='))
                )

    tkexamp(Pvalue.norm.sim(), lst, plotloc='left')
}


run.Pvalue.binom.sim <- function() {
    lst <- list( Sim=list( n=list('numentry', init=100),
                   p=list('numentry', init=0.5),
                   B=list('numentry',init=1000),
                   alpha=list('numentry', init=0.05)
                 ),
                 Test=list( test=list('radiobuttons',
                            values=c('exact','approx'),
                              init='exact'),
                            p0=list('numentry', init=0.5),
                            alternative=list('radiobuttons',
                              values=c('!=','<','>'), init='!='))
                )

    tkexamp(Pvalue.binom.sim(), lst, plotloc='left')
}

