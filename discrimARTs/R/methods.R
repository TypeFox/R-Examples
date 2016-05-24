print.discrimARTs <- function(x, ...) {
    .dat <- x[c('method', 'convergence','neglogLik','pars.init','MLE.est')]
    .dat$pars.init <- unlist(.dat$pars.init)
    print(.dat)
}

plot.discrimARTs <- function(x, npoints=1e3, main=NULL, xlab='Measured trait', 
    legend=T, legend.x=0.9, legend.y=0.9, legend.digits=3, legend.fontsize=8, ...) {
    ## Input data for plotting
    .dat <- x$input
    .quants <- seq(from=min(.dat), to=max(.dat), length.out=npoints)
    .est <- as.list(x$MLE.est)
    if (x$method == 'normal' ) {
        dist1 <- with(.est, { 
            (1-mix.prob) * dnorm( .quants, mean=dist1.par1, sd=dist1.par2)
        })
        dist2 <- with(.est, { 
            mix.prob * dnorm( .quants, mean=dist2.par1, sd=dist2.par2)
        })
    } else if (x$method == 'facing.gamma') {
        dist1 <- with(.est, { 
            (1-mix.prob) * dgamma( .quants - x$lower, shape=dist1.par1, scale=dist1.par2)
        })
        dist2 <- with(.est, { 
            mix.prob * dgamma( x$upper - .quants, shape=dist2.par1, scale=dist2.par2)
        })
    } else { stop('Plotting not implemented for method %s', x$method)}

    dist.mix <- dist1 + dist2
    ## Histogram of original observations,
    ## Make sure ylim doesn't cut off distribs
    hist(.dat, freq=FALSE, ylim=c(0, max(dist.mix)), main=main, xlab=xlab)
    ## Original observations
    points(.dat, rep(0, length(.dat)))
    ## over-plot individual dists and mix
    lines( .quants, dist1, lty=3)
    lines( .quants, dist2, lty=3)
    lines( .quants, dist.mix )
    if (legend) {
        ## collapse estimated parameter names and values into one string 
        ## with each par separated by newlines
        .param.text <- paste(paste( names(.est), round(as.numeric(.est), digits=legend.digits), sep='='), collapse='\n')
        if( x$method == 'facing.gamma') {
            ## Add upper and lower
            .bounds.text <- paste(paste( c('Lower bound', 'Upper bound'), c(x$lower, x$upper), sep='='), collapse='\n')
            .param.text <- paste(.param.text, .bounds.text, sep='\n')
        }
        ## Make a "legend" including negative log likelihood and parameter estimates
        grid.text(x=legend.x, y=legend.y, sprintf('Negative logLik = %s\nMLE Parameter Estimates:\n%s', round(x$neglogLik, digits=legend.digits), .param.text), just=c('right', 'top'), gp=gpar(fontsize=legend.fontsize))
    }
    return()
    }
