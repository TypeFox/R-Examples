plot.rh <-
function(x, ...){
    age <- x$age
    year <- x$year
    par.old <- par(no.readonly=T) 
    layout(matrix(c(rep(1,4),rep(2,4),rep(3,4),rep(4,5), rep(5,7)), 2, 12,
      byrow = T))
    l <- 0.8
    par(mar = c(4, 5, 4, 1.5) + 0.1)
    plot(age, x$ax, type='l', ylab=expression(alpha[x]), xlab='Age')
    title(main = 'Main age effects', line=l)
    plot(age, x$bx, type='l', ylab=expression(beta[x]^(1)), xlab='Age')
    title(main = 'Period Interaction effects', line=l)
    plot(age, x$bx0, type='l', ylab=expression(beta[x]^(0)), xlab='Age')
    title(main = 'Cohort Interaction effects', line=l)
    par(mar = c(5, 4, 2.5, 1.5) + 0.1)
    plot(x$kt, ylab=substitute(kappa[t][' '] (a), list(a=x$adj)), xlab='Calendar year')
    title(main='Period effects', line=l)
    plot(x$itx, ylab=substitute(iota[t-x][' '] (a), list(a=x$adj)), xlab='Year of birth')
    title(main='Cohort effects', line=l)
    tit <- paste(c('Age', if (any(bool(x$kt))) 'Period',
            if (any(bool(x$itx))) 'Cohort'), collapse='-')
    title(paste(tit, 'LC Regression for', x$label, sqb(names(x)[4])),
      outer=T, line=-1.3, cex.main=1.4, font.main=4)
    invisible(par(par.old))  
}
