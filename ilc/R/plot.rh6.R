plot.rh6 <-
function(obj){
    age <- obj$age
    year <- obj$year
    par.old <- par(no.readonly=T) 
    layout(matrix(c(1,2,3,3), 2, 2, byrow=T))
    # title line:
    l <- 0.8
    par(mar = c(4, 5, 4, 1.5) + 0.1)
    plot(age, obj$ax, type='l', ylab=expression(alpha[x]), xlab='Age')
    title(main = 'Main effects', line=l)
    plot(age, obj$bx, type='l', ylab=expression(beta[x]^(1)), xlab='Age')
    title(main = 'Interaction effects', line=l)
    par(mar = c(5, 4, 2.5, 1.5) + 0.1)  
    plot(obj$kt, ylab=substitute(kappa[t][' '] (a), list(a=obj$adj)), xlab='Calendar year')
    title(main='Period effects', line=l)
    title(paste('Standard LC Regression for', obj$label, sqb(names(obj)[4])),
      outer=T, line=-1.3, cex.main=1.1, font.main=4)
    invisible(par(par.old))  
}
