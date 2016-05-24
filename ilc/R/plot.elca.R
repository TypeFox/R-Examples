plot.elca <-
function(x, ...){
    age <- x$age
    year <- x$year
    n <- length(year); k <- length(age); g <- length(x$ag)
    ax <- with(x, array(ax+ag[gl(g, k)], dim=c(k,g)))
    col <- rainbow(g, start=0.05)
    txt <- paste('X', names(x$ag), sep='')
    txt[nchar(txt)==2] <- paste(txt[nchar(txt)==2], '')
    txt <- paste(txt, rb(round(x$ag, 2)))
    xlim <- range(age); xlim[1] <- xlim[1]-2; xlim[2] <- xlim[2]+2
    par.old <- par(no.readonly=T) 
    layout(matrix(c(1,1,2,3), 2, 2, byrow = F))
    l <- 0.8
    par(mar = c(4, 4, 4, 1.5) + 0.1)
    matplot(age, ax, type='l', col=col, ylab=expression(alpha['x,g']), xlab='Age',
      lty=seq(g), xlim=xlim, lwd=1.2)
    title(main = 'Main effects', line=l)
    legend(coord('UL'), legend=txt, lty=seq(g), y.intersp=0.95, x.intersp=0.5, # cex=0.95,
      text.col=col, col=col, lwd=1.2)
    text(age[k], ax[k,1], label=' - X1', col=col[1], adj=c(0, 0.5), cex=0.75)
    plot(age, x$bx, type='l', ylab = expression(beta[x]), xlab = 'Age')
    title(main='Interaction effects', line=l)
    par(mar = c(4, 4, 2.5, 1.5) + 0.1)
    plot(x$kt, ylab=substitute(kappa[t][' '] (a), list(a=x$adj)), xlab='Calendar year')
    title(main='Period effects', line=l)
    title(paste('Adjusted LC Regression for', x$label),
      outer=T, line=-1.3, cex.main=1.2, font.main=4)
    invisible(par(par.old))  
}
