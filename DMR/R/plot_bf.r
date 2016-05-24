plot_bf <- function(m){
    mod <- m$Best$Model
    n <- nrow(mod$model)
    llik <- m$LogLik
    len <- length(llik)
    d <- ifelse(class(mod) == 'glm', 0, 1)
    bic <- sapply(1:len, function(x) -2*llik[x] + (len - x + d + 1)*log(n))
    bic1 <- exp(-0.5*(bic - min(bic)))
    plot(bic1,length(bic):1, pch = '*', col = 'firebrick', cex = 2, main = 'Approximate Bayes Factors', xlab = 'BF', ylab = '', axes = F)
    axis(2, labels = 1:length(bic), las = 2, at = length(bic):1)
    axis(1, labels = seq(0,1, 0.2), las = 1, at = seq(0,1, 0.2))
    box()
    lines(bic1,length(bic):1, lty = 2, col = 1 )
    abline(v = c(1/3, 1/10), col = 'chartreuse3')
    grid()
    #abline(h = seq(1.5, 12.5, 1), col = 'gray77', lty = 3, xpd = T)
}
