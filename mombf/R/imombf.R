###
### imombf.R
###

imombf <- function(lm1,
                   coef,
                   g,
                   prior.mode,
                   nu=1,
                   theta0,
                   method='adapt',
                   nquant=100,
                   B=10^5) {
    UseMethod("imombf")
}

