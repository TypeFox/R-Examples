
# Parameters and weighted sums

prof.gb2 <- function(x, shape1, scale, w=rep(1, length(x))){
sw <- sum(w)
y <- (x/scale)^shape1
slog    <- sum(w*log(y))/sw                  # = m
sloga   <- sum(w*log(1+y))/sw
sloga.a <- sum(w*log(y)*y/(1+y))/sw
r <- sum(w*y/(1+y))/sw
s <- 1/(sloga.a-r*slog)
p <- r*s
q <- (1-r)*s
return(c(r, s, p, q, slog, sloga))
}

# Profile log-likelihood of a, b

proflogl.gb2 <- function(x, shape1, scale, w=rep(1, length(x))){
pars <- prof.gb2(x, shape1, scale, w)
pll <- -lbeta(pars[3],pars[4]) + log(abs(shape1)/scale) + (pars[3]-1/shape1)*pars[5] - (pars[3]+pars[4])*pars[6]
return(pll)
}

# Scores for the profile log-likelihood

profscores.gb2 <- function(x, shape1, scale, w=rep(1, length(x))){
sw <- sum(w)
y <- (x/scale)^shape1
pars <- prof.gb2(x, shape1, scale, w)
dr.da <- (1/shape1)*sum(w*y*log(y)/(1+y)^2)/sw
ds.da <-  -(pars[2]^2/shape1)*sum(w*(y*log(y)/(1+y)^2 + y/(1+y))*(log(y)-pars[5]))/sw
dr.db <- -(shape1/scale)*sum(w*y/(1+y)^2)/sw
ds.db <- (pars[2]^2)*(shape1/scale)*sum(w*y*(log(y)-pars[5])/(1+y)^2)/sw
dpll.dr <- pars[2]*(pars[5] - digamma(pars[3]) + digamma(pars[4]))
dpll.ds <- pars[1]*(pars[5] - digamma(pars[3])) + digamma(pars[2]) - (1-pars[1])*digamma(pars[4]) - pars[6]
dpll.da <- dpll.dr*dr.da + dpll.ds*ds.da
dpll.db <- dpll.dr*dr.db + dpll.ds*ds.db
return(c(dpll.da,dpll.db))
}