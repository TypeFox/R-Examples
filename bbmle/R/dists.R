
snorm <- function(mean,sd) {
  list(title="Normal",
       mean=mean,sd=sd,
       median=mean,
       mode=mean,
       variance=sd^2,
       sd=sd)
}

sbinom <- function(size,prob) {
  list(title="Binomial",
       prob=prob,size=size,
       mean=prob*size,
       median=qbinom(0.5,size,prob),
       mode=NA,
       variance=size*prob*(1-prob),
       sd=sqrt(size*prob*(1-prob)),
       formula="x*log(prob)+(size-x)*log(1-prob)")
}

sbeta <- function(shape1,shape2) {
  list(title="Beta",
       shape1=shape1,shape2=shape2,
       mean=shape1/(shape1+shape2),
       median=qbeta(0.5,shape1,shape2),
       mode=NA,
       variance=shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1)),
       sd=sqrt(shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1))))
}

snbinom <- function(size,prob,mu) {
    if (missing(mu) && !missing(prob)) {
        mupar <- FALSE
        mu = NA ## FIXME
        warning("STUB in snbinom: calc. mu as a function of prob")
    }
    if (!missing(mu) && missing(prob)) {
        mupar <- TRUE
        prob = size/(size+mu)
    }
    v <- if (mupar) mu+mu^2/size else size*(1-prob)/prob^2
    list(title="Negative binomial",
         prob=prob,mu=mu,size=size,
         mean=if (mupar) mu else size*(1-prob)/prob,
         median= if (mupar) qnbinom(0.5,mu=mu,size) else qnbinom(0.5,prob=prob,size),
         mode=NA,
         variance=v,
         sd=sqrt(v))
}

spois <- function(lambda) {
  list(title="Poisson",
       lambda=lambda,
       mean=lambda,
       median=qpois(0.5,lambda),
       mode=NA,
       variance=lambda,
       sd=sqrt(lambda))      
}

sbetabinom <- function(size,prob,theta) {
  list(title="Beta-binomial",
       prob=prob,size=size,theta=theta,
       mean=prob*size,
       median=NA, ## qbetabinom(0.5,size,prob),
       mode=NA,
       variance=size*prob*(1-prob)/theta,
       sd=sqrt(size*prob*(1-prob)))
}

sgamma <- function(shape,rate=1,scale=1/rate) {
    if (missing(rate)) rate <- 1/scale
    list(title="Gamma",
         mean=shape/rate,sd=sqrt(shape)/rate,
         median=NA,
         mode=NA,
         variance=shape/rate^2)
}
