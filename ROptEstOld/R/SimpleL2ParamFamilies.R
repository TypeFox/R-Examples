##################################################################
## Binomial family
##################################################################
BinomFamily <- function(size = 1, prob = 0.5, trafo){ 
    name <- "Binomial family"
    distribution <- Binom(size = size, prob = prob)
    if(prob == 0.5)
        distrSymm <- SphericalSymmetry(SymmCenter = size*prob)
    else
        distrSymm <- NoSymmetry()
    param <- ParamFamParameter(name = "probability of success",  
                            main = prob, trafo = trafo)
    props <- c("The Binomial family is symmetric with respect to prob = 0.5;", 
               "i.e., d(Binom(size, prob))(k)=d(Binom(size,1-prob))(size-k)")
    fct <- function(x){ (x-size*prob)/(prob*(1-prob)) }
    body(fct) <- substitute({ (x-size*prob)/(prob*(1-prob)) },
                        list(size = size, prob = prob))
    L2deriv <- EuclRandVarList(RealRandVariable(Map = list(fct), Domain = Reals()))
    L2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = size*prob)) 
    L2derivDistr <- UnivarDistrList((distribution - size*prob)/(prob*(1-prob)))
    if(prob == 0.5)
        L2derivDistrSymm <- DistrSymmList(SphericalSymmetry(SymmCenter = 0))
    else
        L2derivDistrSymm <- DistrSymmList(NoSymmetry())
    FisherInfo <- PosDefSymmMatrix(matrix(size/(prob*(1-prob))))

    L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, props = props, 
        L2deriv = L2deriv, L2derivSymm = L2derivSymm, 
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo = FisherInfo)
}

##################################################################
## Poisson family
##################################################################
PoisFamily <- function(lambda = 1, trafo){ 
    name <- "Poisson family"
    distribution <- Pois(lambda = lambda)
    distrSymm <- NoSymmetry()
    param <- ParamFamParameter(name = "positive mean",  
                            main = lambda, trafo = trafo)
    props <- character(0)
    fct <- function(x){ x/lambda-1 }
    body(fct) <- substitute({ x/lambda-1 }, list(lambda = lambda))
    L2deriv <- EuclRandVarList(RealRandVariable(Map = list(fct), Domain = Reals()))
    L2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = lambda))
    L2derivDistr <- UnivarDistrList(distribution/lambda - 1)
    L2derivDistrSymm <- DistrSymmList(NoSymmetry())
    FisherInfo <- PosDefSymmMatrix(matrix(1/lambda))

    L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, props = props, 
        L2deriv = L2deriv, L2derivSymm = L2derivSymm, 
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo = FisherInfo)
}

##################################################################
## Normal location family
##################################################################
NormLocationFamily <- function(mean = 0, sd = 1, trafo){ 
    name <- "normal location family"
    distribution <- Norm(mean = mean, sd = sd)
    distrSymm <- SphericalSymmetry(SymmCenter = mean)
    param <- ParamFamParameter(name = "location", main = mean, trafo = trafo)
    props <- c("The normal location family is invariant under",
               "the group of transformations 'g(x) = x + mean'",
               "with location parameter 'mean'")
    fct <- function(x){ (x - mean)/sd^2 }
    body(fct) <- substitute({ (x - mean)/sd^2 }, list(mean = mean, sd = sd))                
    L2deriv <- EuclRandVarList(RealRandVariable(Map = list(fct), Domain = Reals()))
    L2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = mean))
    L2derivDistr <- UnivarDistrList(Norm(mean=0, sd=1/sd))
    L2derivDistrSymm <- DistrSymmList(SphericalSymmetry(SymmCenter = 0))
    FisherInfo <- PosDefSymmMatrix(matrix(1/sd^2))

    L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, props = props, 
        L2deriv = L2deriv, L2derivSymm = L2derivSymm, 
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo = FisherInfo)
}

##################################################################
## Gumbel location family
##################################################################
GumbelLocationFamily <- function(loc = 0, scale = 1, trafo){ 
    name <- "Gumbel location family"
    distribution <- Gumbel(loc = loc, scale = scale)
    distrSymm <- NoSymmetry()
    param <- ParamFamParameter(name = "location", main = loc, trafo = trafo)
    props <- c("The Gumbel location family is invariant under",
               "the group of transformations 'g(x) = x + loc'",
               "with location parameter 'loc'")
    fct <- function(x){ (1 - exp(-(x-loc)/scale))/scale }
    body(fct) <- substitute({ (1 - exp(-(x-loc)/scale))/scale }, 
                         list(loc = loc, scale = scale))                
    L2deriv <- EuclRandVarList(RealRandVariable(Map = list(fct), Domain = Reals()))
    L2derivSymm <- FunSymmList(NonSymmetric())
    L2derivDistr <- UnivarDistrList((1-Exp(rate=1))/scale)
    L2derivDistrSymm <- DistrSymmList(NoSymmetry())
    FisherInfo <- PosDefSymmMatrix(matrix(1/scale^2))

    L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, props = props, 
        L2deriv = L2deriv, L2derivSymm = L2derivSymm, 
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo = FisherInfo)
}

##################################################################
## Normal scale family
##################################################################
NormScaleFamily <- function(sd = 1, mean = 0, trafo){ 
    name <- "normal scale family"
    distribution <- Norm(mean = mean, sd = sd)
    distrSymm <- SphericalSymmetry(SymmCenter = mean)
    param <- ParamFamParameter(name = "scale", main = sd, trafo = trafo)
    props <- c("The normal scale family is invariant under",
               "the group of transformations 'g(y) = sd*y'",
               "with scale parameter 'sd'")
    fct <- function(x){ (((x-mean)/sd)^2 - 1)/sd }
    body(fct) <- substitute({ (((x-mean)/sd)^2 - 1)/sd }, list(sd = sd, mean = mean))                
    L2deriv <- EuclRandVarList(RealRandVariable(Map = list(fct), Domain = Reals()))
    L2derivSymm <- FunSymmList(EvenSymmetric(SymmCenter = mean))
    L2derivDistr <- UnivarDistrList((Chisq(df = 1, ncp = 0)-1)/sd)
    L2derivDistrSymm <- DistrSymmList(NoSymmetry())
    FisherInfo <- PosDefSymmMatrix(matrix(2/sd^2))

    L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, props = props, 
        L2deriv = L2deriv, L2derivSymm = L2derivSymm, 
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo = FisherInfo)
}

##################################################################
## Exponential scale family
##################################################################
ExpScaleFamily <- function(rate = 1, trafo){ 
    name <- "Exponential scale family"
    distribution <- Exp(rate = rate)
    distrSymm <- NoSymmetry()
    param <- ParamFamParameter(name = "scale", main = 1/rate, trafo = trafo)
    props <- c("The Exponential scale family is invariant under",
               "the group of transformations 'g(y) = scale*y'",
               "with scale parameter 'scale = 1/rate'")
    fct <- function(x){ (rate*x - 1)*rate }
    body(fct) <- substitute({ (rate*x - 1)*rate }, list(rate = rate))                
    L2deriv <- EuclRandVarList(RealRandVariable(Map = list(fct), Domain = Reals()))
    L2derivSymm <- FunSymmList(EvenSymmetric(SymmCenter = 1/rate))
    L2derivDistr <- UnivarDistrList((Exp(rate=1)-1)*rate)
    L2derivDistrSymm <- DistrSymmList(NoSymmetry())
    FisherInfo <- PosDefSymmMatrix(matrix(rate^2))

    L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, props = props, 
        L2deriv = L2deriv, L2derivSymm = L2derivSymm, 
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo = FisherInfo)
}

##################################################################
## Lognormal scale family
##################################################################
LnormScaleFamily <- function(meanlog = 0, sdlog = 1, trafo){ 
    name <- "lognormal scale family"
    distribution <- Lnorm(meanlog = meanlog, sdlog = sdlog)
    distrSymm <- NoSymmetry()
    param <- ParamFamParameter(name = "scale", main = exp(meanlog), trafo = trafo)
    props <- c("The Lognormal scale family is invariant under",
               "the group of transformations 'g(y) = scale*y'",
               "with scale parameter 'scale = exp(meanlog)'")
    fct <- function(x){ exp(-meanlog)*(log(x)-meanlog)/sdlog^2 }
    body(fct) <- substitute({ exp(-meanlog)*(log(x) - meanlog)/sdlog^2 },
                       list(meanlog = meanlog, sdlog = sdlog))                
    L2deriv <- EuclRandVarList(RealRandVariable(Map = list(fct), Domain = Reals()))
    L2derivSymm <- FunSymmList(NonSymmetric())
    L2derivDistr <- UnivarDistrList(Norm(mean=0, sd=exp(-meanlog)/sdlog^2))
    L2derivDistrSymm <- DistrSymmList(SphericalSymmetry(SymmCenter = 0))
    FisherInfo <- PosDefSymmMatrix(matrix(exp(-meanlog)^2/sdlog^2))

    L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, props = props, 
        L2deriv = L2deriv, L2derivSymm = L2derivSymm, 
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo = FisherInfo)
}

##################################################################
## Gamma family
##################################################################
GammaFamily <- function(scale = 1, shape = 1, trafo){ 
    name <- "Gamma family"
    distribution <- Gammad(scale = scale, shape = shape)
    distrSymm <- NoSymmetry()
    param <- ParamFamParameter(name = "scale and shape",  
                        main = c(scale, shape), trafo = trafo)
    props <- c("The Gamma family is scale invariant via the parametrization",
               "'(nu,shape)=(log(scale),shape)'")
    fct1 <- function(x){ (x/scale - shape)/scale }
    body(fct1) <- substitute({ (x/scale - shape)/scale },
                        list(scale = scale, shape = shape))
    fct2 <- function(x){ (log(x/scale) - digamma(shape)) }
    body(fct2) <- substitute({ (log(x/scale) - digamma(shape)) },
                        list(scale = scale, shape = shape))
    L2deriv <- EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals()))
    L2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = scale*shape), NonSymmetric())
    L2derivDistr <- UnivarDistrList((Gammad(scale = 1, shape = shape) - shape)/scale, 
                                         (log(Gammad(scale = 1, shape = shape)) - digamma(shape)))
    L2derivDistrSymm <- DistrSymmList(NoSymmetry(), NoSymmetry())
    FisherInfo <- PosDefSymmMatrix(matrix(c(shape/scale^2, 1/scale, 
                                            1/scale, trigamma(shape)), ncol=2))

    L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, props = props, 
        L2deriv = L2deriv, L2derivSymm = L2derivSymm, 
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo = FisherInfo)
}

##################################################################
## Normal location and scale family
##################################################################
NormLocationScaleFamily <- function(mean = 0, sd = 1, trafo){ 
    name <- "normal location and scale family"
    distribution <- Norm(mean = mean, sd = sd)
    distrSymm <- SphericalSymmetry(SymmCenter = mean)
    param <- ParamFamParameter(name = "location and scale", main = c(mean, sd), trafo = trafo)
    props <- c("The normal location and scale family is invariant under",
               "the group of transformations 'g(x) = sd*x + mean'",
               "with location parameter 'mean' and scale parameter 'sd'")
    fct1 <- function(x){ (x - mean)/sd^2 }
    body(fct1) <- substitute({ (x - mean)/sd^2 }, list(mean = mean, sd = sd))                
    fct2 <- function(x){ (((x-mean)/sd)^2 - 1)/sd }
    body(fct2) <- substitute({ (((x-mean)/sd)^2 - 1)/sd }, list(sd = sd, mean = mean))                
    L2deriv <- EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals()))
    L2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = mean), EvenSymmetric(SymmCenter = mean))
    L2derivDistr <- UnivarDistrList(Norm(mean=0, sd=1/sd), (Chisq(df = 1, ncp = 0)-1)/sd)
    L2derivDistrSymm <- DistrSymmList(SphericalSymmetry(), NoSymmetry())
    FisherInfo <- PosDefSymmMatrix(matrix(c(1/sd^2, 0, 0, 2/sd^2), ncol=2))

    L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, props = props, 
        L2deriv = L2deriv, L2derivSymm = L2derivSymm, 
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo = FisherInfo)
}
