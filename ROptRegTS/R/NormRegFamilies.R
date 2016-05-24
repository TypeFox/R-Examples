## generating function
NormLinRegFamily <- function(theta, trafo, RegDistr = Norm(), 
                             RegSymm, Reg2Mom){
    if(missing(theta))
        theta <- numeric(dimension(img(RegDistr)))

    lth <- length(theta)
    if(lth != dimension(img(RegDistr)))
        stop("'theta' has wrong dimension")

    name <- "L2 differentiable linear regression family"
    distribution <- LMCondDistribution(Error = Norm(),
                                theta = theta, scale = 1, intercept = 0) 
    symmfun <- function(cond){ cond%*%theta }
    body(symmfun) <- substitute({ cond%*%theta }, list(theta = theta))
    distrSymm <- SphericalSymmetry(SymmCenter = symmfun)

    param <- ParamFamParameter(name = paste("Parameter of", name), 
                        main = theta, nuisance = NULL, trafo = trafo)
    props <- c("The linear regression family is invariant under",
               "the group of transformations 'g(x,y) = (x, t(x)theta + y)'",
               "with regression parameter 'theta'")

    ErrorDistr <- Norm()
    ErrorSymm <- SphericalSymmetry(SymmCenter = 0)

    fct <- function(x){ as.vector(x[1:k]*(x[k+1] - x[1:k] %*% theta)) }
    body(fct) <- substitute({ as.vector(x[1:k]*(x[k+1] - x[1:k] %*% theta)) },
                        list(k = lth, theta = theta))
    L2deriv <- EuclRandVarList(EuclRandVariable(Map = list(fct), 
                        Domain = EuclideanSpace(dimension = lth + 1),
                        dimension = trunc(lth)))

    Regressor <- EuclRandVariable(Map = list(function(x){x}), 
                        Domain = EuclideanSpace(dimension = trunc(lth)),
                        dimension = trunc(lth))
    if(missing(RegSymm)) RegSymm <- NoSymmetry()

    ErrorL2deriv <- EuclRandVarList(RealRandVariable(Map = list(function(x){x}), Domain = Reals()))
    ErrorL2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = 0))
    ErrorL2derivDistr <- UnivarDistrList(Norm())
    ErrorL2derivDistrSymm <- DistrSymmList(SphericalSymmetry(SymmCenter = 0))

    if(missing(Reg2Mom))
        FisherInfo <- PosDefSymmMatrix(E(RegDistr, function(x){x %*% t(x)}))
    else
        FisherInfo <- PosDefSymmMatrix(Reg2Mom)

    new("L2RegTypeFamily", name = name, distribution = distribution, distrSymm = distrSymm, 
        param = param, props = props, L2deriv = L2deriv, ErrorDistr = ErrorDistr, 
        RegDistr = RegDistr, ErrorSymm = ErrorSymm, RegSymm = RegSymm, Regressor = Regressor, 
        ErrorL2deriv = ErrorL2deriv, ErrorL2derivSymm = ErrorL2derivSymm, 
        ErrorL2derivDistr = ErrorL2derivDistr, ErrorL2derivDistrSymm = ErrorL2derivDistrSymm, 
        FisherInfo = FisherInfo)
}

## generating function
NormLinRegScaleFamily <- function(theta, scale = 1, trafo, RegDistr = Norm(), 
                                  RegSymm, Reg2Mom, nuisance = FALSE){
    if(missing(theta))
        theta <- numeric(dimension(img(RegDistr)))

    lth <- length(theta)
    if(lth != dimension(img(RegDistr)))
        stop("'theta' has wrong dimension")

    name <- "L2 differentiable linear regression family with unknown scale"
    distribution <- LMCondDistribution(Error = Norm(),
                                theta = theta, scale = scale, intercept = 0) 
    symmfun <- function(cond){ cond%*%theta }
    body(symmfun) <- substitute({ cond%*%theta }, list(theta = theta))
    distrSymm <- SphericalSymmetry(SymmCenter = symmfun)

    if(nuisance)
        param <- ParamFamParameter(name = paste("Parameter of", name), 
                            main = theta, nuisance = scale, trafo = trafo)
    else
        param <- ParamFamParameter(name = paste("Parameter of", name), 
                            main = c(theta, scale), nuisance = NULL, trafo = trafo)

    props <- c("The linear regression family with unknown scale is invariant", 
               "under the group of transformations 'g(x,y) = (x, t(x)theta + scale y)'",
               "with regression parameter 'theta' and scale parameter 'scale'")

    ErrorDistr <- Norm()
    ErrorSymm <- SphericalSymmetry(SymmCenter = 0)

    fct1 <- function(x){ as.vector(x[1:k]*(x[k+1] - x[1:k] %*% theta)/sd^2) }
    body(fct1) <- substitute({ as.vector(x[1:k]*(x[k+1] - x[1:k] %*% theta)/sd^2) },
                        list(k = lth, theta = theta, sd = scale))
    fct2 <- function(x){ as.vector((((x[k+1] - x[1:k] %*% theta)/sd)^2 - 1)/sd) }
    body(fct2) <- substitute({ as.vector((((x[k+1] - x[1:k] %*% theta)/sd)^2 - 1)/sd) }, 
                        list(k = lth, theta = theta, sd = scale))

    L2deriv <- EuclRandVarList(EuclRandVariable(Map = list(fct1), 
                                         Domain = EuclideanSpace(dimension = lth+1),
                                         dimension = trunc(lth)),
                               EuclRandVariable(Map = list(fct2),
                                         Domain = EuclideanSpace(dimension = lth+1),
                                         dimension = 1))

    Regressor <- EuclRandVariable(Map = list(function(x){x}), 
                        Domain = EuclideanSpace(dimension = trunc(lth)),
                        dimension = trunc(lth))
    if(missing(RegSymm)) RegSymm <- NoSymmetry()

    fct1 <- function(x){ x/sd }
    body(fct1) <- substitute({ x/sd }, list(sd = scale))
    fct2 <- function(x){ (x^2 - 1)/sd }
    body(fct2) <- substitute({ (x^2 - 1)/sd }, list(sd = scale))
    ErrorL2deriv <- EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals()))
    ErrorL2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = 0), EvenSymmetric(SymmCenter = 0))
    ErrorL2derivDistr <- UnivarDistrList(Norm(sd = scale), (Chisq()-1)/scale)
    ErrorL2derivDistrSymm <- DistrSymmList(SphericalSymmetry(SymmCenter = 0), NoSymmetry())

    FisherInfo <- matrix(0, ncol = lth + 1, nrow = lth + 1)
    if(missing(Reg2Mom)){
        K <- PosDefSymmMatrix(E(RegDistr, function(x){x %*% t(x)}))
        FisherInfo[1:lth, 1:lth] <- K
        FisherInfo[(lth+1), (lth+1)] <- 2
        FisherInfo <- PosDefSymmMatrix(FisherInfo/scale^2)
    }else{
        if(ncol(Reg2Mom) != lth || nrow(Reg2Mom) != lth)
            stop("'Reg2Mom' has wrong dimension")
        FisherInfo[1:lth, 1:lth] <- Reg2Mom
        FisherInfo[(lth+1), (lth+1)] <- 2
        FisherInfo <- PosDefSymmMatrix(FisherInfo/scale^2)
    }

    new("L2RegTypeFamily", name = name, distribution = distribution, distrSymm = distrSymm, 
        param = param, props = props, L2deriv = L2deriv, ErrorDistr = ErrorDistr, 
        RegDistr = RegDistr, ErrorSymm = ErrorSymm, RegSymm = RegSymm, Regressor = Regressor, 
        ErrorL2deriv = ErrorL2deriv, ErrorL2derivSymm = ErrorL2derivSymm, 
        ErrorL2derivDistr = ErrorL2derivDistr, ErrorL2derivDistrSymm = ErrorL2derivDistrSymm, 
        FisherInfo = FisherInfo)
}

## generating function
NormLinRegInterceptFamily <- function(theta, intercept = 0, trafo, RegDistr = Norm(), 
                                      RegSymm, Reg2Mom, nuisance = FALSE){
    mu <- E(RegDistr)
    if(!identical(all.equal(mu, numeric(length(mu))), TRUE))
        stop("'RegDistr' has expectation != 0")

    if(missing(theta))
        theta <- numeric(dimension(img(RegDistr)))

    lth <- length(theta)
    if(lth != dimension(img(RegDistr)))
        stop("'theta' has wrong dimension")

    name <- "L2 differentiable linear regression family with unknown intercept"
    distribution <- LMCondDistribution(Error = Norm(),
                                theta = theta, scale = 1, intercept = intercept) 
    symmfun <- function(cond){ cond%*%theta + mu}
    body(symmfun) <- substitute({ cond%*%theta + mu}, list(theta = theta, mu = intercept))
    distrSymm <- SphericalSymmetry(SymmCenter = symmfun)

    if(nuisance)
        param <- ParamFamParameter(name = paste("Parameter of", name), 
                            main = theta, nuisance = intercept, trafo = trafo)
    else
        param <- ParamFamParameter(name = paste("Parameter of", name), 
                            main = c(theta, intercept), nuisance = NULL, trafo = trafo)

    props <- c("The linear regression family with unknown intercept is invariant", 
               "under the group of transformations 'g(x,y) = (x, intercept + t(x)theta + y)'",
               "with regression parameter 'theta' and intercept parameter 'intercept'")

    ErrorDistr <- Norm()
    ErrorSymm <- SphericalSymmetry(SymmCenter = 0)

    fct1 <- function(x){ as.vector(x[1:k]*(x[k+1] - x[1:k] %*% theta - mu)) }
    body(fct1) <- substitute({ as.vector(x[1:k]*(x[k+1] - x[1:k] %*% theta - mu)) },
                        list(k = lth, theta = theta, mu = intercept))
    fct2 <- function(x){ as.vector(x[k+1] - x[1:k] %*% theta - mu) }
    body(fct2) <- substitute({ as.vector(x[k+1] - x[1:k] %*% theta - mu) }, 
                        list(k = lth, theta = theta, mu = intercept))

    L2deriv <- EuclRandVarList(EuclRandVariable(Map = list(fct1), 
                                         Domain = EuclideanSpace(dimension = lth+1),
                                         dimension = trunc(lth)),
                               EuclRandVariable(Map = list(fct2),
                                         Domain = EuclideanSpace(dimension = lth+1),
                                         dimension = 1))

    Regressor <- EuclRandVariable(Map = list(function(x){x}), 
                        Domain = EuclideanSpace(dimension = trunc(lth)),
                        dimension = trunc(lth))
    if(missing(RegSymm)) RegSymm <- NoSymmetry()

    fct1 <- function(x){ x }
    ErrorL2deriv <- EuclRandVarList(RealRandVariable(Map = list(fct1, fct1), Domain = Reals()))
    ErrorL2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = 0), OddSymmetric(SymmCenter = 0))
    ErrorL2derivDistr <- UnivarDistrList(Norm(), Norm())
    ErrorL2derivDistrSymm <- DistrSymmList(SphericalSymmetry(SymmCenter = 0), 
                                           SphericalSymmetry(SymmCenter = 0))

    FisherInfo <- matrix(0, ncol = lth + 1, nrow = lth + 1)
    if(missing(Reg2Mom)){
        K <- PosDefSymmMatrix(E(RegDistr, function(x){x %*% t(x)}))
        FisherInfo[1:lth, 1:lth] <- K
        FisherInfo[(lth+1), (lth+1)] <- 1
        FisherInfo <- PosDefSymmMatrix(FisherInfo)
    }else{
        if(ncol(Reg2Mom) != lth || nrow(Reg2Mom) != lth)
            stop("'Reg2Mom' has wrong dimension")
        FisherInfo[1:lth, 1:lth] <- Reg2Mom
        FisherInfo[(lth+1), (lth+1)] <- 1
        FisherInfo <- PosDefSymmMatrix(FisherInfo)
    }

    new("L2RegTypeFamily", name = name, distribution = distribution, distrSymm = distrSymm, 
        param = param, props = props, L2deriv = L2deriv, ErrorDistr = ErrorDistr, 
        RegDistr = RegDistr, ErrorSymm = ErrorSymm, RegSymm = RegSymm, Regressor = Regressor, 
        ErrorL2deriv = ErrorL2deriv, ErrorL2derivSymm = ErrorL2derivSymm, 
        ErrorL2derivDistr = ErrorL2derivDistr, ErrorL2derivDistrSymm = ErrorL2derivDistrSymm, 
        FisherInfo = FisherInfo)
}
