##################################################################
## Binomial family
##################################################################
BinomFamily <- function(size = 1, prob = 0.5, trafo){ 
    name <- "Binomial family"
    distribution <- Binom(size = size, prob = prob)
    if(.isEqual(prob,0.5))
        distrSymm <- SphericalSymmetry(SymmCenter = size*prob)
    else
        distrSymm <- NoSymmetry()
    param0 <- prob
    names(param0) <- "prob"
    param1 <- size
    names(param1) <- "size"
    if(missing(trafo)) trafo <- matrix(1, dimnames = list("prob","prob"))
    param <- ParamFamParameter(name = "probability of success",  
                               main = param0, 
                               fixed = param1, 
                               trafo = trafo)
    modifyParam <- function(theta){ Binom(size = size, prob = theta) }
    body(modifyParam) <- substitute({ Binom(size = size, prob = theta) }, list(size = size))
    props <- c("The Binomial family is symmetric with respect to prob = 0.5;", 
               "i.e., d(Binom(size, prob))(k)=d(Binom(size,1-prob))(size-k)")
    
    startPar <- function(x,...) c(.Machine$double.eps,1-.Machine$double.eps)
    makeOKPar <- function(param) {if(param<=0) return(.Machine$double.eps)
                                  if(param>=1) return(1-.Machine$double.eps)
                                  return(param)}
    L2deriv.fct <- function(param){
                   prob <- main(param)
                   fct <- function(x){}
                   body(fct) <- substitute({ (x-size*prob)/(prob*(1-prob)) },
                                list(size = size, prob = prob))
                   return(fct)}
    L2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = size*prob)) 
    L2derivDistr <- UnivarDistrList((distribution - size*prob)/(prob*(1-prob)))
    if(.isEqual(prob,0.5))
        L2derivDistrSymm <- DistrSymmList(SphericalSymmetry(SymmCenter = 0))
    else
        L2derivDistrSymm <- DistrSymmList(NoSymmetry())
    FisherInfo.fct <- function(param){
                       prob <- main(param)
                       PosDefSymmMatrix(matrix(size/(prob*(1-prob)),
                           dimnames=list("prob","prob")))}

    FisherInfo <- FisherInfo.fct(param)
    res <- L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, modifyParam = modifyParam,
        props = props, L2deriv.fct = L2deriv.fct, L2derivSymm = L2derivSymm,
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo.fct = FisherInfo.fct, FisherInfo = FisherInfo,
        startPar = startPar, makeOKPar = makeOKPar, 
        .returnClsName = "BinomFamily")
    if(!is.function(trafo))
       f.call <- substitute(BinomFamily(size = s, prob = p,
  	                     trafo = matrix(Tr, dimnames = list("prob","prob"))),
  	                     list(s = size, p = prob, Tr = trafo))    
    else
       f.call <- substitute(BinomFamily(size = s, prob = p,
  	                     trafo = Tr), list(s = size, p = prob, Tr = trafo))    
    
    res@fam.call <- f.call
    return(res)
}

##################################################################
## Poisson family
##################################################################
PoisFamily <- function(lambda = 1, trafo){ 
    name <- "Poisson family"
    distribution <- Pois(lambda = lambda)
    distrSymm <- NoSymmetry()
    param0 <- lambda
    names(param0) <- "lambda"
    if(missing(trafo)) trafo <- matrix(1, dimnames = list("lambda","lambda"))
    param <- ParamFamParameter(name = "positive mean",
                               main = param0, 
                               trafo = trafo)
    modifyParam <- function(theta){ Pois(lambda = theta) }
    props <- character(0)
    startPar <- function(x,...) c(.Machine$double.eps,max(x))
    makeOKPar <- function(param) {if(param<=0) return(.Machine$double.eps)
                                  return(param)}
    L2deriv.fct <- function(param){
                   lambda <- main(param)
                   fct <- function(x){}
                   body(fct) <- substitute({ x/lambda-1 }, 
                                list(lambda = lambda))
                   return(fct)}
    L2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = lambda))
    L2derivDistr <- UnivarDistrList(distribution/lambda - 1)
    L2derivDistrSymm <- DistrSymmList(NoSymmetry())
    FisherInfo.fct <- function(param){
                   lambda <- main(param)
                   PosDefSymmMatrix(matrix(1/lambda,
                           dimnames=list("lambda","lambda")))}

    FisherInfo <- FisherInfo.fct(param)
    res <- L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, modifyParam = modifyParam,
        props = props, L2deriv.fct = L2deriv.fct, L2derivSymm = L2derivSymm,
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo.fct = FisherInfo.fct, FisherInfo = FisherInfo,
        startPar = startPar, makeOKPar = makeOKPar, 
        .returnClsName = "PoisFamily")
    if(!is.function(trafo))
       f.call <- substitute(PoisFamily(lambda = l,
                         trafo = matrix(Tr, dimnames = list("lambda","lambda"))),
                         list(l = lambda, Tr = trafo))
    else
       f.call <- substitute(PoisFamily(lambda = l, trafo = Tr),
                         list(l = lambda, Tr = trafo))

    res@fam.call <- f.call
    return(res)
}

##################################################################
## NegBinomial family
##################################################################
NbinomFamily <- function(size = 1, prob = 0.5, trafo){ 
    name <- "Negative Binomial family"
    distribution <- Nbinom(size = size, prob = prob)
    distrSymm <- NoSymmetry()
    param0 <- prob
    names(param0) <- "prob"
    param1 <- size
    names(param1) <- "size"
    if(missing(trafo)) trafo <- matrix(1, dimnames = list("prob","prob"))
    param <- ParamFamParameter(name = "probability of success",  
                               main = param0, 
                               fixed = param1, 
                               trafo = trafo)
    modifyParam <- function(theta){ Nbinom(size = size, prob = theta) }
    body(modifyParam) <- substitute({ Nbinom(size = size, prob = theta) }, list(size = size))
    props <- ""
    
    startPar <- function(x,...) c(.Machine$double.eps,1-.Machine$double.eps)
    makeOKPar <- function(param) {if(param<=0) return(.Machine$double.eps)
                                  if(param>=1) return(1-.Machine$double.eps)
                                  return(param)}
    L2deriv.fct <- function(param){
                   prob <- main(param)
                   fct <- function(x){}
                   body(fct) <- substitute({ (size/prob- x/(1-prob)) },
                                list(size = size, prob = prob))
                   return(fct)}
    L2derivSymm <- FunSymmList(NonSymmetric()) 
    L2derivDistr <- UnivarDistrList((size/prob- distribution/(1-prob)))
    L2derivDistrSymm <- DistrSymmList(NoSymmetry())
    FisherInfo.fct <- function(param){
                       prob <- main(param)
                       PosDefSymmMatrix(matrix(size/(prob^2*(1-prob)),
                           dimnames=list("prob","prob")))}

    FisherInfo <- FisherInfo.fct(param)
    res <- L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, modifyParam = modifyParam,
        props = props, L2deriv.fct = L2deriv.fct, L2derivSymm = L2derivSymm,
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo.fct = FisherInfo.fct, FisherInfo = FisherInfo,
        startPar = startPar, makeOKPar = makeOKPar, 
        .returnClsName = "NbinomFamily")
    if(!is.function(trafo))
       f.call <- substitute(NbinomFamily(size = s, prob = p,
  	                     trafo = matrix(Tr, dimnames = list("prob","prob"))),
  	                     list(s = size, p = prob, Tr = trafo))    
    else
       f.call <- substitute(NbnomFamily(size = s, prob = p,
  	                     trafo = Tr), list(s = size, p = prob, Tr = trafo))    
    
    res@fam.call <- f.call
    return(res)
}


NbinomwithSizeFamily <- function(size = 1, prob = 0.5, trafo,
                withL2derivDistr = TRUE){
    name <- "Negative Binomial family"
    distribution <- Nbinom(size = size, prob = prob)
    distrSymm <- NoSymmetry()
    param0 <- c(size,prob)
    names(param0) <- nms <- c("size","prob")
    if(missing(trafo)) trafo <- matrix(c(1,0,0,1),2,2, dimnames = list(c("size","prob"),c("size","prob")))
    param <- ParamFamParameter(name = "NegBinomParameter",  
                               main = param0, 
                               trafo = trafo)
    modifyParam <- function(theta){ Nbinom(size = theta[1], prob = theta[2]) }
    body(modifyParam) <- substitute({ Nbinom(size = theta[1], prob = theta[2]) })
    props <- ""
    
    startPar <- function(x,...){ param1 <- c(1,0.5)
                                 names(param1) <- c("size","prob")
                                 return(param1)}
    makeOKPar <- function(param) {if(param["prob"]<=0) param["prob"] <- .Machine$double.eps
                                  if(param["prob"]>=1) param["prob"] <- (1-.Machine$double.eps)
                                  param["size"] <- min(1e-8, param["size"])
                                  return(param)}
    L2deriv.fct <- function(param){
                   prob <- main(param)["prob"]
                   size <- main(param)["size"]
                   fct1 <- function(x){}
                   fct2 <- function(x){}
                   body(fct2) <- substitute({ (size/prob- x/(1-prob)) },
                                list(size = size, prob = prob))
                   body(fct1) <- substitute({ digamma(x+size)-digamma(size)+log(prob)},
                                list(size = size, prob = prob))
                   return(list(fct1, fct2))}

    L2derivSymm <- FunSymmList(NonSymmetric(), NonSymmetric())
    L2derivDistr <- NULL
    if(withL2derivDistr)
       L2derivDistr <- UnivarDistrList( digamma(distribution+size)-
                                        digamma(size)+log(prob),
                                    (size/prob- distribution/(1-prob)))
    L2derivDistrSymm <- DistrSymmList(NoSymmetry(), NoSymmetry())

    FisherInfo.fct <- function(param){
                   prob <- main(param)["prob"]
                   size <- main(param)["size"]
                   xn <- 0:min(max(support(distribution)),
                               qnbinom(1e-6,size=size,prob=prob,lower.tail=FALSE),
                               1e5)
                   I11 <- -sum((trigamma(xn+size)-trigamma(size))*dnbinom(xn,size=size,prob=prob))
                   I12 <- -1/prob
                   I22 <- size/prob^2/(1-prob)
                   PosDefSymmMatrix(matrix(c(I11,I12,I12,I22),2,2,
                           dimnames=list(nms,nms)))}

    FisherInfo <- FisherInfo.fct(param)
    res <- L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, modifyParam = modifyParam,
        props = props, L2deriv.fct = L2deriv.fct, L2derivSymm = L2derivSymm,
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo.fct = FisherInfo.fct, FisherInfo = FisherInfo,
        startPar = startPar, makeOKPar = makeOKPar, 
        .returnClsName = "NbinomwithSizeFamily")
    if(!is.function(trafo))
       f.call <- substitute(NbinomwithSizeFamily(size = s, prob = p,
  	                     trafo = matrix(Tr, dimnames = DN)),
  	                     list(s = size, p = prob, Tr = trafo,
                         DN = list(nms,nms)))    
    else
       f.call <- substitute(NbnomwithSizeFamily(size = s, prob = p,
  	                     trafo = Tr), list(s = size, p = prob, Tr = trafo))    
    
    res@fam.call <- f.call
    return(res)
}

NbinomMeanSizeFamily <- function(size = 1, mean = .5, trafo,
                                 withL2derivDistr = TRUE){
    name <- "Negative Binomial family"
    prob.0 <- size/(size+mean)
    distribution <- Nbinom(size = size, prob = size/(size+mean))
    distrSymm <- NoSymmetry()
    param0 <- c(size,mean)
    names(param0) <- nms <- c("size","mean")
    if(missing(trafo)) trafo <- matrix(c(1,0,0,1),2,2, dimnames = list(nms,nms))
    param <- ParamFamParameter(name = "probability of success",  
                               main = param0, 
                               trafo = trafo)
    modifyParam <- function(theta){ Nbinom(size = theta[1], prob = theta[1]/(theta[1]+theta[2])) }
    body(modifyParam) <- substitute({ Nbinom(size = theta[1], prob = theta[1]/(theta[1]+theta[2])) })
    props <- ""
    
    startPar <- function(x,...){ param1 <- c(1,0.5)
                                 names(param1) <- c("size","mean")
                                 return(param1)}
    makeOKPar <- function(param) {if(param["mean"]<=0) param["mean"] <- .Machine$double.eps
                                  if(param["mean"]>=1) param["mean"] <- (1-.Machine$double.eps)
                                  param["size"] <- min(1e-8, param["size"])
                                  return(param)}
    L2deriv.fct <- function(param){
                   size.0 <- main(param)["size"]
                   mean.0 <- main(param)["mean"]
                   prob.0 <- size.0/(size.0+mean.0)
                   
                   fct1 <- function(x){}
                   fct1.2 <- function(x){}
                   fct2 <- function(x){}
                   body(fct1) <- substitute({ digamma(x+size)-digamma(size)+log(prob.2)},
                                list(size = size.0, prob.2 = prob.0))
                   body(fct1.2)<- substitute({ (size/prob.2- x/(1-prob.2)) },
                                list(size = size.0, prob.2 = prob.0))
                   body(fct2)  <- substitute({ (1/prob.2-1)* fct1(x) - size/prob.2^2 * fct1.2(x)},
                                list(size = size.0, prob.2 = prob.0))
                   return(list(fct1, fct2))}
    L2derivSymm <- FunSymmList(NonSymmetric(), NonSymmetric())

    .di1 <- function(x)  digamma(x+size)-digamma(size)+log(prob.0)
    .di2 <- function(x) .di1(x)*(1/prob.0-1)- (size/prob.0- x/(1-prob.0))*size/prob.0^2 
    .supp1 <- support(distribution)
    .supp0 <- .di2(.supp1)
    .prob1 <- aggregate(data.frame(prob(as(distribution,"DiscreteDistribution"))),
                 by=list(round(.supp0,5)),sum)[,2]
    .Di2 <- DiscreteDistribution( supp=.supp0, prob=.prob1, .withArith = TRUE)
    L2derivDistr <- NULL
    if(withL2derivDistr)
       L2derivDistr <- UnivarDistrList( digamma(distribution+size)-digamma(size)+log(prob.0),
                                     .Di2)
                                     
    L2derivDistrSymm <- DistrSymmList(NoSymmetry(), NoSymmetry())
    FisherInfo.fct <- function(param){
                   mean <- main(param)["mean"]
                   size <- main(param)["size"]
                   prob.0 <- size/(size+mean)
                   xn <- 0:min(max(support(distribution)),
                               qnbinom(1e-6,size=size,prob=prob.0,lower.tail=FALSE),
                               1e5)
                   I11 <- -sum((trigamma(xn+size)-trigamma(size))*dnbinom(xn,size=size,prob=prob.0))
                   I12 <- -1/prob.0
                   I22 <- size/prob.0^2/(1-prob.0)
                   D.m <- matrix(c(1,1/prob.0-1,0,-size/prob.0^2),2,2)
                   ma  <- D.m%*%matrix(c(I11,I12,I12,I22),2,2)%*%t(D.m)
                   dimnames(ma) <- list(nms,nms)
                   PosDefSymmMatrix(ma)}

    FisherInfo <- FisherInfo.fct(param)
    res <- L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, modifyParam = modifyParam,
        props = props, L2deriv.fct = L2deriv.fct, L2derivSymm = L2derivSymm,
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo.fct = FisherInfo.fct, FisherInfo = FisherInfo,
        startPar = startPar, makeOKPar = makeOKPar, 
        .returnClsName = "NbinomMeanSizeFamily")
    if(!is.function(trafo)){
       f.call <- substitute(NbinomMeanSizeFamily(size = s, mean = m,
  	                     trafo = matrix(Tr, dimnames = DN)),
  	                     list(s = size, m = mean, Tr = trafo, DN = list(nms,nms)))    
    }else{
       f.call <- substitute(NbinomMeanSizeFamily(size = s, mean = m,
  	                     trafo = Tr), list(s = size, m= mean, Tr = trafo))    
    
    }
    res@fam.call <- f.call
    return(res)
}

##################################################################
## Gamma family
##################################################################
GammaFamily <- function(scale = 1, shape = 1, trafo, withL2derivDistr = TRUE){
    name <- "Gamma family"
    distribution <- Gammad(scale = scale, shape = shape)
    distrSymm <- NoSymmetry()
    param0 <- c(scale, shape)
    names(param0) <- nms <- c("scale", "shape")
    if(missing(trafo)) {trafo <- diag(2); dimnames(trafo) <-list(nms,nms)}
    param <- ParamFamParameter(name = "scale and shape",  
                        main = param0, trafo = trafo,
                               withPosRestr = TRUE,
                               .returnClsName ="ParamWithScaleAndShapeFamParameter")
    modifyParam <- function(theta){ Gammad(scale = theta[1], shape = theta[2]) }
    props <- c("The Gamma family is scale invariant via the parametrization",
               "'(nu,shape)=(log(scale),shape)'")
    startPar <- function(x,...){ x <- pmax(0,x)
                              E1 <- mean(x)
                              V <- var(x)
                              st <- c(V/E1,E1^2/V)
                              names(st) <- nms
                              return(st)               
                              }
    makeOKPar <- function(param) {param <- abs(param)
                                  return(param)}
    L2deriv.fct <- function(param){
                   scale <- main(param)[1]
                   shape <- main(param)[2]
                   fct1 <- function(x){}
                   fct2 <- function(x){}
                   body(fct1) <- substitute({ (x/scale - shape)/scale },
                        list(scale = scale, shape = shape))
                   body(fct2) <- substitute({ log(x/scale) - digamma(shape) },
                        list(scale = scale, shape = shape))
                   return(list(fct1, fct2))}
    L2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = scale*shape),
                               NonSymmetric())
    L2derivDistr <- NULL
    if(withL2derivDistr)
       L2derivDistr <- UnivarDistrList((Gammad(scale = 1, shape = shape) -
                                            shape)/scale, 
                                    (log(Gammad(scale = 1, shape = shape)) - 
                                            digamma(shape)))
    L2derivDistrSymm <- DistrSymmList(NoSymmetry(), NoSymmetry())
    FisherInfo.fct <- function(param){
                   scale <- main(param)[1]
                   shape <- main(param)[2]
                   PosDefSymmMatrix(matrix(c(shape/scale^2, 1/scale, 
                                            1/scale, trigamma(shape)), ncol=2,
                           dimnames=list(nms,nms)))}

    FisherInfo <- FisherInfo.fct(param)
    L2Fam <- new("GammaFamily")
    L2Fam@name <- name
    L2Fam@distribution <- distribution
    L2Fam@distrSymm <- distrSymm
    L2Fam@param <- param
    L2Fam@modifyParam <- modifyParam
    L2Fam@props <- props
    L2Fam@L2deriv.fct <- L2deriv.fct
    L2Fam@L2derivSymm <- L2derivSymm
    L2Fam@L2derivDistr <- L2derivDistr
    L2Fam@L2derivDistrSymm <- L2derivDistrSymm
    L2Fam@FisherInfo.fct <- FisherInfo.fct
    L2Fam@FisherInfo <- FisherInfo
    L2Fam@startPar <- startPar
    L2Fam@makeOKPar <- makeOKPar
    L2Fam@scaleshapename <- c("scale"="scale","shape"="shape")
    
    L2deriv <- EuclRandVarList(RealRandVariable(L2deriv.fct(param),
                               Domain = Reals()))

    if(!is.function(trafo))
       f.call <- substitute(GammaFamily(scale = s1, shape = s2,
  	                           trafo = matrix(Tr, ncol = 2, dimnames = DN)),
  	                     list(s1 = scale, s2 = shape, Tr = trafo,
  	                          DN = dimnames(trafo)))
    else
       f.call <- substitute(GammaFamily(scale = s1, shape = s2, trafo = Tr),
  	                     list(s1 = scale, s2 = shape, Tr = trafo))

    L2Fam@fam.call <- f.call

    L2Fam@LogDeriv <- function(x) -(shape-1)/x + 1/scale
    L2Fam@L2deriv <- L2deriv


    return(L2Fam)
}

##################################################################
## Beta family   :: new  08/08 P.R.
##################################################################
BetaFamily <- function(shape1 = 1, shape2 = 1, trafo, withL2derivDistr = TRUE){
    name <- "Beta family"
    distribution <- Beta(shape1=shape1, shape2 = shape2)
    distrSymm <- NoSymmetry()
    param0 <- c(shape1, shape2)
    names(param0) <- nms <- c("shape1", "shape2")
    if(missing(trafo)) {trafo <- diag(2); dimnames(trafo) <-list(nms,nms)}
    param <- ParamFamParameter(name = "shape1 and shape2",  
                        main = param0, trafo = trafo)
    modifyParam <- function(theta){ Beta(shape1 = theta[1], shape2 = theta[2]) }
    makeOKPar <- function(param) {param <- pmax(.Machine$double.eps,param)
                                  return(param)}
    props <- c("The Beta family is invariant in the following sense",
               "if (x_i)~Beta(s1,s2) then (1-x_i)~Beta(s2,s1)")
    startPar <- function(x,...){ x <- pmax(0,pmin(1,x))
                              E1 <- mean(x)
                              V <- var(x)
                              D <- E1*(1-E1)/V-1
                              st <- c(E1*D,(1-E1)*D)
                              names(st) <- nms
                              return(st)
                              }
    L2deriv.fct <- function(param){
                   shape1 <- main(param)[1]
                   shape2 <- main(param)[2]
                   fct1 <- function(x){}
                   fct2 <- function(x){}
                   body(fct1) <- substitute({log(x)-digamma(shape1)+
                                             digamma(shape1+shape2)},
                        list(shape1 = shape1, shape2 = shape2))
                   body(fct2) <- substitute({log(1-x)-digamma(shape2)+
                                             digamma(shape1+shape2)},
                        list(shape1 = shape1, shape2 = shape2))
                   return(list(fct1, fct2))}
    L2derivSymm <- FunSymmList(NonSymmetric(), NonSymmetric())
    L2derivDistr <- NULL
    if(withL2derivDistr)
       L2derivDistr <- UnivarDistrList(log(Beta(shape1 = shape1, shape2 = shape2))-
                                        digamma(shape1)+digamma(shape1+shape2), 
                                    log(Beta(shape1 = shape2, shape2 = shape1))-
                                        digamma(shape2)+digamma(shape1+shape2))
    L2derivDistrSymm <- DistrSymmList(NoSymmetry(), NoSymmetry())
    FisherInfo.fct <- function(param){
                   shape1 <- main(param)[1]
                   shape2 <- main(param)[2]
                   FI <- diag(trigamma(main(param)))-trigamma(sum(main(param))) 
                   dimnames(FI) <- list(nms,nms)
                   PosDefSymmMatrix(FI)}

    FisherInfo <- FisherInfo.fct(param)
    res <- L2ParamFamily(name = name, distribution = distribution, 
        distrSymm = distrSymm, param = param, modifyParam = modifyParam,
        props = props, L2deriv.fct = L2deriv.fct, L2derivSymm = L2derivSymm,
        L2derivDistr = L2derivDistr, L2derivDistrSymm = L2derivDistrSymm,
        FisherInfo.fct = FisherInfo.fct, FisherInfo = FisherInfo,
        startPar = startPar, makeOKPar = makeOKPar, 
        .returnClsName = "BetaFamily")
    if(!is.function(trafo))
       f.call <- substitute(BetaFamily(shape1 = s1, shape2 = s2,
                                   trafo = matrix(Tr, ncol = 2, dimnames = DN)),
                        list(s1 = shape1, s2 = shape2, Tr = trafo,
                             DN = dimnames(trafo)))
    else
       f.call <- substitute(BetaFamily(shape1 = s1, shape2 = s2, trafo = Tr),
                        list(s1 = shape1, s2 = shape2, Tr = trafo))
    res@fam.call <- f.call
    return(res)
}


################################################################################
## Group Models with central distribution Norm(0,1)
################################################################################

##################################################################
## Normal location family
##################################################################
NormLocationFamily <- function(mean = 0, sd = 1, trafo){ 
    if(missing(trafo)) trafo <- matrix(1, dimnames=list("mean","mean"))
    modParam <- function(theta){}
    body(modParam) <- substitute({ Norm(mean = theta, sd = scale) },
                                 list(scale = sd))
    res <- L2LocationFamily(loc = mean, name = "normal location family",
                     locname = c("loc"="mean"),
                     centraldistribution = Norm(mean = 0, sd = sd),
                     modParam = modParam,
                     LogDeriv = function(x) x/sd^2,
                     L2derivDistr.0 = Norm(mean = 0, sd = 1/sd),
                     distrSymm = SphericalSymmetry(SymmCenter = mean),
                     L2derivSymm = FunSymmList(OddSymmetric(SymmCenter = mean)), 
                     L2derivDistrSymm = DistrSymmList(SphericalSymmetry()),
                     FisherInfo.0 = matrix(1/sd^2, dimnames = list("mean","mean")),
                     trafo = trafo, .returnClsName = "NormLocationFamily")
    if(!is.function(trafo))
       f.call <- substitute(NormLocationFamily(mean = m, sd = s,
                                trafo = matrix(Tr, dimnames=list("mean","mean"))),
                         list(m = mean, s = sd, Tr = trafo))
    else
       f.call <- substitute(NormLocationFamily(mean = m, sd = s, trafo = Tr),
                         list(m = mean, s = sd, Tr = trafo))
    res@fam.call <- f.call
    return(res)
}

##################################################################
## Normal scale family
##################################################################
NormScaleFamily <- function(sd = 1, mean = 0, trafo){ 
    if(missing(trafo)) trafo <- matrix(1, dimnames=list("scale","scale"))
    modParam <- function(theta){}
    body(modParam) <- substitute({ Norm(mean = loc, sd = theta) },
                                 list(loc = mean))
    res <- L2ScaleFamily(loc = mean, scale = sd, name = "normal scale family",
                  locscalename = c("loc"="mean", "scale"="sd"), 
                  modParam = modParam, 
                  LogDeriv = function(x) x,
                  L2derivDistr.0 = (Chisq(df = 1, ncp = 0)-1)/sd,
                  distrSymm = SphericalSymmetry(SymmCenter = mean),
                  L2derivSymm = FunSymmList(EvenSymmetric(SymmCenter = mean)),
                  L2derivDistrSymm = DistrSymmList(NoSymmetry()),                  
                  FisherInfo.0 = matrix(2, dimnames = list("sd", "sd")),
                  trafo = trafo, .returnClsName = "NormScaleFamily")
    if(!is.function(trafo))
       f.call <- substitute(NormScaleFamily(sd = s, mean = m,
                            trafo = matrix(Tr, dimnames=list("sd","sd"))),
                         list(s = sd, m = mean, Tr = trafo))
    else
       f.call <- substitute(NormScaleFamily(sd = s, mean = m, trafo = Tr),
                         list(s = sd, m = mean, Tr = trafo))
    res@fam.call <- f.call
    return(res)
}

##################################################################
## Normal location and scale family
##################################################################
NormLocationScaleFamily <- function(mean = 0, sd = 1, trafo){ 
    lsname <- c("loc"="mean", "scale"="sd")
    if(missing(trafo)) {trafo <- diag(2) 
                        dimnames(trafo) <- list(lsname,lsname)}
    res <- L2LocationScaleFamily(loc = mean, scale = sd, 
              name = "normal location and scale family", 
              locscalename = lsname, 
              modParam = function(theta) Norm(mean = theta[1], sd = theta[2]),
              LogDeriv = function(x) x,
              L2derivDistr.0 = list( Norm(mean = 0, sd=1/sd), 
                                    (Chisq(df = 1, ncp = 0)-1)/sd),
              FisherInfo.0 = matrix(c(1,0,0,2),2,2, 
                                           dimnames = list(lsname, lsname)),
              distrSymm = SphericalSymmetry(SymmCenter = mean),
              L2derivSymm = FunSymmList(OddSymmetric(SymmCenter = mean), 
                                        EvenSymmetric(SymmCenter = mean)),
              L2derivDistrSymm = DistrSymmList(SphericalSymmetry(), 
                                               NoSymmetry()),
              trafo = trafo, .returnClsName = "NormLocationScaleFamily")
    if(!is.function(trafo))
       f.call <- substitute(NormLocationScaleFamily(mean = m, sd = s,
  	                               trafo = matrix(Tr, ncol = 2, dimnames = DN)),
  	                   list(m = mean, s = sd, Tr = trafo, DN = dimnames(trafo)))
    else
       f.call <- substitute(NormLocationScaleFamily(mean = m, sd = s, trafo = Tr),
  	                   list(m = mean, s = sd, Tr = trafo))
    res@fam.call <- f.call
    return(res)
}

###############################################################################
## other location and / or scale models
###############################################################################

##################################################################
## Exponential scale family
##################################################################
ExpScaleFamily <- function(scale = 1, trafo){ 
    if(missing(trafo)) trafo <- matrix(1, dimnames = list("scale","scale"))
    res <- L2ScaleFamily(loc = 0, scale = scale, name = "Exponential scale family", 
                  centraldistribution = Exp(rate = 1),
                  locscalename = c("loc"="", "scale"="scale"), 
                  modParam = function(theta) Exp(rate = 1/theta),
                  LogDeriv = function(x) 1,
                  L2derivDistr.0 = (Exp(rate = 1)-1)/scale,
                  FisherInfo.0 = matrix(1, dimnames = list("scale","scale")), 
                  distrSymm = NoSymmetry(), 
                  L2derivSymm = FunSymmList(OddSymmetric(SymmCenter = scale)), 
                  L2derivDistrSymm = DistrSymmList(NoSymmetry()),
                  trafo = trafo, .returnClsName = "ExpScaleFamily")
    if(!is.function(trafo))
       f.call <- substitute(ExpScaleFamily(scale = s,
                         trafo = matrix(Tr, dimnames = list("scale","scale"))),
                         list(s = scale, Tr = trafo))
    else
      f.call <- substitute(ExpScaleFamily(scale = s, trafo = Tr),
                         list(s = scale, Tr = trafo))

    par0 <- res@param
    par0@fixed <- NULL
    res@param <- par0
    
    res@fam.call <- f.call
    return(res)
}


##################################################################
## Lognormal scale family
##################################################################
LnormScaleFamily <- function(meanlog = 0, sdlog = 1, trafo){ 
    if(missing(trafo)) trafo <- matrix(1, dimnames = list("scale","scale"))
    modParam <- function(theta){}
    body(modParam) <- substitute({ Lnorm(meanlog = log(theta), sdlog = sd1) },
                                 list(sd1 = sdlog))
    res <- L2ScaleFamily(loc = 0, scale = exp(meanlog),  
                  name = "lognormal scale family", 
                  locscalename = c("loc"="", "scale"="meanlog"), 
                  centraldistribution = Lnorm(meanlog = 0, sdlog = sdlog),
                  modParam = modParam,
                  LogDeriv = function(x) log(x)/x/sdlog^2 + 1/x,
#                  L2derivDistr.0 = AbscontDistribution(r=function(n){
#                    x <- rlnorm(n); (log(x)-1)/x}),
# wrong in my opinion
# (x/scale*LogDeriv(x/scale) - 1)/scale = (log(x/scale)/sdlog^2 + 1 - 1)/scale
#                                       = log(x/scale)/sdlog^2/scale
#                                       = (log(x) - meanlog)/sdlog/(sdlog*scale)
# now x ~ Lnorm(meanlog, sdlog)
# => log(x) ~ Norm(meanlog, sdlog^2)
# => (log(x) - meanlog)/sdlog ~ Norm(0, 1)
                  L2derivDistr.0 = Norm(mean=0, sd=1/sdlog/exp(meanlog)),
                  FisherInfo.0 = matrix(1/exp(2*meanlog)/sdlog^2,
                                        dimnames = list("scale","scale")), 
                  distrSymm = NoSymmetry(), 
                  L2derivSymm = FunSymmList(NonSymmetric()), 
                  L2derivDistrSymm = DistrSymmList(SphericalSymmetry(SymmCenter = 0)),
                  trafo = trafo, .returnClsName = "LnormScaleFamily")
    if(!is.function(trafo))
       f.call <- substitute(LnormScaleFamily(meanlog = m, sdlog = s,
                          trafo = matrix(Tr, dimnames = list("scale","scale"))),
                         list(m = meanlog, s = sdlog, Tr = trafo))
    else   
       f.call <- substitute(LnormScaleFamily(meanlog = m, sdlog = s, trafo = Tr),
                         list(m = meanlog, s = sdlog, Tr = trafo))
    res@fam.call <- f.call
    return(res)
}




##################################################################
## Cauchy location scale family
##################################################################
CauchyLocationScaleFamily <- function(loc = 0, scale = 1, trafo){ 
    if(missing(trafo)) {trafo <- diag(2) 
                        dimnames(trafo) <- list(c("loc","scale"),
                                                c("loc","scale"))}
    res <- L2LocationScaleFamily(loc = loc, scale = scale, 
                  name = "Cauchy Location and scale family", 
                  centraldistribution = Cauchy(),
                  LogDeriv = function(x)  2*x/(x^2+1),  
                  distrSymm = SphericalSymmetry(SymmCenter = loc),
                  L2derivSymm = FunSymmList(OddSymmetric(SymmCenter = loc), 
                                            EvenSymmetric(SymmCenter = loc)),
                  L2derivDistrSymm = DistrSymmList(SphericalSymmetry(), 
                                                   NoSymmetry()),
                  L2derivDistr.0 = UnivarDistrList(Arcsine(),abs(Arcsine())),
                  FisherInfo.0 = matrix(c(1,0,0,1)/2,2,2, 
                                           dimnames = list(c("loc","scale"),
                                                           c("loc","scale"))),
                  trafo = trafo, .returnClsName = "CauchyLocationScaleFamily")
    
    if(!is.function(trafo))
       f.call <- substitute(CauchyLocationScaleFamily(loc = l, scale = s,
  	                         trafo = matrix(Tr, ncol = 2, dimnames = DN)),
  	                 list(l = loc, s = scale, Tr = trafo, DN = dimnames(trafo)))
    else
       f.call <- substitute(CauchyLocationScaleFamily(loc = l, scale = s, 
                            trafo = Tr), list(l = loc, s = scale, Tr = trafo))
    res@fam.call <- f.call
    return(res)
}


#####################################
#####################################
#### normal models with nuisance
#####################################
#####################################

##################################################################
## Normal location family  with unknown scale
##################################################################
NormLocationUnknownScaleFamily <- function(mean = 0, sd = 1, trafo){ 
    if(missing(trafo)) {trafo <- diag(1) 
                        dimnames(trafo) <- list("mean","mean")}
    lsname <- c("loc"="mean", "scale"="sd")
    res <- L2LocationUnknownScaleFamily(loc = mean, scale = sd, 
                     name = "normal location family with unknown scale (as nuisance)",
                     locscalename = lsname, 
                     modParam = function(theta) Norm(mean = theta[1], sd = theta[2]),
                     LogDeriv = function(x) x,
                     L2derivDistr.0 = list( Norm(mean = 0, sd=1/sd), 
                                    (Chisq(df = 1, ncp = 0)-1)/sd),
                     FisherInfo.0 = matrix(c(1,0,0,2),2,2, 
                                           dimnames = list(lsname,
                                                           lsname)),
                     distrSymm = SphericalSymmetry(SymmCenter = mean),
                     L2derivSymm = FunSymmList(OddSymmetric(SymmCenter = mean), 
                                               EvenSymmetric(SymmCenter = mean)),
                     L2derivDistrSymm = DistrSymmList(SphericalSymmetry(), 
                                                      NoSymmetry()),
                     trafo = trafo)
    if(!is.function(trafo))
       f.call <- substitute(NormLocationUnknownScaleFamily(mean = m, sd = s,
                              trafo = matrix(Tr, dimnames = list("mean","mean"))),
                         list(m = mean, s = sd, Tr = trafo))
    else
       f.call <- substitute(NormLocationUnknownScaleFamily(mean = m, sd = s,
                              trafo = Tr), list(m = mean, s = sd, Tr = trafo))
    res@fam.call <- f.call
    return(res)
}

##################################################################
## Normal scale family  with unknown location
##################################################################
NormScaleUnknownLocationFamily <- function(sd = 1, mean = 0, trafo){ 
    if(missing(trafo)) {trafo <- diag(1) 
                        dimnames(trafo) <- list("sd","sd")}
    lsname <- c("loc"="mean", "scale"="sd")
    res <- L2ScaleUnknownLocationFamily(loc = mean, scale = sd, 
                  name = "normal scale family with unknown location (as nuisance)", 
                  locscalename = lsname, 
                  modParam = function(theta) Norm(mean = theta[2], sd = theta[1]),
                  LogDeriv = function(x) x,
                  L2derivDistr.0 = list( Norm(mean = 0, sd=1/sd), 
                                       (Chisq(df = 1, ncp = 0)-1)/sd),
                  FisherInfo.0 = matrix(c(1,0,0,2),2,2, 
                                        dimnames = list(lsname,
                                                        lsname)),
                  distrSymm = SphericalSymmetry(SymmCenter = mean),
                  L2derivSymm = FunSymmList(OddSymmetric(SymmCenter = mean), 
                                            EvenSymmetric(SymmCenter = mean)),
                  L2derivDistrSymm = DistrSymmList(SphericalSymmetry(), 
                                                   NoSymmetry()),
                  trafo = trafo)
    if(!is.function(trafo))
       f.call <- substitute(NormScaleUnknownLocationFamily(sd = s, mean = m,
                              trafo = matrix(Tr, dimnames = list("sd","sd"))),
                         list(m = mean, s = sd, Tr = trafo))
    else
       f.call <- substitute(NormScaleUnknownScaleFamily(sd = s, mean = m, 
                              trafo = Tr), list(m = mean, s = sd, Tr = trafo))
    res@fam.call <- f.call
    return(res)
}
