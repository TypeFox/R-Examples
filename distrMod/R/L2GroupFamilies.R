##################################################################
## L2 location family
##################################################################
L2LocationFamily <- function(loc = 0, name, centraldistribution = Norm(),
                             locname = "loc", modParam,
                             LogDeriv, L2derivDistr.0,
                             FisherInfo.0, 
                             distrSymm, L2derivSymm, L2derivDistrSymm,
                             trafo, .returnClsName = NULL){
    if(missing(name))
       name <- "L2 location family"

    if(!length(locname)==1) stop("argument 'locname' must be of length 1.")
    names(locname) <- "loc"
    
    distribution <- centraldistribution + loc

    if(missing(distrSymm)){
        distrSymm <- SphericalSymmetry(SymmCenter = loc)
    }else{
        if(!is(distrSymm, "NoSymmetry")){
            if(!is(distrSymm@SymmCenter, "numeric"))
                stop("slot 'SymmCenter' of 'distrSymm' has to be of class 'numeric'")
            if(length(distrSymm@SymmCenter) != 1)
                stop("slot 'SymmCenter' of 'distrSymm' has wrong dimension")
        }
    }

    makeOKPar <- function(param) param
    startPar <- function(x,...) c(min(x),max(x))
    param0 <- loc
    names(param0) <- locname
    if(missing(trafo))  {trafo <- matrix(1)
                         dimnames(trafo) <- list(locname,locname)}
    param <- ParamFamParameter(name = "loc", main = param0, trafo = trafo)
    if(missing(modParam))
        modParam <- function(theta){ centraldistribution + theta }
    props <- c(paste("The", name, "is invariant under"),
               "the group of transformations 'g(x) = x + loc'",
               "with location parameter 'loc'")
    if(missing(LogDeriv)) LogDeriv <- .getLogDeriv(centraldistribution)
    L2deriv.fct <- function(param){
                   loc <- main(param)
                   fct <- function(x){}
                   body(fct) <- substitute({ LogDeriv(x - loc) }, list(loc = loc))
                   return(fct)}
    L2deriv <- EuclRandVarList(RealRandVariable(list(L2deriv.fct(param)), Domain = Reals())) 


    L2derivDistr <- if(missing(L2derivDistr.0))
        imageDistr(RandVar = L2deriv, distr = distribution) else 
        UnivarDistrList(L2derivDistr.0)

    if(missing (L2derivDistrSymm)){ 
        L2derivDistrSymm <- DistrSymmList(NoSymmetry())
    }else{
        if(!length(L2derivDistrSymm) == 1) 
            stop("wrong length of argument L2derivDistrSymm")
    }
    if(missing (L2derivSymm)){ 
        L2derivSymm <- FunSymmList(NonSymmetric())
    }else{
        if(!length(L2derivSymm) == 1) 
            stop("wrong length of argument L2derivSymm")
    }

    FI0 <- if(missing(FisherInfo.0))
           E(centraldistribution, fun = function(x) LogDeriv(x)^2,
             useApply = FALSE) else FisherInfo.0

    FI0 <- matrix(FI0,1,1, dimnames = list(locname,locname))
    FisherInfo.fct <- function(param) PosDefSymmMatrix(FI0)

    if(is.function(trafo))  
        Tr0 <- trafo
    else if(is.matrix(trafo)){
            Tr0 <- trafo
            if(is.null(dimnames(trafo)))
               dimnames(trafo) <- list(locname,locname) 
         }else{
            Tr0 <- matrix(trafo, dimnames = list(locname,locname)) 
         }                    
    
    f.call <- substitute(L2LocationFamily(loc = l,
                  name = N,
                  centraldistribution = D0,
                  locname = lN,
                  modParam = mP,
                  LogDeriv = lD,
                  L2derivDistr.0 = L2D0,
                  FisherInfo.0 = F.0,
                  distrSymm = DSymm,
                  L2derivSymm = L2Symm,
                  L2derivDistrSymm = L2DSymm,
                  trafo = Tr,
                  .returnClsName = rtn),
             list(l = loc,
                  N = name,
                  D0 = centraldistribution,
                  lN = locname,
                  mP = modParam,
                  lD = LogDeriv,
                  L2D0 = L2derivDistr[[1]],
                  F.0 = FI0,
                  DSymm = distrSymm,
                  L2Symm = L2derivSymm,
                  L2DSymm = L2derivDistrSymm,
                  Tr = Tr0,   
                  rtn = .returnClsName))

    if(is.null(.returnClsName))
       .returnClsName <- "L2LocationFamily"
    L2Fam <- new(.returnClsName)
    L2Fam@name <- name
    L2Fam@locscalename <- locname
    L2Fam@distribution <- distribution
    L2Fam@distrSymm <- distrSymm
    L2Fam@param <- param
    L2Fam@startPar <- startPar
    L2Fam@makeOKPar <- makeOKPar
    L2Fam@modifyParam <- modParam
    L2Fam@fam.call <- f.call
    L2Fam@props <- props
    L2Fam@LogDeriv <- LogDeriv
    L2Fam@L2deriv.fct <- L2deriv.fct
    L2Fam@L2deriv <- L2deriv
    L2Fam@L2derivSymm <- L2derivSymm
    L2Fam@L2derivDistr <- L2derivDistr
    L2Fam@L2derivDistrSymm <- L2derivDistrSymm
    L2Fam@FisherInfo.fct <- FisherInfo.fct
    L2Fam@FisherInfo <- FisherInfo.fct(param)

    return(L2Fam)
}


##################################################################
## L2 scale family
##################################################################
L2ScaleFamily <- function(scale = 1, loc = 0, name, centraldistribution = Norm(),
                          locscalename = c("loc", "scale"), modParam,
                          LogDeriv, L2derivDistr.0,
                          FisherInfo.0,
                          distrSymm, L2derivSymm, L2derivDistrSymm,
                          trafo, .returnClsName = NULL){
    if(length(scale) != 1 || !is.numeric(scale))
        stop("scale has to be a numeric of length 1")
    if(scale < 0)
        stop("scale has to be positive")
    if(missing(name))
       name <- "L2 scale family"
    if((length(locscalename)<1)||(length(locscalename)>2)) 
        stop("argument 'locscalename' must be of length 1 or 2.")
    if(length(locscalename)==1){ 
        locscalename <- c(locscalename, "loc")
        names(locscalename) <- c("scale", "loc")
    }else{
        if(!all(names(locscalename)%in% c("loc","scale")) || 
            is.null(names(locscalename)))
           names(locscalename) <- c("loc", "scale")   
    }
    scalename <- locscalename["scale"]
    distribution <- scale*centraldistribution + loc

    if(missing(distrSymm)){
        distrSymm <- SphericalSymmetry(SymmCenter = loc)
    }else{
        if(!is(distrSymm, "NoSymmetry")){
            if(!is(distrSymm@SymmCenter, "numeric"))
                stop("slot 'SymmCenter' of 'distrSymm' has to be of class 'numeric'")
            if(length(distrSymm@SymmCenter) != 1)
                stop("slot 'SymmCenter' of 'distrSymm' has wrong dimension")
        }
    }

    param0 <- scale
    names(param0) <- locscalename["scale"]
    param1 <- loc
    names(param1) <- locscalename["loc"]
    
    startPar <- function(x,...) c(.Machine$double.eps,max(x)-min(x))
    makeOKPar <- function(param) abs(param)+.Machine$double.eps
    if(missing(trafo))  {trafo <- matrix(1)
                         dimnames(trafo) <- list(scalename,scalename)}
    param <- ParamFamParameter(name = "scale", main = param0, 
                               fixed = param1, trafo = trafo,
                               .returnClsName ="ParamWithScaleFamParameter")
    if(missing(modParam)){
        modParam <- function(theta){}
        body(modParam) <- substitute({ theta*centraldistribution+loc },
                                        list(loc = loc))
    }
    props <- c(paste("The", name, "is invariant under"),
               "the group of transformations 'g(y) = scale*y'",
               "with scale parameter 'scale'")
    if(missing(LogDeriv)) LogDeriv <- .getLogDeriv(centraldistribution)
    L2deriv.fct <- function(param){
                   scale <- main(param)
                   fct <- function(x){}
                   body(fct) <- substitute({ ((x - loc)/scale*LogDeriv((x - loc)/scale)-1)/scale },
                                             list(loc = loc, scale = scale))
                   return(fct)}
    L2deriv <- EuclRandVarList(RealRandVariable(list(L2deriv.fct(param)), Domain = Reals())) 


    L2derivDistr <- if(missing(L2derivDistr.0))
        imageDistr(RandVar = L2deriv, distr = distribution) else 
        UnivarDistrList(L2derivDistr.0)

    if(missing (L2derivDistrSymm)){ 
        L2derivDistrSymm <- DistrSymmList(NoSymmetry())
    }else{
        if(!length(L2derivDistrSymm) == 1) 
            stop("wrong length of argument L2derivDistrSymm")
    }
    if(missing (L2derivSymm)){ 
        L2derivSymm <- FunSymmList(NonSymmetric())
    }else{
        if(!length(L2derivSymm) == 1) 
            stop("wrong length of argument L2derivSymm")
    }

    FI0 <- if(missing(FisherInfo.0)) 
           E(centraldistribution, fun = function(x) (x*LogDeriv(x)-1)^2,
             useApply = FALSE) else FisherInfo.0

    FI0 <- matrix(FI0,1,1, dimnames = list(scalename,scalename))

    FisherInfo.fct <- function(param){
                   scale <- main(param)
                   PosDefSymmMatrix(FI0/scale^2)}

    if(is.function(trafo))  
        Tr0 <- trafo
    else if(is.matrix(trafo)){
            Tr0 <- trafo
            if(is.null(dimnames(trafo)))
               dimnames(trafo) <- list(scalename,scalename) 

         }else{
            Tr0 <- matrix(trafo, dimnames = list(scalename,scalename) ) 
         }                    

    f.call <- substitute(L2ScaleFamily(scale = s,
                           loc = l,
                           name = N,
                           centraldistribution = D0,
                           locscalename = lN,
                           modParam = mP,
                           LogDeriv = lD,
                           L2derivDistr.0 = L2D0,
                           FisherInfo.0 = F.0,
                           distrSymm = DSymm,
                           L2derivSymm = L2Symm,
                           L2derivDistrSymm = L2DSymm,
                           trafo = Tr,
                           .returnClsName = rtn),
                      list(s = scale,
                           l = loc,
                           N = name,
                           D0 = centraldistribution,
                           lN = locscalename,
                           mP = modParam,
                           lD = LogDeriv,
                           L2D0 = L2derivDistr[[1]],
                           F.0 = FI0,
                           DSymm = distrSymm,
                           L2Symm = L2derivSymm,
                           L2DSymm = L2derivDistrSymm,
                           Tr = Tr0,   
                           rtn = .returnClsName))

    if(is.null(.returnClsName))
       .returnClsName <- "L2ScaleFamily"
    L2Fam <- new(.returnClsName)
    L2Fam@name <- name
    L2Fam@locscalename <- locscalename
    L2Fam@distribution <- distribution
    L2Fam@distrSymm <- distrSymm
    L2Fam@param <- param
    L2Fam@modifyParam <- modParam
    L2Fam@startPar <- startPar
    L2Fam@makeOKPar <- makeOKPar
    L2Fam@fam.call <- f.call
    L2Fam@props <- props
    L2Fam@LogDeriv <- LogDeriv
    L2Fam@L2deriv.fct <- L2deriv.fct
    L2Fam@L2deriv <- L2deriv
    L2Fam@L2derivSymm <- L2derivSymm
    L2Fam@L2derivDistr <- L2derivDistr
    L2Fam@L2derivDistrSymm <- L2derivDistrSymm
    L2Fam@FisherInfo.fct <- FisherInfo.fct
    L2Fam@FisherInfo <- FisherInfo.fct(param)

    return(L2Fam)
}


##################################################################
## L2 location and scale family
##################################################################
L2LocationScaleFamily <- function(loc = 0, scale = 1, name, 
                             centraldistribution = Norm(),
                             locscalename = c("loc", "scale"), modParam,
                             LogDeriv, L2derivDistr.0,
                             FisherInfo.0, 
                             distrSymm, L2derivSymm, L2derivDistrSymm,
                             trafo, .returnClsName = NULL){
    if(length(scale) != 1 || !is.numeric(scale))
        stop("scale has to be a numeric of length 1")
    if(scale < 0)
        stop("scale has to be positive")
    if(missing(name))
       name <- "L2 location and scale family"

    distribution <- scale*centraldistribution+loc

    if(!length(locscalename)==2) 
        stop("argument 'locscalename' must be of length 2.")
    if(!all(names(locscalename)%in% c("loc","scale")) || 
        is.null(names(locscalename)))
           names(locscalename) <- c("loc", "scale")   

    if(missing(distrSymm)){
        distrSymm <- SphericalSymmetry(SymmCenter = loc)
    }else{
        if(!is(distrSymm, "NoSymmetry")){
            if(!is(distrSymm@SymmCenter, "numeric"))
                stop("slot 'SymmCenter' of 'distrSymm' has to be of class 'numeric'")
            if(length(distrSymm@SymmCenter) != 1)
                stop("slot 'SymmCenter' of 'distrSymm' has wrong dimension")
        }
    }

    mad.const <- 1/ if (is(distrSymm, "NoSymmetry")) 
                        mad(centraldistribution) else q(centraldistribution)(.75)
    
    param0 <- c(loc, scale)
    names(param0) <- locscalename
    if(missing(trafo))  {trafo <- diag(2)
                         dimnames(trafo) <- list(locscalename,
                                                 locscalename)}
    param <- ParamFamParameter(name = "location and scale", main = param0,
                               trafo = trafo,
                               .returnClsName ="ParamWithScaleFamParameter")
    
    startPar <- function(x,...) {
                   st <- c(median(x),mad(x, constant=mad.const))
                   names(st) <- locscalename
                   return(st)}
    makeOKPar <- function(param) {
                    st <- c(param[1],abs(param[2])+.Machine$double.eps)
                    names(st) <- locscalename
                   return(st)}
    if(missing(modParam))
        modParam <- function(theta){theta[2]*centraldistribution+theta[1] }
    props <- c(paste("The", name, "is invariant under"),
               "the group of transformations 'g(x) = scale*x + loc'",
               "with location parameter 'loc' and scale parameter 'scale'")

    if(missing(LogDeriv)) LogDeriv <- .getLogDeriv(centraldistribution)
    L2deriv.fct <- function(param){
                   nmsL <- names(main(param))
                   lnm <- if(locscalename["loc"] %in% nmsL)
                             locscalename["loc"] else 1
                   snm <- if(locscalename["scale"] %in% nmsL)
                             locscalename["scale"] else 2          
                   mean <- main(param)[lnm]
                   sd <-   main(param)[snm]
                   fct1 <- function(x){}
                   fct2 <- function(x){}
                   body(fct1) <- substitute({ LogDeriv((x - loc)/scale)/scale },
                                             list(loc = mean, scale = sd))
                   body(fct2) <- substitute({ 
                        ((x - loc)/scale * LogDeriv((x - loc)/scale)-1)/scale },
                                             list(loc = mean, scale = sd))
                   return(list(fct1, fct2))}

    L2deriv <- EuclRandVarList(RealRandVariable(L2deriv.fct(param), 
                               Domain = Reals())) 

    L2derivDistr <- if (missing(L2derivDistr.0))
         imageDistr(RandVar = L2deriv, distr = distribution) else
         UnivarDistrList(L2derivDistr.0[[1]],L2derivDistr.0[[2]])

    if(missing (L2derivDistrSymm)){ 
        L2derivDistrSymm <- DistrSymmList(NoSymmetry(),NoSymmetry())
    }else{
        if(!length(L2derivDistrSymm) == 2) 
            stop("wrong length of argument L2derivDistrSymm")
    }
    if(missing (L2derivSymm)){ 
        L2derivSymm <- FunSymmList(NonSymmetric(),NonSymmetric())
    }else{
        if(!length(L2derivSymm) == 2) 
            stop("wrong length of argument L2derivSymm")
    }

    if(missing(FisherInfo.0)){
        FI11 <- E(centraldistribution, fun = function(x) LogDeriv(x)^2,
                  useApply = FALSE)
        FI22 <- E(centraldistribution, fun = function(x) (x*LogDeriv(x)-1)^2,
               useApply = FALSE)
        if( is(distrSymm, "SphericalSymmetry") ){
            FI12 <- 0
        }else{
            FI12 <- E(centraldistribution, fun = function(x) x*LogDeriv(x)^2,
                      useApply = FALSE)
        }
        FI0 <- matrix(c(FI11,FI12,FI12,FI22),2,2)
    }else{ 
        FI0 <- FisherInfo.0 
    }

    FI0 <- matrix(FI0,2,2,dimnames=list(locscalename,locscalename))

    FisherInfo.fct <- function(param){
                   nmsI <- names(main(param))
                   if(locscalename["scale"] %in% nmsI)
                       scale <- main(param)[locscalename["scale"]]
                   else
                       scale <- main(param)[2]
                   PosDefSymmMatrix(FI0/scale^2)}

    if(is.function(trafo))  
        Tr0 <- trafo
    else if(is.matrix(trafo)){
            Tr0 <- trafo
            if(is.null(dimnames(trafo)))
               dimnames(trafo) <- list(locscalename,locscalename) 

         }else{
            Tr0 <- matrix(trafo, dimnames = list(locscalename,locscalename) ) 
         }                    
f.call <- substitute(L2LocationScaleFamily(loc = l,
               scale = s,
               name = N,
               centraldistribution = D0,
               locscalename = lN,
               modParam = mP,
               LogDeriv = lD,
               L2derivDistr.0 = L2D0,
               FisherInfo.0 = F.0,
               distrSymm = DSymm,
               L2derivSymm = L2Symm,
               L2derivDistrSymm = L2DSymm,
               trafo = Tr,
               .returnClsName = rtn,
               .withMDE = FALSE),
           list(s = scale,
               l = loc,
               N = name,
               D0 = centraldistribution,
               lN = locscalename,
               mP = modParam,
               lD = LogDeriv,
               L2D0 = as(L2derivDistr, "list"),
               F.0 = FI0,
               DSymm = distrSymm,
               L2Symm = L2derivSymm,
               L2DSymm = L2derivDistrSymm,
               Tr = Tr0,
               rtn = .returnClsName))

    if(is.null(.returnClsName))
       .returnClsName <- "L2LocationScaleFamily"
    L2Fam <- new(.returnClsName)
    L2Fam@name <- name
    L2Fam@locscalename <- locscalename
    L2Fam@distribution <- distribution
    L2Fam@distrSymm <- distrSymm
    L2Fam@param <- param
    L2Fam@startPar <- startPar
    L2Fam@makeOKPar <- makeOKPar
    L2Fam@modifyParam <- modParam
    L2Fam@fam.call <- f.call
    L2Fam@props <- props
    L2Fam@LogDeriv <- LogDeriv
    L2Fam@L2deriv.fct <- L2deriv.fct
    L2Fam@L2deriv <- L2deriv
    L2Fam@L2derivSymm <- L2derivSymm
    L2Fam@L2derivDistr <- L2derivDistr
    L2Fam@L2derivDistrSymm <- L2derivDistrSymm
    L2Fam@FisherInfo.fct <- FisherInfo.fct
    L2Fam@FisherInfo <- FisherInfo.fct(param)
    L2Fam@.withMDE <- FALSE
    return(L2Fam)
}

##################################################################
## L2 location with unknown scale (as nuisance) family
##################################################################
L2LocationUnknownScaleFamily <- function(loc = 0, scale = 1, name, 
                             centraldistribution = Norm(),
                             locscalename = c("loc", "scale"), modParam,
                             LogDeriv, L2derivDistr.0,
                             FisherInfo.0, 
                             distrSymm, L2derivSymm, L2derivDistrSymm,
                             trafo, .returnClsName = NULL){
    if(length(scale) != 1 || !is.numeric(scale))
        stop("scale has to be a numeric of length 1")
    if(scale < 0)
        stop("scale has to be positive")
    if(missing(name))
       name <- "L2 location with unknown scale (as nuisance) family"

    if(!length(locscalename)==2) 
        stop("argument 'locscalename' must be of length 2.")
    if(!all(names(locscalename)%in% c("loc","scale")) || 
        is.null(names(locscalename)))
           names(locscalename) <- c("loc", "scale")   

    distribution <- scale*centraldistribution+loc

    if(missing(distrSymm)){
        distrSymm <- SphericalSymmetry(SymmCenter = loc)
    }else{
        if(!is(distrSymm, "NoSymmetry")){
            if(!is(distrSymm@SymmCenter, "numeric"))
                stop("slot 'SymmCenter' of 'distrSymm' has to be of class 'numeric'")
            if(length(distrSymm@SymmCenter) != 1)
                stop("slot 'SymmCenter' of 'distrSymm' has wrong dimension")
        }
    }

    param0 <- c(loc, scale)
    names(param0) <- locscalename
    startPar <- function(x,...) {
                   st <- c(median(x),mad(x))
                   names(st) <- locscalename
                   return(st)}
    makeOKPar <- function(param) {
                    st <- c(param[1],abs(param[2])+.Machine$double.eps)
                    names(st) <- locscalename
                   return(st)}
    if(missing(trafo))  {trafo <- matrix(1)
                         dimnames(trafo) <- list(locscalename["loc"],
                                                 locscalename["loc"])}
    param <- ParamFamParameter(name = "location and scale", main = param0[1],
                               nuisance = param0[2], trafo = trafo,
                               .returnClsName ="ParamWithScaleFamParameter")
    if(missing(modParam))
        modParam <- function(theta){theta[2]*centraldistribution+theta[1] }
    props <- c(paste("The", name, "is invariant under"),
               "the group of transformations 'g(x) = scale*x + loc'",
               "with location parameter 'loc' and scale parameter 'scale'")

    if(missing(LogDeriv)) LogDeriv <- .getLogDeriv(centraldistribution)
    L2deriv.fct <- function(param){
                   mean <- main(param)
                   sd <-   nuisance(param)
                   fct1 <- function(x){}
                   fct2 <- function(x){}
                   body(fct1) <- substitute({ LogDeriv((x - loc)/scale)/scale },
                                             list(loc = mean, scale = sd))
                   body(fct2) <- substitute({ 
                        ((x - loc)/scale * LogDeriv((x - loc)/scale)-1)/scale },
                                             list(loc = mean, scale = sd))
                   return(list(fct1, fct2))}

    L2deriv <- EuclRandVarList(RealRandVariable(L2deriv.fct(param), 
                               Domain = Reals())) 


    L2derivDistr <- if (missing(L2derivDistr.0))
         imageDistr(RandVar = L2deriv, distr = distribution) else
         UnivarDistrList(L2derivDistr.0[[1]],L2derivDistr.0[[2]])

    if(missing (L2derivDistrSymm)){ 
        L2derivDistrSymm <- DistrSymmList(NoSymmetry(),NoSymmetry())
    }else{
        if(!length(L2derivDistrSymm) == 2) 
            stop("wrong length of argument L2derivDistrSymm")
    }
    if(missing (L2derivSymm)){ 
        L2derivSymm <- FunSymmList(NonSymmetric(),NonSymmetric())
    }else{
        if(!length(L2derivSymm) == 2) 
            stop("wrong length of argument L2derivSymm")
    }

    if(missing(FisherInfo.0)){
        FI11 <- E(centraldistribution, fun = function(x) LogDeriv(x)^2,
                  useApply = FALSE)
        FI22 <- E(centraldistribution, fun = function(x) (x*LogDeriv(x)-1)^2,
               useApply = FALSE)
        if( is(distrSymm, "SphericalSymmetry") ){
            FI12 <- 0
        }else{
            FI12 <- E(centraldistribution, fun = function(x) x*LogDeriv(x)^2,
                      useApply = FALSE)
        }
        FI0 <- matrix(c(FI11,FI12,FI12,FI22),2,2)
    }else{ 
        FI0 <- FisherInfo.0 
    }
    FI0 <- matrix(FI0,2,2,dimnames=list(names(param0),names(param0)))

    FisherInfo.fct <- function(param){
                   scale <- nuisance(param)
                   PosDefSymmMatrix(FI0/scale^2)}

    if(is.function(trafo))  
        Tr0 <- trafo
    else if(is.matrix(trafo)){
            Tr0 <- trafo
            if(is.null(dimnames(trafo)))
               dimnames(trafo) <- list(names(param0),names(param0)) 

         }else{
            Tr0 <- matrix(trafo, dimnames = list(names(param0),names(param0)) ) 
         }                    
    f.call <- substitute(L2LocationUnknownScaleFamily(loc = l,
                 scale = s,
                 name = N,
                 centraldistribution = D0,
                 locscalename = lN,
                 modParam = mP,
                 LogDeriv = lD,
                 L2derivDistr.0 = L2D0,
                 FisherInfo.0 = F.0,
                 distrSymm = DSymm,
                 L2derivSymm = L2Symm,
                 L2derivDistrSymm = L2DSymm,
                 trafo = Tr,                 
                 .returnClsName = rtn),
             list(s = scale,
                 l = loc,
                 N = name,
                 D0 = centraldistribution,
                 lN = locscalename,
                 mP = modParam,
                 lD = LogDeriv,
                 L2D0 = as(L2derivDistr, "list"),
                 F.0 = FI0,
                 DSymm = distrSymm,
                 L2Symm = L2derivSymm,
                 L2DSymm = L2derivDistrSymm,
                 Tr = Tr0,   
                 rtn = .returnClsName))

    if(is.null(.returnClsName))
       .returnClsName <- "L2LocationScaleFamily"
    L2Fam <- new(.returnClsName)
    L2Fam@name <- name
    L2Fam@locscalename <- locscalename
    L2Fam@distribution <- distribution
    L2Fam@distrSymm <- distrSymm
    L2Fam@param <- param
    L2Fam@modifyParam <- modParam
    L2Fam@startPar <- startPar
    L2Fam@makeOKPar <- makeOKPar
    L2Fam@fam.call <- f.call
    L2Fam@props <- props
    L2Fam@LogDeriv <- LogDeriv
    L2Fam@L2deriv.fct <- L2deriv.fct
    L2Fam@L2deriv <- L2deriv
    L2Fam@L2derivSymm <- L2derivSymm
    L2Fam@L2derivDistr <- L2derivDistr
    L2Fam@L2derivDistrSymm <- L2derivDistrSymm
    L2Fam@FisherInfo.fct <- FisherInfo.fct
    L2Fam@FisherInfo <- FisherInfo.fct(param)

    return(L2Fam)
}

##################################################################
## L2 scale with unknown location (as nuisance) family
##################################################################
L2ScaleUnknownLocationFamily <- function(loc = 0, scale = 1, name, 
                             centraldistribution = Norm(),
                             locscalename = c("loc", "scale"), modParam,
                             LogDeriv, L2derivDistr.0,
                             FisherInfo.0, 
                             distrSymm, L2derivSymm, L2derivDistrSymm,
                             trafo, .returnClsName = NULL){
    if(length(scale) != 1 || !is.numeric(scale))
        stop("scale has to be a numeric of length 1")
    if(scale < 0)
        stop("scale has to be positive")
    if(missing(name))
       name <- "L2 scale with unknown location (as nuisance) family"

    if(!length(locscalename)==2) 
        stop("argument 'locscalename' must be of length 2.")
    if(!all(names(locscalename)%in% c("loc","scale")) || 
        is.null(names(locscalename)))
           names(locscalename) <- c("loc", "scale")   

    distribution <- scale*centraldistribution+loc

    if(missing(distrSymm)){
        distrSymm <- SphericalSymmetry(SymmCenter = loc)
    }else{
        if(!is(distrSymm, "NoSymmetry")){
            if(!is(distrSymm@SymmCenter, "numeric"))
                stop("slot 'SymmCenter' of 'distrSymm' has to be of class 'numeric'")
            if(length(distrSymm@SymmCenter) != 1)
                stop("slot 'SymmCenter' of 'distrSymm' has wrong dimension")
        }
    }

    param0 <- c(scale, loc)
    names(param0) <- locscalename[c("scale", "loc")]
    startPar <- function(x,...) {
                   st <- c(median(x),mad(x))
                   names(st) <- locscalename
                   return(st)}
    makeOKPar <- function(param) {
                    st <- c(param[locscalename["loc"]],
                            abs(param[locscalename["scale"]])+.Machine$double.eps)
                    names(st) <- locscalename
                   return(st)}
    if(missing(trafo))  {trafo <- matrix(1)
                         dimnames(trafo) <- list("scale","scale")}
    param <- ParamFamParameter(name = "scale and location", main = param0[1],
                               nuisance = param0[2], trafo = trafo,
                               .returnClsName ="ParamWithScaleFamParameter")
    if(missing(modParam))
        modParam <- function(theta){theta[1]*centraldistribution+theta[2] }
    props <- c(paste("The", name, "is invariant under"),
               "the group of transformations 'g(x) = scale*x + loc'",
               "with location parameter 'loc' and scale parameter 'scale'")

    if(missing(LogDeriv)) LogDeriv <- .getLogDeriv(centraldistribution)
    L2deriv.fct <- function(param){
                   mean <- nuisance(param)
                   sd <-   main(param)
                   fct1 <- function(x){}
                   fct2 <- function(x){}
                   body(fct1) <- substitute({ LogDeriv((x - loc)/scale)/scale },
                                             list(loc = mean, scale = sd))
                   body(fct2) <- substitute({ 
                        ((x - loc)/scale * LogDeriv((x - loc)/scale)-1)/scale },
                                             list(loc = mean, scale = sd))
                   return(list(fct1, fct2))}

    L2deriv <- EuclRandVarList(RealRandVariable(L2deriv.fct(param), 
                               Domain = Reals())) 


    L2derivDistr <- if (missing(L2derivDistr.0))
         imageDistr(RandVar = L2deriv, distr = distribution) else
         UnivarDistrList(L2derivDistr.0[[1]],L2derivDistr.0[[2]])

    if(missing (L2derivDistrSymm)){ 
        L2derivDistrSymm <- DistrSymmList(NoSymmetry(),NoSymmetry())
    }else{
        if(!length(L2derivDistrSymm) == 2) 
            stop("wrong length of argument L2derivDistrSymm")
    }
    if(missing (L2derivSymm)){ 
        L2derivSymm <- FunSymmList(NonSymmetric(),NonSymmetric())
    }else{
        if(!length(L2derivSymm) == 2) 
            stop("wrong length of argument L2derivSymm")
    }

    if(missing(FisherInfo.0)){
        FI11 <- E(centraldistribution, fun = function(x) LogDeriv(x)^2,
                  useApply = FALSE)
        FI22 <- E(centraldistribution, fun = function(x) (x*LogDeriv(x)-1)^2,
               useApply = FALSE)
        if( is(distrSymm, "SphericalSymmetry") ){
            FI12 <- 0
        }else{
            FI12 <- E(centraldistribution, fun = function(x) x*LogDeriv(x)^2,
                      useApply = FALSE)
        }
        FI0 <- matrix(c(FI11,FI12,FI12,FI22),2,2)
    }else{ 
        FI0 <- FisherInfo.0 
    }
    FI0 <- matrix(FI0,2,2,dimnames=list(names(param0),names(param0)))

    FisherInfo.fct <- function(param){
                   scale <- main(param)
                   PosDefSymmMatrix(FI0/scale^2)}

    if(is.function(trafo))  
        Tr0 <- trafo
    else if(is.matrix(trafo)){
            Tr0 <- trafo
            if(is.null(dimnames(trafo)))
               dimnames(trafo) <- list(names(param0),names(param0)) 

         }else{
            Tr0 <- matrix(trafo, dimnames = list(names(param0),names(param0) )) 
         }                    

    f.call <- substitute(L2ScaleUnknownLocationFamily(loc = l,
                 scale = s,
                 name = N,
                 centraldistribution = D0,
                 locscalename = lN,
                 modParam = mP,
                 LogDeriv = lD,
                 L2derivDistr.0 = L2D0,
                 FisherInfo.0 = F.0,
                 distrSymm = DSymm,
                 L2derivSymm = L2Symm,
                 L2derivDistrSymm = L2DSymm,
                 trafo = Tr,
                 .returnClsName = rtn),
             list(s = scale,
                 l = loc,
                 N = name,
                 D0 = centraldistribution,
                 lN = locscalename,
                 mP = modParam,
                 lD = LogDeriv,
                 L2D0 = as(L2derivDistr, "list"),
                 F.0 = FI0,
                 DSymm = distrSymm,
                 L2Symm = L2derivSymm,
                 L2DSymm = L2derivDistrSymm,
                 Tr = Tr0,
                 rtn = .returnClsName))
    if(is.null(.returnClsName))
       .returnClsName <- "L2LocationScaleFamily"
    L2Fam <- new(.returnClsName)
    L2Fam@name <- name
    L2Fam@locscalename <- locscalename
    L2Fam@distribution <- distribution
    L2Fam@distrSymm <- distrSymm
    L2Fam@param <- param
    L2Fam@startPar <- startPar
    L2Fam@makeOKPar <- makeOKPar
    L2Fam@modifyParam <- modParam
    L2Fam@fam.call <- f.call
    L2Fam@props <- props
    L2Fam@LogDeriv <- LogDeriv
    L2Fam@L2deriv.fct <- L2deriv.fct
    L2Fam@L2deriv <- L2deriv
    L2Fam@L2derivSymm <- L2derivSymm
    L2Fam@L2derivDistr <- L2derivDistr
    L2Fam@L2derivDistrSymm <- L2derivDistrSymm
    L2Fam@FisherInfo.fct <- FisherInfo.fct
    L2Fam@FisherInfo <- FisherInfo.fct(param)

    return(L2Fam)
}
