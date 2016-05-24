makeAbscontDistribution <- function(object, gaps = NULL,
                       param = NULL, img = NULL,
                   withgaps = getdistrOption("withgaps"),
                   ngrid = getdistrOption("DefaultNrGridPoints"),
                   ep = getdistrOption("TruncQuantile")){
   if(!is(object,"UnivariateDistribution"))
        stop("Argument 'object' must be of class 'UnivariateDistribution'")
   if(missing(img)) img0 <- img(object)
   if(is.null(img)) img0 <- img(object)
   pfun <- p(object)
   low0 <- q(object)(0)*1.001
   up0 <- q(object)(1)*1.001
   low1 <- getLow(object,ep)*1.001
   up1 <- getUp(object,ep)*1.001
   wS <- object@.withSim
   wA <- object@.withArith
   lE <- .lowerExact(object)
   loE <- .logExact(object)
   AbscontDistribution(p = pfun, gaps = gaps, param = param, img = img0,
                   .withSim = wS, .withArith = wA, .lowerExact = lE,
                   .logExact = loE, withgaps = withgaps, low1 = low1, up1 = up1,
                   low = low0, up = up0, ngrid = ngrid,
                   ep = ep)
   }

   
