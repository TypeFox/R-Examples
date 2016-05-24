getRadius <- function(IC, risk, neighbor, L2Fam){
   if(!is(IC, "HampIC")) stop("'IC' must be of class 'HampIC'.")
   if(!is(risk,"asGRisk")) stop("'risk' must be of class 'asGRisk'.")
   if(!is(neighbor,"UncondNeighborhood"))
         stop("'neighbor' must be of class 'UncondNeighborhood'.")
   if (missing(L2Fam)) L2Fam <- force(eval(IC@CallL2Fam))
   if(!is(L2Fam,"L2ParamFamily")) stop("'L2Fam' must be of class 'L2ParamFamily'.")
   L2deriv.0 <- L2Fam@L2deriv
   if(numberOfMaps(L2deriv.0)==1){ ## L2derivDim <- L2Fam@L2deriv
      z <- cent(IC)/as.vector(stand(IC))
      c0 <- clip(IC)/abs(as.vector(stand(IC)))
      symm <- FALSE
      if(is(L2Fam@L2derivDistrSymm[[1]], "SphericalSymmetry"))
         symm <- L2Fam@L2derivDistrSymm[[1]]@SymmCenter == 0
      r <- getInfRad(clip = c0, L2deriv = L2Fam@L2derivDistr[[1]],
                     risk = risk, neighbor = neighbor, biastype = biastype(risk),
                     cent = z, symm = symm, trafo = trafo(L2Fam@param))
   }else{
      if(!is(L2Fam@distribution, "UnivariateDistribution"))
         stop("not yet implemented")
      if((length(L2Fam@L2deriv) == 1) & is(L2Fam@L2deriv[[1]], "RealRandVariable")){
                    L2deriv <- L2Fam@L2deriv[[1]]
      }else{
                    L2deriv <- diag(dimension(L2Fam@L2deriv)) %*% L2Fam@L2deriv
      }
      z <- solve(stand(IC),cent(IC))
      r <- getInfRad(clip = clip(IC), L2deriv = L2deriv,
                     risk = risk, neighbor = neighbor, biastype = biastype(risk),
                     Distr = L2Fam@distribution, stand = stand(IC),
                     cent = z , trafo = trafo(L2Fam@param))
   }
   return(r)
}
