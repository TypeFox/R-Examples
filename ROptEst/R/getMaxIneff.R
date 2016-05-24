getMaxIneff <- function(IC, neighbor, biastype = symmetricBias(), 
                        normtype = NormType(), z.start = NULL, 
                        A.start = NULL, maxiter = 50, 
                        tol = .Machine$double.eps^0.4,
                        warn = TRUE, verbose = NULL){
            if(!is(IC,"IC")) 
               stop("Argument IC must be of class 'IC'.")

            ow <- options("warn")
            on.exit(options(ow))

            sb <- .getSB(IC,neighbor)
            si <- sb$s^2
            bi <- sb$b^2
            
            risk <- asBias(biastype=biastype, normtype=normtype)

            L2Fam <- eval(IC@CallL2Fam)
            trafo <- trafo(L2Fam@param)
            symm <- L2Fam@L2derivDistrSymm[[1]]
            Finfo <- L2Fam@FisherInfo
            L2derivDim <- numberOfMaps(L2Fam@L2deriv)
            
            FI0 <- trafo%*%solve(Finfo)%*%t(trafo)
            std <- if(is(normtype,"QFNorm")) 
                      QuadForm(normtype) else diag(nrow(trafo))
            s0 <- sum(diag(std%*%FI0))
            
            if(L2derivDim==1){
              L2deriv <- L2Fam@L2derivDistr[[1]]
              b0 <- getInfRobIC(L2deriv = L2deriv, risk = risk, 
                                neighbor = neighbor, symm = symm, 
                                trafo = trafo, maxiter = maxiter,  
                                tol = tol, warn = warn, Finfo = Finfo, 
                                verbose = verbose)$risk$asBias$value^2
            }else{ 
              if(is(L2Fam@distribution, "UnivariateDistribution")){
                 if((length(L2Fam@L2deriv) == 1) & is(L2Fam@L2deriv[[1]], "RealRandVariable")){
                    L2deriv <- L2Fam@L2deriv[[1]]
                    L2derivSymm <- L2Fam@L2derivSymm
                    L2derivDistrSymm <- L2Fam@L2derivDistrSymm
                }else{
                    L2deriv <- diag(dimension(L2deriv)) %*% L2deriv
                    L2deriv <- RealRandVariable(Map = L2deriv@Map, 
                                                Domain = L2deriv@Domain)
                    nrvalues <- numberOfMaps(L2deriv)
                    if(numberOfMaps(L2deriv) != nrvalues){
                        L1 <- vector("list", nrvalues)
                        L2 <- vector("list", nrvalues)
                        for(i in 1:nrvalues){
                            L1[[i]] <- NonSymmetric()
                            L2[[i]] <- NoSymmetry()
                        }
                        L2derivSymm <- new("FunSymmList", L1)
                        L2derivDistrSymm <- new("DistrSymmList", L2)
                    }
                }
                if(!warn) options(warn = -1)

                b0 <- getInfRobIC(L2deriv = L2deriv, neighbor = neighbor,
                            risk = risk,  Distr = L2Fam@distribution, 
                            DistrSymm = L2Fam@distrSymm, L2derivSymm = L2derivSymm,
                            L2derivDistrSymm = L2derivDistrSymm, 
                            trafo = trafo, z.start = z.start, A.start = A.start, 
                            maxiter = maxiter, tol = tol, warn = warn, 
                            Finfo = Finfo,
                            verbose = verbose)$risk$asBias$value^2
              }else{
                stop("not yet implemented")
              }

            } 
         return(max(si/s0,bi/b0))
}                               
             

  