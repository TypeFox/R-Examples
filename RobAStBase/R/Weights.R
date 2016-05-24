setMethod("name", "RobWeight", function(object) object@name)
setReplaceMethod("name", "RobWeight", 
    function(object, value){
        object@name <- value
        object
    })


setMethod("clip", "BoundedWeight", function(x1) x1@clip)
setReplaceMethod("clip", "BoundedWeight", 
    function(object, value){
        object@clip <- value
        object
    })


setMethod("stand", "BdStWeight", function(object) object@stand)
setReplaceMethod("stand", "BdStWeight", 
    function(object, value){
        object@stand <- value
        object
    })


setMethod("cent", "HampelWeight", function(object) object@cent)
setReplaceMethod("cent", "HampelWeight", 
    function(object, value){
        object@cent <- value
        object
    })

setMethod("weight", "RobWeight", function(object) object@weight)
setReplaceMethod("weight", "RobWeight", 
    function(object, value){
        object@weight <- value
        object
    })

setMethod("getweight",
          signature(Weight = "HampelWeight", neighbor = "ContNeighborhood",
                    biastype = "BiasType"),# normtype = "NormType"),
          function(Weight, neighbor, biastype, normW)
               {A <- stand(Weight)
                b <- clip(Weight)
                z <- cent(Weight)
                function(x){
                   y <- A%*%(x-z)
                   norm0 <- fct(normW)(y) 
                   ind2 <- (norm0 < b/2)
                   norm1 <- ind2*b/2 + (1-ind2)*norm0
                   ind1 <- (norm0 < b)
                   ind1 + (1-ind1)*b/norm1
                   }
                }
          )


setMethod("getweight",
          signature(Weight = "HampelWeight", neighbor = "ContNeighborhood",
                    biastype = "onesidedBias"),#  norm = "missing"),
          function(Weight, neighbor, biastype, ...)
               {A <- stand(Weight)
                b <- clip(Weight)
                z <- cent(Weight)
                function(x){
                   y <- as.numeric(as.matrix(A)%*%(x-z))*sign(biastype)
                   norm1 <- pmax(y,b/2)
                   pmin(1,b/norm1) 
                   }
                }
          )

setMethod("getweight",
          signature(Weight = "HampelWeight", neighbor = "ContNeighborhood",
                    biastype = "asymmetricBias"),# norm = "missing"),
          function(Weight, neighbor, biastype, ...)
               {A <- stand(Weight)
                b <- clip(Weight)
                b1 <- b/nu(biastype)[1]
                b2 <- b/nu(biastype)[2]
                z <- cent(Weight)
                function(x){
                   y <- as.numeric(as.matrix(A)%*%(x-z))
                   norm1 <- pmax(-y,b1/2)
                   norm2 <- pmax(y,b2/2)
                   pmin(1,b1/norm1,b2/norm2) 
                   }
                }
          )


setMethod("getweight",
          signature(Weight = "BdStWeight", neighbor = "TotalVarNeighborhood",
                    biastype = "BiasType"),#  norm = "missing"),
          function(Weight, neighbor, biastype, ...)
               {A <- stand(Weight)
                b <- clip(Weight)
                b1 <- -b[1]
                b2 <- b[2]
                function(x){
                   y <- as.numeric(as.matrix(A)%*%x)
                   norm1 <- pmax(-y,b1/2)
                   norm2 <- pmax(y,b2/2)
                   pmin(1,b1/norm1,b2/norm2) 
                   }
                }
          )

setMethod("minbiasweight",
          signature(Weight = "HampelWeight", neighbor = "ContNeighborhood",
                    biastype = "BiasType"),#  norm = "NormType"),
          function(Weight, neighbor, biastype, normW)
               {A <- stand(Weight)
                b <- clip(Weight)
                z <- cent(Weight)
                function(x){
                   y <- A%*%(x-z)
                   norm0 <- fct(normW)(y) 
                   ind <- 1-.eq(norm0)
                   ind*b/(norm0+1-ind)
                   }
                }
          )


setMethod("minbiasweight",
          signature(Weight = "HampelWeight", neighbor = "ContNeighborhood",
                    biastype = "asymmetricBias"),#  norm = "missing"),
          function(Weight, neighbor, biastype, ...)
               {A <- stand(Weight)
                b <- clip(Weight)
                b1 <- -b[1]
                b2 <- b[2]
                z <- cent(Weight)
                function(x){
                   y <- as.numeric(as.matrix(A)%*%(x-z))
                   indp <- (y>0)
                   ind0 <- .eq(y)
                   indm <- (y<0)
                   indm*b1/(y+ind0) + indp*b2/(y+ind0)
                   }
                }
          )

setMethod("minbiasweight",
          signature(Weight = "HampelWeight", neighbor = "ContNeighborhood",
                    biastype = "onesidedBias"),#  norm = "missing"),
          function(Weight, neighbor, biastype, ...)
               {A <- stand(Weight)
                b <- clip(Weight)
                z <- cent(Weight)
                function(x){
                   y <- as.numeric(as.matrix(A)%*%(x-z))
                   ind <- (y*sign(biastype) >0)
                   ind0 <- .eq(y)
                   ind*b/(y+ind0)+(1-ind)
                   }
                }
          )


setMethod("minbiasweight",
          signature(Weight = "BdStWeight", neighbor = "TotalVarNeighborhood",
                    biastype = "BiasType"),
          function(Weight, neighbor, biastype, ...)
               {A <- stand(Weight)
                b <- clip(Weight)
                b1 <- b[1]
                b2 <- b[2]
                function(x){
                   y <- as.numeric(as.matrix(A)%*%(x))
                   indp <- (y>0)
                   ind0 <- .eq(y)
                   indm <- (y<0)
                   indm*b1/(y+ind0) + indp*b2/(y+ind0)
                   }
                }
          )
