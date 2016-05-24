
################################
##
## Class: Laplace or Double Exponential distribution
##
################################


## wrapped access methods
setMethod("rate", "DExp", function(object) rate(param(object)))

## wrapped replace methods
setMethod("rate<-", "DExp", function(object, value) new("DExp", rate = value))

DExp <- function(rate = 1) {D <- new("DExp", rate = rate)
                            D@Symmetry <- SphericalSymmetry(0)
                            D}

setMethod("*", c("DExp","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            if(isTRUE(all.equal(e2,0))) 
               return(new("Dirac", location = 0, .withArith = TRUE))
            else  
               return(new("DExp", rate = rate(e1) / abs(e2), .withArith = TRUE))
          })

