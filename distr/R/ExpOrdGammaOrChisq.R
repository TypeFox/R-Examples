######## ExpOrGammaOrChisq - methods

setMethod("*", c("ExpOrGammaOrChisq","numeric"),
          function(e1, e2){
            if(is(e1,"Gammad")){
               if(isTRUE(all.equal(e2,0)))
                  return(new("Dirac", location = 0, .withArith = TRUE))
               if(e2 > 0)
                  return(new("Gammad", shape = shape(e1),
                              scale = scale(e1) * e2, .withArith = TRUE))
               return(getMethod("*",c("AbscontDistribution","numeric"))(
                                      e1 = Gammad(shape = shape(e1),
                                                  scale = scale(e1) * (-e2)), 
                                      e2 = -1)) 
             }                         
             else return(getMethod("*",
                            c("AbscontDistribution","numeric"))(e1, e2) )
          })

setMethod("+", c("ExpOrGammaOrChisq","ExpOrGammaOrChisq"),
          function(e1,e2){
            if(is(e1,"Gammad")&&is(e2,"Gammad"))
               {e10 <- as(e1,"Gammad")
                e20 <- as(e2,"Gammad")
                newshape <- shape(e10) + shape(e20)
                if(is.logical(all.equal(scale(e10),scale(e20))))
                   return(new("Gammad", shape = newshape,
                               scale = scale(e10), .withArith = TRUE))
               }
            return(as(e1, "AbscontDistribution") + e2)
          })
