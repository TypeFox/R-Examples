## thanks to Martin Morgan

setMethod("[",
          signature=signature(x = "haploList",i = "ANY",j = "missing"),
          function(x,i,j,...,drop = TRUE){
            ## update and return x
            initialize(x,x@.Data[i])
          })

#from Gentleman 2008, p 109
setMethod("c",
          signature=signature(x = "haploList"),
          function(x, ..., recursive = FALSE){
            if(nargs() == 2)
              cpair(x,...)
            else if (nargs() > 2)
              cpair(x,c(...))
            else
              x
          })

setGeneric("cpair",function(x,y)standardGeneric("cpair"))

setMethod("cpair",signature = signature(x = "haploList",y = "haploList"),
          function(x,y){
            if(!identical(attributes(x),attributes(y)))
              stop("x and y can not be paired")
            initialize(x,x@.Data <- c(x@.Data,y@.Data))
          })

setMethod("cpair",signature = signature(x = "haploList",y = "haplotype"),
          function(x,y){
            initialize(x,x@.Data <- c(x@.Data,y))
          })

setMethod("cpair",signature = signature(x = "haploList",y = "list"),
          function(x,y){
            y <- haploList(list = y)
            initialize(x,x@.Data <- c(x@.Data,y))
          })
          

                     
          
                    
           
          


