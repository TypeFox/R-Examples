##Constructor and Initializer of class 'yuima.data'

# we convert objects to "zoo" internally

 
setMethod("initialize", "yuima.data",
           function(.Object, original.data, delta=NULL){
             .Object@original.data <- original.data
             if(is.list(original.data) && is.zoo(original.data[[1]])) {
               .Object@zoo.data <- original.data
             } else {
               .Object@zoo.data <- as.list(as.zoo(.Object@original.data))
			 }
             if(!is.null(delta)){
                 delta <- rep(delta, length(.Object@zoo.data))
                 for(i in 1:length(.Object@zoo.data)){
                    n <- length(.Object@zoo.data[[i]])
                    t <- seq(0, delta[i]*n, length=n)
                    index(.Object@zoo.data[[i]]) <- t
                 }
             }
             return(.Object)
           })

# utils
onezoo <- function(ydata) {
  dat <- get.zoo.data(ydata)
  dats <- dat[[1]]
  if(length(dat)>1) {
    for(i in 2:(length(dat))) {
      dats <- merge(dats,dat[[i]])
    }
  }
  return(dats)
}

# accessors
setData <-
  function(original.data, delta=NULL){
    return(new("yuima.data", original.data=original.data, delta=delta ))
  }


setGeneric("get.zoo.data",
           function(x)
           standardGeneric("get.zoo.data")
           )
		   
setMethod("get.zoo.data", signature(x="yuima.data"),
          function(x){
            return(x@zoo.data)
          })

# following funcs are basic generic funcs

setGeneric("plot",
           function(x,y,...)
           standardGeneric("plot")
           )
		   

setMethod("plot",signature(x="yuima.data"),
          function(x,y,main="",xlab="index",ylab=names(x@zoo.data),...){
            plot(onezoo(x),main=main,xlab=xlab,ylab=ylab,...)
          }
          )

#setGeneric("time",
#           function(x,...)
#           standardGeneric("time")
#           )
		   
#setMethod("time", signature(x="yuima.data"),
#          function(x,...){
#            return(time(x@zoo.data))
#          }) 


#setGeneric("end",
#           def = function(x,...) standardGeneric("end")
#           )
		   
#setMethod("end", signature(x="yuima.data"),
#          function(x,...){
#            return(end(x@zoo.data))
#          }) 

#setGeneric("start",
#           function(x,...)
#           standardGeneric("start")
#           )
		   
#setMethod("start", signature(x="yuima.data"),
#          function(x,...){
#            return(start(x@zoo.data))
#          }) 

# length is primitive, so no standardGeneric should be defined


setMethod("length", signature(x= "yuima.data"),
          function(x){
		  #  if(is.null(dim(x@zoo.data)))
	      #    return(length(x@zoo.data))
	      #  else 
          #    return(dim(x@zoo.data)[1])
		    result <- numeric()
		    for(i in 1:(length(x@zoo.data))) result <- c(result,length(x@zoo.data[[i]]))
		    return(result)
          }
          ) 
		  
setMethod("dim", signature(x = "yuima.data"),
          function(x){
            return(length(x@zoo.data))
          }
          ) 


# same methods for 'yuima'. Depend on same methods for 'data'
setMethod("get.zoo.data", "yuima",
          function(x){
            return(get.zoo.data(x@data))
          })
setMethod("length", "yuima",
          function(x){
            return(length(x@data))
          }) 
setMethod("dim", "yuima",
          function(x){
           return(dim(x@data))
          }) 


setMethod("plot","yuima",
          function(x,y,xlab=x@model@time.variable,ylab=x@model@solve.variable,...){
		    if(length(x@model@time.variable)==0) {
              plot(x@data,...)
			} else {
              plot(x@data,xlab=xlab,ylab=ylab,...)
			}
          })


##:: yuima.data obj cbind ( implementation 08/18 )
setGeneric("cbind.yuima",
           function(x, ...)
           standardGeneric("cbind.yuima")
           )

setMethod("cbind.yuima", signature(x="yuima"),
          function(x, ...){
            ##:: init
            y.list <- list(x, ...)
            y.num <- length(y.list)

            ##:: bind yuima.data in yuima
            yd.tmp <- y.list[[1]]@data
            for(idx in 2:y.num){
              ##:: error check
              if( class(y.list[[idx]])!="yuima"){
                stop("arg ", idx, " is not yuima-class")
              }
              ##:: bind
              yd.tmp <- cbind.yuima(yd.tmp, y.list[[idx]]@data)
            }

            ##:: substitute yuima.data
            x@data <- yd.tmp

            ##:: return result
            return(x)
          }
          )

setMethod("cbind.yuima", signature(x="yuima.data"),
          function(x, ...){
            ##:: init
            yd.list <- list(x, ...)
            yd.num <- length(yd.list)
            
            ##:: bind yuima.data (original.data)
            od.tmp <- yd.list[[1]]@original.data
            for(idx in 2:yd.num){
              ##:: error check
              if( class(yd.list[[idx]])!="yuima.data" ){
                stop("arg ", idx, " is not yuima.data-class.")
              }
              ##:: bind
              od.tmp <- cbind(od.tmp, yd.list[[idx]]@original.data)
            }
            ##:: return result
            return(new("yuima.data", original.data=od.tmp))
          }
          )

##:: END ( yuima.data obj cbind )
