#
# vim:set ff=unix expandtab ts=2 sw=2:
BoundInFlux.new=function#constructor of the class BoundInFlux.
### A BoundInFlux is nothing more than an usual R-function of one argument augmented by the lower and upper boundary of the interval where it is defined.
(t_start, ##<<A number marking the begin of the time domain where the function is valid
 t_end,   ##<<A number the end of the time domain where the function is valid
 f        ##<<The time dependent function definition (a function in R's sense)
 ){
   warning("The function is deprecated, constructors are now called as the class. To get rid of this warning please use BoundInFlux() in the future.")
   obj=BoundInFlux(f,t_start,t_end) 
return(obj)
### An object of class BoundInFlux that can be used to describe models.
}
#-----------------------------------------------------------------------------------

### defines a time dependent inputrate as function of time and including the domain where the function is well defined. This can be used to avoid interpolations out of range when mixing different time dependent data sets
setClass(
   Class="BoundInFlux",
   contains="InFlux",
   slots=list(
	starttime="numeric"
    ,
	endtime="numeric"
    ,
    map="function"
    ,
    lag="numeric"
   )
)
#------------------------------ constructors------------------------------------------
setMethod(
    f="initialize",
    signature="BoundInFlux",
    definition=function # internal constructor
    ### This mehtod is intended for internal use only, it may change with the internal representation of the class. In user code please use the generic constructor \code{\link{BoundInFlux}} instead.
    (.Object,starttime=numeric(),endtime=numeric(),map=function(t){t},lag=0){
    #cat("-initializer at work-\n")
    .Object@starttime=starttime
    .Object@endtime=endtime
    .Object@map=map
    .Object@lag=lag
    return(.Object)
    }
)
setMethod(
  f="BoundInFlux",
  signature=c(
    map="function",
    starttime="numeric",
    endtime="numeric",
    lag="numeric",
    interpolation="missing"
  ),
  definition=function # constructor 
  ### the method constructs an object from its basic ingredients
  (
    map,
    starttime,
    endtime,
    lag
    ){
    new("BoundInFlux",map=map,starttime=starttime,endtime=endtime,lag=lag)
  }
)
setMethod(
  f="BoundInFlux",
  signature=c(map="function",starttime="numeric",endtime="numeric",lag="missing",interpolation="missing"),
  definition=function # constructor 
  ### the method constructs an object from its basic ingredients
  (map,starttime,endtime){
    BoundInFlux(map=map,starttime=starttime,endtime=endtime,lag=0)
  }
)
setMethod(
  f="BoundInFlux",
  signature=c(map="data.frame",starttime="missing",endtime="missing",lag="numeric",interpolation="function"),
  definition=function #constructor
  ### This function is another constructor of the class BoundInFlux.
  (
    map ,##<<A data frame; the first column is interpreted as time
    lag, ##<< lag time
    interpolation ##<< function used for interpolation
  ){
     t=map[,1]  
     y=map[,2]  
     o=order(t)
     tyo=cbind(t[o],y[o])
     #to=tyo[,1]+lag# account for the lag time
     to=tyo[,1] # since lag is also part of the result
     yo=tyo[,2]
     t_start=min(to)
     t_start=min(t)
     t_end=max(t)
     interpol=interpolation(to,yo)
     obj <- BoundInFlux(map=interpol,starttime=t_start,endtime=t_end,lag=lag) 
  return(obj)
  ### An object of class BoundInFlux that contains the interpolation function and the limits of the time range where the function is valid. Note that the limits change according to the time lag
  }
)

setMethod(
  f="BoundInFlux",
  signature=c(map="data.frame",starttime="missing",endtime="missing",lag="missing",interpolation="missing"),
  definition=function #constructor
  ### This function is another constructor of the class BoundInFlux.
  (map##<<A data frame; the first column is interpreted as time
   ){
     obj=BoundInFlux(map=map,lag=0,interpolation=splinefun) 
  return(obj)
  ### An object of class BoundInFlux that contains the interpolation function and the limits of the time range where the function is valid. Note that the limits change according to the time lag
  }
)
#-------------------------------------------other methods---------------------------------------------------
setMethod(
    f="as.character",
    signature="BoundInFlux",
    definition=function # convert BoundInFlux Objects to something printable.
    ### returns a string describing the object
    (x, ##<< An object 
     ...
     ){
        return(
            paste( class(x),
                  "(\n map=",
                  x@map,
                  "(\n starttime=",
                  x@starttime,
                  "\n endtime=",
                  x@endtime,
                  ")",
                  sep=""
            )
        )
    }
)    
setMethod(
    f="getTimeRange",
    signature="BoundInFlux",
    definition=function # time domain of the function
    ### The method returns a vector containing the start and end time where the intepolation is valid.
    ( object){
        return(
               c("t_min"=object@starttime,"t_max"=object@endtime))
    }
)
setMethod(
    f="getFunctionDefinition",
    signature="BoundInFlux",
    definition=function(object){
    ### extract the function definition (the R-function) from the BoundInFlux 
        return(object@map)
    }
)
