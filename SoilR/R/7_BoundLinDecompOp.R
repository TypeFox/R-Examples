#
# vim:set ff=unix expandtab ts=2 sw=2:

setClass(# a decomposition operator described by a matrix valued function of time
    Class="BoundLinDecompOp",
    contains="DecompOp",   
    slots=list(
    map="function"
    ,
    lag="numeric"
    ,
    starttime="numeric"
    ,
    endtime="numeric"
    ) 
    
   )
#---------------------------------------------------------------------
setMethod(
    f="initialize",
    ### 
    signature="BoundLinDecompOp",
    definition=function #initialize called by (new)
    ### This is the internal constructor for objects of this class
    ### It is called by statements of the form new("..) 
    ### which in turn is called by all other constructors 
    ### which may have more convienient interfaces
    (.Object,
    starttime=numeric(),
    endtime=numeric(),
    map=function(t){t},
    lag=0
    ){
    #cat("- initializer of Class Time Map at work - \n")
    .Object@starttime=starttime
    .Object@endtime=endtime
    .Object@map=map
    .Object@lag=lag
    return(.Object)
    }
)
#---------------------------------------------------------------------
setMethod(
      f="BoundLinDecompOp",
      signature=c(map="function",starttime="numeric",endtime="numeric",lag="numeric"),
      definition=function # a constructor 
      ### This method creates a BoundLinDecompOp from a timedependent function and its domain 
      (
        map, ##<<a function
        starttime,##<< the begin of the time domain
        endtime,##<< the end of the time domain
        lag ##<< lag time
      ){
      return(new("BoundLinDecompOp",starttime=starttime,endtime=endtime,map=map,lag=lag))
    }
)
#---------------------------------------------------------------------
setMethod(
      f="BoundLinDecompOp",
      signature=c(
        map="function",
        starttime="numeric",
        endtime="numeric",
        lag="missing"
        ),
      definition=function # a constructor 
      ### This method creates a BoundLinDecompOp from a timedependent function and its domain 
      (
        map,
        starttime,
        endtime
      ){
      return(BoundLinDecompOp(starttime=starttime,endtime=endtime,map=map,lag=0))
    }
)
#---------------------------------------------------------------------
#There will be a new Class LinDecompOp (without Bound ) for this 
#setMethod(
#      f="BoundLinDecompOp",
#      ### 
#      signature=c(map="function",starttime="missing",endtime="missing",lag="numeric"),
#      definition=function # a constructor for a single function without limits  
#      ### This method creates a BoundLinDecompOp from a timedependent function and a lag only
#      (map,lag=0){
#      return(BoundLinDecompOp(starttime=-Inf,endtime=Inf,map=map))
#    }
#)
#---------------------------------------------------------------------
setMethod(
    f="getTimeRange",
    signature="BoundLinDecompOp",
    definition=function # ask for the boundaries of the underlying time interval
    ### The method returns the time range of the given object 
    ### It is ( probably mostly ) used internally to make sure that 
    ### time dependent functions retrieved from data are not
    ### used outside the interval where they are valid. 
    
    (object 
    ){
        return( c("t_min"=object@starttime,"t_max"=object@endtime))
        ### a vector of length two \code{ c(t_min,t_max) }
        ### containing start and end time of the time interval 
        ### for which the object has been defined.
    }
)
#---------------------------------------------------------------------
setMethod(
    f="getFunctionDefinition",
    signature="BoundLinDecompOp",
    definition=function(object){
    ### extract the function definition (the R-function) 
        return(object@map)
    }
)

