#
# vim:set ff=unix expandtab ts=2 sw=2:
TimeMap.new=function#basic constructor of the class TimeMap.
### A TimeMap is nothing more than an usual R-function of one argument augmented by the lower and upper boundary of the interval where it is defined.
(t_start, ##<<A number marking the begin of the time domain where the function is valid
 t_end,   ##<<A number the end of the time domain where the function is valid
 f        ##<<The time dependent function definition (a function in R's sense)
 ){
   warning("This function is going deprecated for 2 reasons:\n
      1.) There are new more specialized classes to replace it.( BoundedInFlux,BoundedLinDecompOp)\n
      2.) Constructors for SoilR classes have been renamed consistently to the name of the class (NameOfClass( ) instead of NameOfClass.new() )
   ")

   obj=new(Class="TimeMap",t_start,t_end,f) 
return(obj)
### An object of class TimeMap that can be used to describe models.
}
##########################################################################

### defines a (time dependent) mapping including the function definition and the ### domain where the function is well define.  This can be used to avoid interpolations out of range when mixing different time dependent data sets
setClass(
   Class="TimeMap",
   #contains="UnlimitedTimeMap",
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
##########################################################################
setMethod(
    f="initialize",
    signature="TimeMap",
    definition=function #initialize called by (new)
    ### This is the standard constructor for objects of this class
    ### It is called by statements of the form 
    ### \code{new("TimeMap",start,end,f,lag)}
    (.Object,
    starttime=numeric(),
    endtime=numeric(),
    map=function(t){t},
    lag=0
    ){
   warning("This function is going to be deprecated :\n
      There are new more specialized classes to replace it.( BoundedInFlux,BoundedLinDecompOp)\n
      2.) There are Constructors for SoilR classes, named consistently to the name of the class (NameOfClass( ). For stable code rather use those instead of new(\"NameOfClass\",...) )
   ")
    #cat("- initializer of Class Time Map at work - \n")
    .Object@starttime=starttime
    .Object@endtime=endtime
    .Object@map=map
    .Object@lag=lag
    return(.Object)
    }
)
##########################################################################
setMethod(
    f="as.character",
    signature="TimeMap",
    definition=function #convert TimeMap Objects to something printable.
    ### This method is needed to print a TimeMap object.
    (x, ##<<An Object of class time map
     ...
     ){
        return(
            paste( class(x),
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
##########################################################################
setMethod(
    f="getTimeRange",
    signature="TimeMap",
    definition=function # ask for the boundaries of the underlying time interval
    ### The method returns the time range of the given object 
    ### It is probably mostly used internally to make sure that 
    ### time dependent functions retrieved from data are not
    ### used outside the interval where they are valid. 
    
    (object ##<< An object of class TimeMap or one that inherits from TimeMap
    ){
        return( c("t_min"=object@starttime,"t_max"=object@endtime))
        ### a vector of length two \code{ c(t_min,t_max) }
        ### containing start and end time of the time interval 
        ### for which the TimeMap object has been defined.
    }
)
###########################################################################
setMethod(
    f="getFunctionDefinition",
    signature="TimeMap",
    definition=function(object){
    ### extract the function definition (the R-function) 
        return(object@map)
    }
)

##########################################################################
TimeMap.from.Dataframe=function
### This function is another constructor of the class TimeMap.
(dframe, ##<<A data frame containing exactly two columns:
## the first one is interpreted as time
lag=0, ##<< a scalar describing the time lag. Positive Values shift the argument of the interpolation function forward in time. (retard its effect)
interpolation=splinefun ##<<A function that  returns a function  the default is splinefun. Other possible values are the linear interpolation approxfun or any self made function with the same interface.
 ){
   warning("This function will be deprecated for 2 reaseons: \n 1.) (frequently generic) constructors are now called like the classes they produce objects of.  \n 2.) There are specialized replacements for class TimeMap (BoundFc, BoundLinDecompOp). ")
   t=dframe[,1]  
   y=dframe[,2]  
   o=order(t)
   tyo=cbind(t[o],y[o])
   to=tyo[,1]+lag# account for the lag time
   yo=tyo[,2]
   t_start=min(to)
   t_start=min(t)
   t_end=max(t)
   interpol=interpolation(to,yo)
   obj=new(Class="TimeMap",t_start,t_end,interpol) 
return(obj)
### An object of class TimeMap that contains the interpolation function and the limits of the time range where the function is valid. Note that the limits change according to the time lag
### this serves as a saveguard for Model which thus can check that all involved functions of time are actually defined for the times of interest  
}
#########################################################
setMethod(
      f="BoundLinDecompOp",
      signature=c(map="TimeMap",starttime="missing",endtime="missing",lag="missing"),
      definition=function # create a BoundLinDecompOp from a TimeMap
      ### The method is used internally to convert TimeMap objects to BoundLinDecompOp where the use of TimeMap is now deprecated.
      (map){
      starttime=map@starttime
      endtime=map@endtime
      map=map@map
      return(BoundLinDecompOp(map,starttime,endtime))
     }
     )
#########################################################
setMethod(
      f="BoundInFlux",
      signature=c("TimeMap","missing","missing","missing","missing"),
      definition=function # convert to BoundInFlux
      ### The method is used internally to convert TimeMap objects to BoundInFlux objects, since the use of TimeMap objects is now deprecated.
      (map){
      starttime=map@starttime
      endtime=map@endtime
      map=map@map
      return(BoundInFlux(map,starttime,endtime))
     }
     )
