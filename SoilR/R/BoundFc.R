#
# vim:set ff=unix expandtab ts=2 sw=2:
correctnessOfBoundFc=function#check for unreasonable parameters or unsupported formats
###  The atmospheric C14 data can be represented in more than one format 
###  The function checks if the user required format is supported at the moment
(object ##<< the object to be tested
)
{
   res=TRUE
   supported_formats=supported14CFractionFormats()
   f=object@format
   print(paste("format=",f))
   if (!any(grepl(f,supported_formats))){
      err_str=cat("The required format:",f," describing the atmospheric c_14 fraction is not supported.\n 
   	     The supported formats are: ",supported_formats,". \n",sep="")
      stop(simpleError(err_str))
      return(res)
   }
}

#---------------------------------------------------------------------------------------------------------
setClass(# Objects containing the atmospheric 14C fraction and the format it is provided in. 
    ### Objects of this class contain a time dependent function describing the Atmospheric \eqn{^{14}C}{14C} fraction and a format description, that allows to use the numeric valuest to be interpreted correctly in subsequent computations.
    Class="BoundFc",
    slots=list(
      starttime="numeric" ,
      endtime="numeric" ,
      map="function" ,
      lag="numeric" ,
      format="character" 
   )
    #,validity=correctnessOfBoundFc #set the validating function
)
#------------------------------------------ ------------------------------------------------------------------
#----------------------------- Constructors ------------------------------------------------------------------
#------------------------------------------ ------------------------------------------------------------------
setMethod(
    f="initialize",
    signature="BoundFc",
    definition=function #internal constructor (new)
    ### The function is probably called internally only but used by all other constructors
    ### It calls a sanity check on its arguments and initialized the object 
    (.Object,map=function(t){t},starttime=numeric(),endtime=numeric(),lag=0,format="Delta14C",interpolation=splinefun){
    #cat("-initializer at work-\n")
    .Object@map=map
    .Object@starttime=starttime
    .Object@endtime=endtime
    .Object@lag=lag
    .Object@format=format
    correctnessOfBoundFc(.Object)
    return(.Object)
    }
)
#----------------------------- Constructors with function-------------------------------------------------
#---------------------------------------------------------------------------------------------------------
setMethod(
  f="BoundFc",
  signature=c(map="function",starttime="numeric",endtime="numeric",lag="numeric",format="character",interpolation="missing"),
  definition=function # constructor
  ### the method constructs an object from a function a timerange where it is valid and a format  
(
  map,        ##<< a function of one argument (time) 
  starttime,  ##<< the point in time from which map is a valid representation 
  endtime,    ##<< the point in time until which map is a valid representation
  lag,        ##<< a scalar describing the time lag. Positive Values shift the argument of the interpolation function forward in time. (retard its effect)
  format     ##<< a string that specifies the format used to represent the atmospheric fraction. Possible values are "Delta14C" which is the default or "afn" the Absolute Fraction Normal representation 
){
  obj=new(
    Class="BoundFc",
    map=map,
    starttime=starttime,
    endtime=endtime,
    lag=lag,
    format=format,
  ) 
return(obj)
}
)
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
setMethod(
  f="BoundFc",
  signature=c(map="function",starttime="numeric",endtime="numeric",lag="missing",format="character",interpolation="missing"),
  definition=function # constructor
  ### wrapper for \code{\link{BoundFc_method__function_numeric_numeric_numeric_character_missing}} with the assumption lag=0
  (
    map, 
    starttime,
    endtime,
    format 
  ){
    return(BoundFc(map,starttime,endtime,lag=0,format))
    ### An object  that contains the interpolation function and the limits of the time range where the function is valid. Note that the limits change according to the time lag
  }
)
#----------------------------- Constructors with dataframe -------------------------------------------------
#---------------------------------------------------------------------------------------------------------
setMethod(
  f="BoundFc",
  signature=c(map="data.frame",starttime="missing",endtime="missing",lag="numeric",format="character",interpolation="function"),
  definition=function # constructor 
  ### the method constructs an object from a dataframe a timelag format using the given interpolating function
(
  map,            ##<<A data frame containing exactly two columns:
                  ## the first one is interpreted as time
                  ## the second one is interpreted as atmospheric C14 fraction in the format mentioned
  lag=0,          ##<< a scalar describing the time lag. Positive Values shift the argument of the interpolation function forward in time. (retard its effect)
  format,         ##<< a string that specifies the format used to represent the atmospheric fraction. Possible values are "Delta14C" which is the default or "afn" the Absolute Fraction Normal representation 
  interpolation   ##<<A function that  returns a function  the default is splinefun. Other possible values are the linear interpolation approxfun or any self made function with the same interface.
){
   t=map[,1]  
   y=map[,2]  
   o=order(t)
   tyo=cbind(t[o],y[o])
   to=tyo[,1]+lag# account for the lag time
   yo=tyo[,2]
   t_start=min(to)
   t_start=min(t)
   t_end=max(t)
   interpol=interpolation(to,yo)
   #obj=new(Class="BoundFc",interpol,t_start,t_end,lag=lag,format=format) 
   obj=BoundFc(map=interpol,starttime=t_start,endtime=t_end,lag=lag,format=format) 
return(obj)
### An object  that contains the interpolation function and the limits of the time range where the function is valid. Note that the limits change according to the time lag
}
)
#---------------------------------------------------------------------------------------------------------
setMethod(
  f="BoundFc",
  signature=c(map="data.frame",starttime="missing",endtime="missing",lag="numeric",format="character",interpolation="missing"),
  definition=function # constructor
  ### wrapper for \code{\link{BoundFc_method__data.frame_missing_missing_numeric_character_missing}} with the assumption interpolation =splinefun
(map, 
lag,
format 
){
   obj=BoundFc(map,lag=lag,format=format,interpolation=splinefun) 
return(obj)
}
)
#---------------------------------------------------------------------------------------------------------
setMethod(
  f="BoundFc",
  signature=c(map="data.frame",starttime="missing",endtime="missing",lag="missing",format="character",interpolation="function"),
  definition=function # constructor 
  ### wrapper for \code{\link{BoundFc_method__data.frame_missing_missing_numeric_character_missing}} with the assumption lag=0
(map, 
format, 
interpolation
){
   obj=BoundFc(map,lag=0,format=format,interpolation=interpolation) 
return(obj)
}
)
#---------------------------------------------------------------------------------------------------------
setMethod(
  f="BoundFc",
  signature=c(map="data.frame",starttime="missing",endtime="missing",lag="missing",format="character",interpolation="missing"),
  definition=function # constructor 
  ### wrapper for \code{\link{BoundFc_method__data.frame_missing_missing_numeric_character_missing}} with the assumption lag=0 interpolation = splinefun
(map, 
format, 
interpolation
){
   obj=BoundFc(map,lag=0,format=format,interpolation=splinefun) 
return(obj)
}
)
#----------------------------- end Constructors ------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------
setMethod(
    f="getFormat",
    signature="BoundFc",
    definition=function# extract the format string
    ### the function just yields the format as a string
	  (object ##<< object  containing imformation about the format that could be Delta14C or AFM (Absolute Fraction Modern) for instance
		){
        return(object@format)
    }
)
#---------------------------------------------------------------------------------------------------------
setMethod(
   f= "Delta14C",
   signature("BoundFc"),
   definition=function# convert to Absolute Fraction Normal values  
   ### convert object containing values in any supported format to the appropriate Absolute Fraction Modern values.
   (
    F ##<< object of containing the values in any format
	 ){
      f=F@format
            targetFormat="Delta14C"
            if (f==targetFormat){
	       # do nothing
	       return(F)
	    }
	    if (f=="AbsoluteFractionModern"){
	     f_afn=F@map
             f_d14C=function(t){
	         fd=Delta14C_from_AbsoluteFractionModern(f_afn(t))
	     return(fd)
	    }
	    D14C=F
	    D14C@map=f_d14C
	    D14C@format=targetFormat
	    return(D14C)
	    } 
      stop("conversion not implemented for this format")
    }	 
)
#---------------------------------------------------------------------------------------------------------
setMethod(
  f= "AbsoluteFractionModern",
      signature("BoundFc"),
      definition=function# convert to Absolute Fraction Normal values  
      ### convert a BoundFc object containing values in any supported format to the appropriate Absolute Fraction Modern values.
	    (F ##<< object containing the values in any format
	    ){
        f=F@format
              targetFormat="AbsoluteFractionModern"
              if (f==targetFormat){
	         # do nothing
	         return(F)
	      }
	      if (f=="Delta14C"){
	       f_d14C=F@map
               f_afn=function(t){
	           fprime=AbsoluteFractionModern_from_Delta14C(f_d14C(t))
	       return(fprime)
	       }
	       AFM_tm=F
	       AFM_tm@map=f_afn
	       AFM_tm@format=targetFormat
	       return(AFM_tm)
	      } 
            stop("conversion not implemented for this format")
     }	 
)

setMethod(
    f="getTimeRange",
    signature="BoundFc",
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
    signature="BoundFc",
    definition=function # function definition 
    ### extract the function definition (the R-function) 
    (object){
        return(object@map)
    }
)

