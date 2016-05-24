#
# vim:set ff=unix expandtab ts=2 sw=2:
correctnessOfConstFc=function#check for unreasonable parameters or unsupported formats
###  14C fraction data can be represented in more than one format 
###  The function checks if the user required format is supported at the moment
(object ##<< the object to be tested
)
{
   res=TRUE
   supported_formats=supported14CFractionFormats()
   f=object@format
   print(paste("format=",f))
   if (!any(grepl(f,supported_formats))){
      err_str=cat("The required format:",f," describing the c_14 fraction is not supported.\n 
   	     The supported formats are: ",supported_formats,". \n",sep="")
      stop(simpleError(err_str))
      return(res)
   }
}


setClass(
   Class="ConstFc",
   representation=representation(
	values="numeric",
	format="character"
   )
)
setMethod(
    f="initialize",
    signature="ConstFc",
    definition=function(.Object,values=numeric(),format="Delta14C"){
    #cat("-initializer at work-\n")
    .Object@values=values
    .Object@format=format
    correctnessOfConstFc(.Object)
    return(.Object)
    }
)
setMethod(
    f="getFormat",
    signature="ConstFc",
    definition=function(# extract the format string
			object ##<< object of class ConstFc containing information aboutn the format that could be Delta14C or AFM (Absolute Fraction Modern) for instance
			){
       ### the function just yields the format as a string
        return(object@format)
    }
)
setMethod(
    f="getValues",
    signature="ConstFc",
    definition=function# extract the format string
			(object ##<< object containing information aboutn the format that could be Delta14C or AFM (Absolute Fraction Modern) for instance
			){
       ### the function just yields the format as a string
        return(object@values)
    }
)
setMethod(
   f= "Delta14C",
      signature("ConstFc"),
      definition=function# convert to Absolute Fraction Normal values  
	(F##<< object containing the values in any format
	){
	### convert a ConstFc object containing values in any supported format to the appropriate Absolute Fraction Modern values.
	f=F@format
        targetFormat="Delta14C"
        if (f==targetFormat){
	   # do nothing
	   return(F)
	}
	if (f=="AbsoluteFractionModern"){
	 f_afn=F@values
	 f_d14C=Delta14C_from_AbsoluteFractionModern(f_afn)
	 D14C=F
	 D14C@values=f_d14C
	 D14C@format=targetFormat
	 return(D14C)
	} 
      stop("conversion not implemented for this format")
      }	 
)

setMethod(
   f= "AbsoluteFractionModern",
      signature("ConstFc"),
      definition=function# convert to Absolute Fraction Normal values  
	(F ##<< object containing the values in any format
	){
	### convert a ConstFc object containing values in any supported format to the appropriate Absolute Fraction Modern values.
	f=F@format
        targetFormat="AbsoluteFractionModern"
        if (f==targetFormat){
	   # do nothing
	   return(F)
	}
	if (f=="Delta14C"){
	 f_d14C=F@values
         f_afn=AbsoluteFractionModern_from_Delta14C(f_d14C)
	 AFM_tm=F
	 AFM_tm@values=f_afn
	 AFM_tm@format=targetFormat
	 return(AFM_tm)
	} 
      stop("conversion not implemented for this format")
      }	 
)
ConstFc <- function # creates an object containing the initial values for the 14C fraction needed to create models in SoilR
    ### The function returns an object of class ConstFc which is a building block for any 14C model in SoilR.
    ### The building blocks of a model have to keep iformation about the formats their data are in, because the high level function dealing wiht the models have to know. This function is actually a convienient wrapper for a call to R's standard constructor new, to hide its complexity from the user.
    (
    values=c(0),  ##<< a numeric vector
    format="Delta14C"   ##<< a character string describing the format e.g. "Delta14C"
    )
    {
    
	F0=new(Class="ConstFc",values=values,format=format)
	return(F0)
	### An object of class ConstFc that contains data and a format description that can later be used to convert the data into other formats if the conversion is implemented.
}
