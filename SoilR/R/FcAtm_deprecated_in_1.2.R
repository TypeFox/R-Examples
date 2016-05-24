#
# vim:set ff=unix expandtab ts=2 sw=2:
FcAtm.from.Dataframe=function
### This function is deprecated constructor of the deprecatied class FcAtm
(dframe, ##<<A data frame containing exactly two columns:
## the first one is interpreted as time
## the secon one is interpreted as atmospheric C14 fraction in the format mentioned
lag=0, ##<< a scalar describing the time lag. Positive Values shift the argument of the interpolation function forward in time. (retard its effect)
interpolation=splinefun, ##<<A function that  returns a function  the default is splinefun. Other possible values are the linear interpolation approxfun or any self made function with the same interface.
format ##<< a string that specifies the format used to represent the atmospheric fracton. Possible values are "Delta14C" which is the default or "afn" the Absolute Fraction Normal representation 
){
   warning("The class FcAtm is deprecated, you can use the generic constructor BoundFc with the same data.frame arguemten instead")
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
   warning(TimeMapWarningBoundFc())
   obj=BoundFc(dframe,lag=lag,format=format) 
return(obj)
### An object of the new class BoundFc that replaces FcAtm 
}
