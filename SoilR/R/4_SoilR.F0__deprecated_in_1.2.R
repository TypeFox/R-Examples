#
# vim:set ff=unix expandtab ts=2 sw=2:
 SoilR.F0.new <- function # creates an object of the deprecated class SoilR.F0 

     (
     values=c(0),  ##<< a numeric vector
     format="Delta14C"   ##<< a character string describing the format e.g. "Delta14C"
     )
     {
      warning(WarningConstFc())
 	F0=ConstFc(values=values,format=format)
 	return(F0)
 	### An object of class SoilR.F0 that contains data and a format description that can later be used to convert the data into other formats if the conversion is implemented.
 }
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
