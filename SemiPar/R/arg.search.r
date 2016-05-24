########## S function: arg.search ##########

# For searching for an argument assignment

# Last changed: 16 NOV 2004

arg.search <- function(string,arg.name)
{

   out <- break.string(string,arg.name)

   # Check that arg.name appears in 
   # string

   if (length(out)==1) present <- FALSE
   if (length(out)>1) present <- TRUE

   # Check for more than one occurence of
   # the argument name

   if (length(out)>2) stop("more than one call to this argument")

   # Now extract argument assignment.

   if (present)
   {
      right.string <- out[2]

      type.arg <- "ordinary"
      comma.found <- FALSE
      for (i in 1:nchar(right.string))
      {
         if (substring(right.string,i,i)==",") comma.found <- TRUE
         if ((substring(right.string,i,i)=="(") & comma.found==FALSE)
            type.arg <- "array"
      }

      if (type.arg=="ordinary")
      {
         out.comma <- break.string(right.string,",")
         arg.assign <- paste(arg.name,out.comma[1],sep="")
      }

      if (type.arg=="array")
      {
         out.left <- break.string(right.string,"(")
         out <- break.string(right.string,")")
         arg.assign <- paste(arg.name,out[1],")",sep="")
      }
   }  
  
   if (!present) arg.assign <- NULL

   return(list(arg=arg.assign,present=present))
}

########## End of arg.search ##########


