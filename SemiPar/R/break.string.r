########## S function: break.string ##########

# For breaking up a string that is delimited
# by the character signified by "sep"

# Last changed: 22/05/98

break.string <- function(string,sep=" ")
{
   # Obtain the length of the separating string

   sep.len <- nchar(sep)

   # Next determine the starting positions of the separators.

   sep.loc <- NULL

   for (i in 1:nchar(string))
      if (substring(string,i,(i+sep.len-1))==sep) sep.loc <- c(sep.loc,i)

   numwords <- length(sep.loc) + 1
   out <- rep(0,numwords)

   if (numwords==1)
      out[1] <- string
   else
   {
      out[1] <- substring(string,1,(sep.loc[1]-1))
      if (numwords > 2) 
      for (i in 2:(numwords-1))
         out[i] <- substring(string,(sep.loc[i-1]+sep.len),(sep.loc[i]-1))
   
      out[numwords] <- substring(string,(sep.loc[numwords-1]+sep.len ),
                                 nchar(string))
                                 
   }

   return(out)
}

########## End of break.string ##########
