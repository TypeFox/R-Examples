########## S function: rm.char ##########

# For deleting a specified character from
# a string.

# Last changed: 27/05/98

rm.char <- function(string,char)
{
   return(paste(break.string(string,char),collapse=""))
}

########## End of rm.char ##########
