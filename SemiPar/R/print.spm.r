########## R-function: print.spm ##########

# For printing the essence of an results of spm()

# Last changed: 21 JAN 2005 by MPW

print.spm <- function(x,...)
{
   object <- x

   cat("\n This is an spm() fit object; a list with components named:\n\n")

   cat("     ",names(object),"\n\n")

   cat(" The names() function can be used to obtain names \n")
   cat(" of components and sub-components.\n\n")

   cat(" The summary() function can be used to summarise\n")
   cat(" the fit aspects of the object.\n\n")
}

######### End of print.spm ##########










