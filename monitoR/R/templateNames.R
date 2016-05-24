# Why don't we just make this a method for the generic names function?
# Two functions for extracting and replacing template names
# Modified: 6 Sept 2015

templateNames <- function(object) return(names(object@templates)) 

`templateNames<-` <-
function(
   object=NULL, 
   value
) {

   if(length(object@templates) == length(value)) 
      names(object@templates) <- value else
      stop('Length of names argument, ', length(value), ', does not match number of templates, ', length(object@templates))

   return(object)

}



