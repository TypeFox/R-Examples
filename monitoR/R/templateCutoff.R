# Two functions: for extracting and replacing template score cutoffs
# Modified: 6 Sept 2015

templateCutoff <-
function(object) {

   if(!'templates' %in% slotNames(object)) stop('object does not have a templates slot, so score.cutoff cannot be checked or set.') 
   score.cutoff <- as.numeric(lapply(object@templates, function(x) x@score.cutoff))
   names(score.cutoff) <- names(object@templates)
   return(score.cutoff)

}


`templateCutoff<-` <-
function(
   object, 
   value 
) {

   if(length(value) == 1 && names(value) == 'default') {
      for(i in names(object@templates)) {
         object@templates[[i]]@score.cutoff <- as.vector(value) # as.vector to drop names
      }
   } else {
      if(any(names(value) == 'default')) {
         for(i in names(object@templates)) {
            object@templates[[i]]@score.cutoff <- as.vector(value[names(value) == 'default']) 
         }
         value <- value[names(value) != 'default']
      }
      if(is.null(names(value))) names(value) <- names(object@templates)[1:length(value)] else 
         if (!all(names(value) %in% names(object@templates))) 
            stop('Name or names given for new cutoffs (', paste(names(value), collapse=', '), ') are not present in the templates (', paste(names(object@templates), collapse=', '), ').')
      for(i in names(value)) {
         object@templates[[i]]@score.cutoff <- as.vector(value[i]) 
      }
   }

   if(class(object) == "detectionList") object <- findDetections(object, keep.verify=TRUE)
   return(object)

}

