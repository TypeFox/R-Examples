# Two functions: for extracting and replacing template paths
# Modified: 6 Sept 2015

templatePath <-
function(object) {

   if(!'templates' %in% slotNames(object)) stop('object does not have a templates slot, so path cannot be checked or set.') 
   clip.path <- as.character(lapply(object@templates, function(x) x@clip.path))
   names(clip.path) <- names(object@templates)
   return(clip.path)

}


`templatePath<-` <-
function(
   object, 
   value 
) {

   if(length(value) == 1 && names(value) == 'default') {
      for(i in names(object@templates)) {
         object@templates[[i]]@clip.path <- as.vector(value) # as.vector to drop names
      }
   } else {
      if(any(names(value) == 'default')) {
         for(i in names(object@templates)) {
            object@templates[[i]]@clip.path <- as.vector(value[names(value) == 'default']) 
         }
         value <- value[names(value) != 'default']
      }
      if(is.null(names(value))) names(value) <- names(object@templates)[1:length(value)] else 
         if (!all(names(value) %in% names(object@templates))) 
            stop('Name or names given for new paths (', paste(names(value), collapse=', '), ') are not present in the templates (', paste(names(object@templates), collapse=', '), ').')
      for(i in names(value)) {
         object@templates[[i]]@clip.path <- as.vector(value[i]) 
      }
   }

   return(object)

}

