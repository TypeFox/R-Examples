# find if a structure has an element called field
isfield  <- function(structure, field){
  any( field %in% names(structure) )
} 
