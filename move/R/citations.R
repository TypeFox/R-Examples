setGeneric("citations", function(obj) standardGeneric("citations"))
setMethod("citations", ".MoveGeneral", function(obj){
  return(obj@citation)
})

setGeneric("citations<-", function(obj, value) standardGeneric("citations<-"))
setReplaceMethod("citations", ".MoveGeneral", function(obj, value){
  if (length(value)!=1) {
    value <- as.character(value[1])
    warning("There were more than one citation entered. Only using first element")
  }
  slot(obj,'citation', check=T) <- value
  obj
})
