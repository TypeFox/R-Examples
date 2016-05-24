setGeneric("clearMemoryManagement",
function(node, ...)
{
  standardGeneric("clearMemoryManagement")
})

setMethod("clearMemoryManagement", "XMLInternalElementNode",
function(node, ...)
{
  .Call("R_clearNodeMemoryManagement", node, PACKAGE = "XML")
})


manageMemory_p =
function(finalizer)
{
   if(is.character(finalizer) || is(finalizer, "externalptr") ||
        inherits(finalizer, c("NativeSymbol", "NativeSymbolInfo")))
      return(TRUE)

   as.logical(finalizer)
 }
