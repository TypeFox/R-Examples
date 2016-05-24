# These are classes and facilities for referring to a DTD for the
# DOCTYPE field of an XML document


# The 4 elements are + or -//creator//name of what is being referenced//language (decribed by ISO639)
# See XML Elements of Style by Simon St. Laurent.

validatePublicIdentifier =
function(object)
{
    els = strsplit(object, "//")[[1]]
    if(length(els) != 4)
      return("a PUBLIC identifier must have 4 parts, separated by //")
    if(! (els[1] %in%  c("+", "-")))
      return("first element of PUBLIC identifier must be + or -")

    TRUE
}  

setClass("DTDPublicIdentifier",
         contains = "character",
         validity = validatePublicIdentifier)

 # name is the node name for the top-level node.
 # system and public identify the DTD.
setClass("Doctype", representation(name = "character",
                                   system = "character",
                                   public = "character"),
                    validity = function(object) {
                        if(length(nchar(object@system)) > 0 && length(object@public) > 0)
                           return("only one of system and public can be specified")

                        if(length(object@public) > 0 && length(object@public) != 2)
                          return("the public part of the Doctype must have exactly 2 elements.")

                        if(length(object@public) > 0) {
                           tmp = validatePublicIdentifier(object@public[1])
                           if(!is.logical(tmp))
                             return(tmp)
                         }
                        
                        TRUE
                      })

Doctype =
function(system = character(), public = character(), name = "")
{
  if(length(public) == 1 && length(system) > 0) {
      public = c(public, system)
      system = character()
  }
  new("Doctype", name = name, system = system, public = public)
}

ddQuote =
function(x)
{
  if(length(x) == 0)
    return(character())

  paste('"', x, '"', sep = "")
}  

setAs("Doctype", "character",
       function(from) {

         extra = character()
         if(sum(nchar(from@public), nchar(from@system))) {
            if(length(from@system))
              extra = c(extra, "SYSTEM", ddQuote(from@system))
            if(length(from@public))
              extra = c(extra, "PUBLIC", ddQuote(from@public))
          }
         paste("<!DOCTYPE", from@name, paste(extra, collapse = " "), ">")
       })



