#---------------------------------------------------------------------------
# CLASS "planordesign"
#  S4 class, typically an output from planor.design.designkey
# SLOTS
#    - design: a dataframe containing the final design
#    - factors: the 'designfactors' object that defines the factors
#    - model: list containing the model and estimate specifications
#    - designkey, nunits, recursive:
# all these slots are data extracted from the  designkey
# from which the object has been built
# Objects of this class are created by calls to planor.design
# METHODS: getDesign, as.data.frame, summary, alias
#---------------------------------------------------------------------------
setClass("planordesign",
         representation(design="data.frame",
                        factors="designfactors",
                        model="list",
                        designkey="list",
                        nunits="numeric",
                        recursive="logical"))
##------------------------------------------------------------------------
## getDesign: Extraction methods for "planordesign"
##------------------------------------------------------------------------
getDesign.planordesign <- function(object){return(object@design)}
setMethod("getDesign", signature(object="planordesign"),
          definition=getDesign.planordesign)
##--------------------------------------------------------------------------

setMethod("[",
          signature(x = "planordesign", i = "ANY", j = "ANY", drop = "ANY"),
          definition=function(x,i,j,...,drop){
            if(missing(i)){ i <- seq(nrow(x@design)) }
            if(missing(j)){ j <- names(x@factors) }
            if(is.character(j)){
              jnames <- names(x@factors) %in% j
              j <- seq(jnames)[jnames]
            }
            if(is.logical(j)){
              j <- seq(j)[j]
            }
            x@design <- x@design[ i, j, drop=FALSE ]
            x@factors <- x@factors[ j ]
            keepmodel <- rep(NA, length(x@model))
            for(m in seq_along(x@model)){
              model.names <- unique(
                         all.vars( rev(as.list(x@model[[m]][[1]]))[[1]] ), # model part
                         all.vars( rev(as.list(x@model[[m]][[2]]))[[1]] )) # estimate part
              keepmodel[m] <- all( model.names %in% j )
            }
            x@model <- x@model[keepmodel]
            x
          })
##
##--------------------------------------------------------------------------
# Method as.data.frame for "planordesign"
# The data.frame is the slot 'design' of the "planordesign" object.
# All the other slots are stored in attributes

as.data.frame.planordesign <- function(x, ...) {
  ret <- x@design
  for (unslot in slotNames(x)) {
    attr(ret, unslot) <- slot(x, unslot)
  }
  return(ret)
} # end as.data.frame.planordesign
##--------------------------------------------------------------------------
# Method summary  for "planordesign"

# Summarises the design properties of a planordesign  object, by
# calling summary.keymatrix on each of its key matrices contained
# in the slot designkey 
# ---------------------------------------------

summary.planordesign <- function(object, fact, block,
                              show="dtbw", save="k", ...){
  onesummary <- function(i, oneobject, ...) {
    cat("\n********** Keymatrix ", i, "**********\n\n")
    return( summary.keymatrix(oneobject, ...))
  }
  sortie <- vector("list", length=length(object@designkey))
  for (i in 1:length(object@designkey)) {
    sortie[[i]] <- onesummary(i, object@designkey[[i]], ...)
  }
   if (save != "")
    return(invisible(sortie))
  else
    return(invisible())

}
##--------------------------------------------------------------------------
# Method alias  for "planordesign"

# Summarises the design properties of a planordesign  object, by
# calling alias on each of its key matrices
# ---------------------------------------------

alias.planordesign <- function(object, model,  fact, block, ...){
  onealias <- function(i, oneobject, ...) {
    cat("\n********** Keymatrix ", i, "**********\n\n")
    return( alias.keymatrix(oneobject, ...))
  }
  sortie <- vector("list", length=length(object@designkey))
  for (i in 1:length(object@designkey)) {
    sortie[[i]] <- onealias(i, object@designkey[[i]], ...)
  }
  return(invisible(sortie))
}
