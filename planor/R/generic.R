
##---------------------------------------------------------------------------
# GENERIC FUNCTIONS
##---------------------------------------------------------------------------

# ---------------------------------------------------
# Generic "bind"
# ---------------------------------------------------
#  Generic function to bind two objects.
# ARGUMENTS
#  - x: first object
#  - y: second object
#  - ...: other arguments
# RETURN
#  An object where the input are binded.
# --------------------------------------
setGeneric("bind",
           function(x,y,...){ value <- standardGeneric("bind") },
           useAsDefault=FALSE)
## --------------------------------------------
## Generic "pick"
## --------------------------------------------
#  Generic function to extract one object from a complex object
# ARGUMENTS
#  - keys : an object from which to pick
#  - selection:  index vector to indicate the selection
# RETURN
# The selected object
## --------------------------------------------------------
setGeneric("pick",
           function(keys,...){
             value <- standardGeneric("pick")
         })
## --------------------------------------------
# Generic "planor.design" 
# --------------------------------------
#  Generic function to build a design
# ARGUMENTS
#   - key : an object from which the design will be built
# RETURN
# An object which  contains the design build from the input.
# --------------------------------------
setGeneric("planor.design",
           function(key, ...){
               value <- standardGeneric("planor.design")
           })
##---------------------------------------------------------------------------
setGeneric("getDesign",
           function(object, ...){
               value <- standardGeneric("getDesign")
           })

setGeneric("alias",
           function(object, model, ...){
             value <- standardGeneric("alias")
             })
