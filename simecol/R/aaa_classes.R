############################################################
#  The S4 classes                                          #
############################################################

## helper classes
setClassUnion("functionOrNULL", c("NULL", "function"))
setClassUnion("functionOrcharacter", c("character", "function"))
setClassUnion("listOrNULL", c("NULL", "list"))
setClassUnion("numericOrlist", c("numeric", "list"))
setClassUnion("listOrdata.frame", c("list", "data.frame"))

## main classes of simecol
setClass("simObj",
         representation(
           main = "function",
           equations = "listOrNULL",
           times     = "numeric",
           init      = "ANY",
           parms     = "ANY",
           inputs    = "ANY",
           solver    = "functionOrcharacter",
           observer  = "functionOrNULL",
           out       = "ANY",
           initfunc  = "functionOrNULL"
         )
)

setClass("odeModel",
         representation(
           parms  = "numericOrlist",
           init   = "numeric",
           observer ="NULL" # observer not possible for ODE models
         ),
         contains = "simObj"
)

setClass("gridModel",
         representation(
           parms  = "list",
           init   = "matrix"
         ),
         contains = "simObj"
)

setClass("rwalkModel",
         representation(
           parms  = "list",
           init   = "ANY" # or dataframeOrMatrix
         ),
         contains = "simObj"
)

setClass("indbasedModel",
         representation(
           parms  = "list",
           init   = "listOrdata.frame"
         ),
         contains = "simObj"
)
