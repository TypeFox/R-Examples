#------------------------------------------------------------------------------
#   This file defines S4 classes that are used in this package, as well as
#   get and set functions.
#   Copyright (C) 2014  Bertram Ieong
#------------------------------------------------------------------------------

setOldClass("data.frame")
setClassUnion("data.frameORvector", c("data.frame", "vector"))

## forest class ##

# define the class
setClass("forest", representation(  n = "integer",
                                    p = "integer",
                                    min_node_obs = "integer",
                                    max_depth = "integer",
                                    numsamps = "integer",
                                    numvars = "integer",
                                    numboots = "integer",
                                    numnodes = "integer",
                                    flattened.nodes = "data.frame",
                                    model = "data.frame",
                                    x = "data.frameORvector",
                                    y = "vector",
                                    fmla = "formula",
                                    depvar.restore.info = "list"
                                    ))

# define the get operator for slots that should be viewable by the user
setMethod(
      f="[",
      signature="forest",
      definition=function(x,i){
            if(i=="n"){return(x@n)}else{}
            if(i=="p"){return(x@p)}else{}
            if(i=="min_node_obs"){return(x@min_node_obs)}else{}
            if(i=="max_depth"){return(x@max_depth)}else{}
            if(i=="numsamps"){return(x@numsamps)}else{}
            if(i=="numvars"){return(x@numvars)}else{}
            if(i=="numboots"){return(x@numboots)}else{}
            if(i=="numnodes"){return(x@numnodes)}else{}
            if(i=="model"){return(x@model)}else{}
            if(i=="x"){return(x@x)}else{}
            if(i=="y"){return(x@y)}else{}
            if(i=="fmla"){return(x@fmla)}else{}
      }
)


## tree class ##

setClass("tree", representation(    n = "integer",
                                    p = "integer",
                                    min_node_obs = "integer",
                                    max_depth = "integer",
                                    numnodes = "integer",
                                    flattened.nodes = "data.frame",
                                    fmla = "formula"
                                    ))

