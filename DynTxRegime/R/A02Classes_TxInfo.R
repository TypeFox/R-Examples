#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                                 CLASS TxInfo                                 #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains information regarding the treatment variable for a single dp        #
#                                                                              #
# ptsSubset  : a character vector, the ith element is the nickname of the tx   #
#              subset available to the ith patient                             #
#                                                                              #
# subsetRule : the original fSet function provided by the user                 #
#                                                                              #
# subsets    : a named list of the tx subsets                                  #
#                                                                              #
# superSet   : a character vector of all possible tx options                   #
#                                                                              #
# txName     : column header of data.frame that contains the tx variable       #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

if(!isClass("function or NULL")){
  setClassUnion("function or NULL", 
                members = c("function","NULL"))
}

if(!isClass("integer or factor")){
  setClassUnion("integer or factor", 
                members = c("integer","factor"))
}

setClass(Class = "TxInfo",
         slots = c( ptsSubset = 'character',
                   subsetRule = 'function or NULL',
                      subsets = 'list',
                     superSet = 'character',
                       txName = 'character') )

setClass(Class = "TxInfoList",
         contains = "List" )

if(!isClass("TxInfo or TxInfoList")){
  setClassUnion("TxInfo or TxInfoList", 
                members = c("TxInfo","TxInfoList"))
}

setMethod(f = "PtsSubset",
          signature = c(object = "TxInfo"), 
          definition = function(object, ...){ return( object@ptsSubset ) } )
          
setMethod(f = "SubsetRule",
          signature = c(object = "TxInfo"), 
          definition = function(object, ...){ return( object@subsetRule ) } )

setMethod(f = "Subsets", 
          signature = c(object = "TxInfo"), 
          definition = function(object, ...){ return( object@subsets ) } )
          
setMethod(f = "SuperSet", 
          signature = c(object = "TxInfo"), 
          definition = function(object, ...){ return( object@superSet ) } )
          
setMethod(f = "TxName", 
          signature = c(object = "TxInfo"),
          definition = function(object, ...){ return( object@txName ) } )


