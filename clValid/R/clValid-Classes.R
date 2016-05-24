####################################################################
## clValid Class Definitions
####################################################################

setClassUnion("numeric or NULL", c("numeric", "NULL"))
setClassUnion("array or NULL", c("array", "NULL"))
setClassUnion("character or array or list or NULL or logical",
              c("character","array","list","NULL","logical"))
setClass("clValid",representation(clusterObjs="list",measures="array",
                                  measNames="character",clMethods="character",
                                  labels="character",
                                  nClust="numeric",validation="character",
                                  metric="character",method="character",neighbSize="numeric",
                                  annotation="character or array or list or NULL or logical",
                                  GOcategory="character",
                                  goTermFreq="numeric",
                                  call="call"))
