## f6-graphs.R --- 
## Author          : Jens Henrik Badsberg, Claus Dethlefsen, Luca La Rocca
## Created On      : Fri Jun 24 10:40:00 2005
## Last Modified By: Luca La Rocca
## Last Modified On: Fri Feb 15 14:48:00 2006
## Update Count    : 80
## Status          : Unknown, Use with caution!
######################################################

## construction and visualization

# constructor method for class 'anyGraph'
setMethod("initialize","anyGraph",
          function(.Object,incidenceList){
            if(!missing(incidenceList))
              .Object@incidenceList<-incidenceList
            return(.Object)
          } # end of function
         ) # end of SetMethod
# a valid 'anyGraph' object is returned

# constructor method for class 'generalGraph'
setMethod("initialize","generalGraph",
          function(.Object,incidenceList,incidenceMatrix){
            if(!missing(incidenceMatrix))
              .Object@incidenceMatrix<-incidenceMatrix
            else if(!missing(incidenceList)){
              proper<-unlist(lapply(incidenceList@E,function(y)
                      (is(y,"undirectedEdge")&&(card(y)>0))||(is(y,"directedEdge")&&(length(y)>1))))
              incidenceList@E<-incidenceList@E[proper]
              .Object@incidenceList<-incidenceList
            } # end of if-else-if
            return(.Object)
          } # end of function
         ) # end of SetMethod
# a valid 'generalGraph' object is returned

# constructor method for class 'multiGraph'
setMethod("initialize","multiGraph",
                 function(.Object,incidenceList,incidenceMatrix,adjacencyList){
                     if(!missing(adjacencyList))
                       .Object@adjacencyList<-adjacencyList
                     else if(!missing(incidenceMatrix)){
                       ordinary<-apply(incidenceMatrix@.Data,1,function(y) sum(y>0)<3)
                       .Object@incidenceMatrix<-new("incidenceMatrix",
                                                    matrix(incidenceMatrix@.Data[ordinary,],
                                                           nrow=sum(ordinary),ncol=card(incidenceMatrix)$v,
                                                           dimnames=list(NULL,names(incidenceMatrix))
                                                          ) # end of matrix
                                                   ) # end of new
                     } else if(!missing(incidenceList)){
                       ordinary<-unlist(lapply(incidenceList@E,function(y)
                                 (is(y,"undirectedEdge")&&(card(y)>0)&&(card(y)<3))||
                                 (is(y,"directedEdge")&&(length(y)>1)&&(card(y)<3))))
                       incidenceList@E<-incidenceList@E[ordinary]
                       .Object@incidenceList<-incidenceList
                     } # end of if-else if-else if
                     return(.Object)
                 } # end of function
                ) # end of SetMethod
# a valid 'multiGraph' object is returned

# constructor method for class 'simpleGraph'
setMethod("initialize","simpleGraph",
          function(.Object,incidenceList,incidenceMatrix,adjacencyList,adjacencyMatrix){
            if(!missing(adjacencyMatrix))
              .Object@adjacencyMatrix<-adjacencyMatrix
            else if(!missing(adjacencyList))
              .Object@adjacencyList<-as(as(adjacencyList,"adjacencyMatrix"),"adjacencyList")
            else if(!missing(incidenceMatrix))
              .Object@incidenceMatrix<-as(as(incidenceMatrix,"adjacencyMatrix"),"incidenceMatrix")
            else if(!missing(incidenceList))
              .Object@incidenceList<-as(as(incidenceList,"adjacencyMatrix"),"incidenceList")
            return(.Object)
          } # end of function
         ) # end of SetMethod
# a valid 'simpleGraph' object is returned

# keeping default 'show' methods

# see file 'f8-interfaces.R' for graphical representation methods

## getting and setting information

# see at the end of the file for getting/setting graph representations

# 'names' method for class 'anyGraph'
setMethod("names", "anyGraph", function(x) names(x@incidenceList))
# take the names from the only available representation

# 'names' method for class 'generalGraph'
setMethod("names", "generalGraph",
          function(x){
            if(!isEmpty(x@incidenceMatrix)) names(x@incidenceMatrix)
            else names(x@incidenceList)
          } # end of function
         ) # end of setAs
# take the names from one of the available representations

# 'names' method for class 'multiGraph'
setMethod("names", "multiGraph",
          function(x){
            if(!isEmpty(x@adjacencyList)) names(x@adjacencyList)
            else if(!isEmpty(x@incidenceMatrix)) names(x@incidenceMatrix)
            else names(x@incidenceList)
          } # end of function
         ) # end of setAs
# take the names from one of the available representations

# 'names' method for class 'simpleGraph'
setMethod("names", "simpleGraph",
          function(x){
            if(!isEmpty(x@adjacencyMatrix)) names(x@adjacencyMatrix)
            else if(!isEmpty(x@adjacencyList)) names(x@adjacencyList)
            else if(!isEmpty(x@incidenceMatrix)) names(x@incidenceMatrix)
            else names(x@incidenceList)
          } # end of function
         ) # end of setAs
# take the names from one of the available representations

# 'names<-' replacement method for class 'anyGraph'
setReplaceMethod("names","anyGraph",
                         function(x,value){
                           if(!isEmpty(x@incidenceList)) names(x@incidenceList)<-value
                           x # returns the possibly modified object
                         } # end of function
                        ) # end of setMethod
# 'names<-' replacement method for class 'generalGraph'
setReplaceMethod("names","generalGraph",
                         function(x,value){
                           if(!isEmpty(x@incidenceMatrix)) names(x@incidenceMatrix)<-value
                           if(!isEmpty(x@incidenceList)) names(x@incidenceList)<-value
                           x # returns the possibly modified object
                         } # end of function
                        ) # end of setMethod
# 'names<-' replacement method for class 'multiGraph'
setReplaceMethod("names","multiGraph",
                         function(x,value){
                           if(!isEmpty(x@adjacencyList)) names(x@adjacencyList)<-value
                           if(!isEmpty(x@incidenceMatrix)) names(x@incidenceMatrix)<-value
                           if(!isEmpty(x@incidenceList)) names(x@incidenceList)<-value
                           x # returns the possibly modified object
                         } # end of function
                        ) # end of setMethod
# 'names<-' replacement method for class 'simpleGraph'
setReplaceMethod("names","simpleGraph",
                         function(x,value){
                           if(!isEmpty(x@adjacencyMatrix)) names(x@adjacencyMatrix)<-value
                           if(!isEmpty(x@adjacencyList)) names(x@adjacencyList)<-value
                           if(!isEmpty(x@incidenceMatrix)) names(x@incidenceMatrix)<-value
                           if(!isEmpty(x@incidenceList)) names(x@incidenceList)<-value
                           x # returns the possibly modified object
                         } # end of function
                        ) # end of setMethod
# set the names in non-empty slots

# 'card' method for class 'anyGraph'
setMethod("card", "anyGraph", function(object,...) card(object@incidenceList))
# returns the number of vertices and edges

# 'card' method for class 'generalGraph'
setMethod("card", "generalGraph",
          function(object,...){
            if(!isEmpty(object@incidenceMatrix)) card(object@incidenceMatrix)
            else card(object@incidenceList)
          } #Êend of function
         ) # end of setMethod
# returns the number of vertices and edges

# 'card' method for class 'multiGraph'
setMethod("card", "multiGraph",
          function(object,...){
            if(!isEmpty(object@adjacencyList)) card(object@adjacencyList)
            else if(!isEmpty(object@incidenceMatrix)) card(object@incidenceMatrix)
            else card(object@incidenceList)
          } #Êend of function
         ) # end of setMethod
# returns the number of vertices and edges

# 'card' method for class 'simpleGraph'
setMethod("card", "simpleGraph",
          function(object,...){
            if(!isEmpty(object@adjacencyMatrix)) card(object@adjacencyMatrix)
            else if(!isEmpty(object@adjacencyList)) card(object@adjacencyList)
            else if(!isEmpty(object@incidenceMatrix)) card(object@incidenceMatrix)
            else card(object@incidenceList)
          } #Êend of function
         ) # end of setMethod
# returns the number of vertices and edges

## property checking

# 'isEmpty' method for class 'anyGraph'
setMethod("isEmpty","anyGraph",function(object,...) isEmpty(object@incidenceList))
# a graph object is empty if all its possible representations are empty

# 'isEmpty' method for class 'generalGraph'
setMethod("isEmpty","generalGraph",
          function(object,...) isEmpty(object@incidenceMatrix)&&isEmpty(object@incidenceList)
         ) # end of setMethod
# a graph object is empty if all its possible representations are empty

# 'isEmpty' method for class 'multiGraph'
setMethod("isEmpty","multiGraph",
          function(object,...)
            isEmpty(object@adjacencyList)&&isEmpty(object@incidenceMatrix)&&isEmpty(object@incidenceList)
         ) # end of setMethod
# a graph object is empty if all its possible representations are empty

# 'isEmpty' method for class 'simpleGraph'
setMethod("isEmpty","simpleGraph",
          function(object,...)
            isEmpty(object@adjacencyMatrix)&&isEmpty(object@adjacencyList)&&
            isEmpty(object@incidenceMatrix)&&isEmpty(object@incidenceList)
         ) # end of setMethod
# a graph object is empty if all its possible representations are empty

# 'isPresent' method for 'edge' in 'anyGraph'
setMethod("isPresent",c(el="edge",ou="anyGraph"), function(el,ou) isPresent(el,ou@incidenceList))
# a 'logical' value answering the question is returned

# 'isPresent' method for 'edge' in 'generalGraph'
setMethod("isPresent",c(el="edge",ou="generalGraph"),
          function(el,ou) isPresent(el,ou@incidenceMatrix)||isPresent(el,ou@incidenceList)
         ) # end of setMethod
# a 'logical' value answering the question is returned

# 'isPresent' method for 'edge' in 'multiGraph'
setMethod("isPresent",c(el="edge",ou="multiGraph"),
          function(el,ou)
            isPresent(el,ou@adjacencyList)||isPresent(el,ou@incidenceMatrix)||isPresent(el,ou@incidenceList)
         ) # end of setMethod
# a 'logical' value answering the question is returned

# 'isPresent' method for 'edge' in 'simpleGraph'
setMethod("isPresent",c(el="edge",ou="simpleGraph"),
          function(el,ou)
            isPresent(el,ou@adjacencyMatrix)||isPresent(el,ou@adjacencyList)||
            isPresent(el,ou@incidenceMatrix)||isPresent(el,ou@incidenceList)
         ) # end of setMethod
# a 'logical' value answering the question is returned

# comparison method for class 'anyGraph'
setMethod("areTheSame",c("anyGraph","anyGraph"),
                 function(x,y) areTheSame(incidenceList(x),incidenceList(y))
         ) # end of setMethod
# a 'logical' value answering the question is returned

# comparison methods for classes 'generalGraph', 'multiGraph' and 'simpleGraph' follow by inheritance

## extraction

# multi extractor method for class 'anyGraph'
setMethod("[","anyGraph",
          function(x,i,j=NA,drop=NA){
            x@incidenceList<-x@incidenceList[i]
            return(x)
          } # end of function
         ) # end of setMethod
# the subgraph induced by 'i' is extracted

# single extractor method for class 'anyGraph'
setMethod("[[","anyGraph",function(x,i,j=NA,drop=NA) names(x)[i])
# the name of the i-th vertex is extracted

# multi extractor method for class 'generalGraph'
setMethod("[","generalGraph",
          function(x,i,j=NA,drop=NA){
            x@incidenceMatrix<-x@incidenceMatrix[i]
            x@incidenceList<-x@incidenceList[i]
            return(x)
          } # end of function
         ) # end of setMethod
# the subgraph induced by 'i' is extracted

# single extractor method for class 'generalGraph'
setMethod("[[","generalGraph",function(x,i,j=NA,drop=NA) names(x)[i])
# the name of the i-th vertex is extracted

# multi extractor method for class 'multiGraph'
setMethod("[","multiGraph",
          function(x,i,j=NA,drop=NA){
            x@adjacencyList<-x@adjacencyList[i]
            x@incidenceMatrix<-x@incidenceMatrix[i]
            x@incidenceList<-x@incidenceList[i]
            return(x)
          } # end of function
         ) # end of setMethod
# the subgraph induced by 'i' is extracted

# single extractor method for class 'multiGraph'
setMethod("[[","multiGraph",function(x,i,j=NA,drop=NA) names(x)[i])
# the name of the i-th vertex is extracted

# multi extractor method for class 'simpleGraph'
setMethod("[","simpleGraph",
          function(x,i,j=NA,drop=NA){
            x@adjacencyMatrix<-x@adjacencyMatrix[i]
            x@adjacencyList<-x@adjacencyList[i]
            x@incidenceMatrix<-x@incidenceMatrix[i]
            x@incidenceList<-x@incidenceList[i]
            return(x)
          } # end of function
         ) # end of setMethod
# the subgraph induced by 'i' is extracted

# single extractor method for class 'simpleGraph'
setMethod("[[","simpleGraph",function(x,i,j=NA,drop=NA) names(x)[i])
# the name of the i-th vertex is extracted

## typecasting

# see below for typecasting between graphs

# see file 'f8-interfaces.R' for typecasting to/from classes
# of other related packages ("mathgaph" and "dynamicGraph")

## operators

# see file 'f7-operators.R' for '+/-/*' methods

## -----------------------------------------------------------
## Typecasting between graphs
## -----------------------------------------------------------

## from bottom to top

# typecasting from 'simpleGraph' to 'multiGraph'
setAs("simpleGraph","multiGraph",
      function(from,to){
        res<-new(to)
        if(!isEmpty(from@adjacencyMatrix))
          res@adjacencyList<-as(from@adjacencyMatrix,"adjacencyList")
        else
          res@adjacencyList<-from@adjacencyList
        res@incidenceMatrix<-from@incidenceMatrix
        res@incidenceList<-from@incidenceList
        return(res)
      } # end of function
     ) # end of setAs
# a 'multiGraph' object is returned

# typecasting from 'simpleGraph' to 'generalGraph'
setAs("simpleGraph","generalGraph",
      function(from,to){
        res<-new(to)
        if(!isEmpty(from@adjacencyMatrix))
          res@incidenceMatrix<-as(from@adjacencyMatrix,"incidenceMatrix")
        else if(!isEmpty(from@adjacencyList))
          res@incidenceMatrix<-as(from@adjacencyList,"incidenceMatrix")
        else
          res@incidenceMatrix<-from@incidenceMatrix
        res@incidenceList<-from@incidenceList
        return(res)
      } # end of function
     ) # end of setAs
# a 'generalGraph' object is returned

# typecasting from 'simpleGraph' to 'anyGraph'
setAs("simpleGraph","anyGraph",
      function(from,to){
        res<-new(to)
        if(!isEmpty(from@adjacencyMatrix))
          res@incidenceList<-as(from@adjacencyMatrix,"incidenceList")
        else if(!isEmpty(from@adjacencyList))
          res@incidenceList<-as(from@adjacencyList,"incidenceList")
        else if(!isEmpty(from@incidenceMatrix))
          res@incidenceList<-as(from@incidenceMatrix,"incidenceList")
        else
          res@incidenceList<-from@incidenceList
        return(res)
      } # end of function
     ) #Êend of setAs
# an 'anyGraph' object is returned

# typecasting from 'multiGraph' to 'generalGraph'
setAs("multiGraph","generalGraph",
      function(from,to){
        res<-new(to)
        if(!isEmpty(from@adjacencyList))
          res@incidenceMatrix<-as(from@adjacencyList,"incidenceMatrix")
        else
          res@incidenceMatrix<-from@incidenceMatrix
        res@incidenceList<-from@incidenceList
        return(res)
      } # end of function
     ) # end of setAs
# a 'generalGraph' object is returned

# typecasting from 'multiGraph' to 'anyGraph'
setAs("multiGraph","anyGraph",
      function(from,to){
        res<-new(to)
        if(!isEmpty(from@adjacencyList))
          res@incidenceList<-as(from@adjacencyList,"incidenceList")
        else if(!isEmpty(from@incidenceMatrix))
          res@incidenceList<-as(from@incidenceMatrix,"incidenceList")
        else
          res@incidenceList<-from@incidenceList
        return(res)
      } # end of function
     ) # end of setAs
# an 'anyGraph' object is returned

# typecasting from 'generalGraph' to 'anyGraph'
setAs("generalGraph","anyGraph",
      function(from,to){
        res<-new(to)
        if (!isEmpty(from@incidenceMatrix))
          res@incidenceList<-as(from@incidenceMatrix,"incidenceList")
        else
          res@incidenceList<-from@incidenceList
        return(res)
      } # end of function
     ) #Êend of setAs
# an 'anyGraph' object is returned

## from top to bottom

# typecasting from 'anyGraph' to 'generalGraph'
setAs("anyGraph","generalGraph",
      function(from,to){
        warning("Coercing anyGraph to generalGraph, possibly loosing information...")
        new(to,incidenceList=from@incidenceList)
      } # end of function
     ) #Êend of setAs
# a 'generalGraph' object is returned

# typecasting from 'anyGraph' to 'multiGraph'
setAs("anyGraph","multiGraph",
      function(from,to){
        warning("Coercing anyGraph to multiGraph, possibly loosing information...")
        new(to,incidenceList=from@incidenceList)
      } #Êend of function
     ) # end of setAs
# a 'multiGraph' object is returned

# typecasting from 'anyGraph' to 'simpleGraph'
setAs("anyGraph","simpleGraph",
      function(from,to){
        warning("Coercing anyGraph to simpleGraph, possibly loosing information...")
        new(to,incidenceList=from@incidenceList)
      } # end of function
     ) # end of setAs
# a 'simpleGraph' object is returned

# typecasting from 'generalGraph' to 'multiGraph'
setAs("generalGraph","multiGraph",
      function(from,to){
        warning("Coercing generalGraph to multiGraph, possibly loosing information...")
        if(!isEmpty(from@incidenceList)){ # build from incidence list
          res<-new(to,incidenceList=from@incidenceList)
          if(!isEmpty(from@incidenceMatrix))
            res@incidenceMatrix<-as(res@incidenceList,"incidenceMatrix")
        }else{ # build from incidence matrix
          res<-new(to,incidenceMatrix=from@incidenceMatrix)
        } # end of if-else
        return(res)
      } # end of function
     ) # end of setAs
# a 'multiGraph' object is returned

# typecasting from 'generalGraph' to 'simpleGraph'
setAs("generalGraph","simpleGraph",
      function(from,to){
        warning("Coercing generalGraph to simpleGraph, possibly loosing information...")
        if(!isEmpty(from@incidenceList)){ # build from incidence list
          res<-new(to,incidenceList=from@incidenceList)
          if(!isEmpty(from@incidenceMatrix))
            res@incidenceMatrix<-as(res@incidenceList,"incidenceMatrix")
        }else{ # build from incidence matrix
          res<-new(to,incidenceMatrix=from@incidenceMatrix)
        } # end of if-else
        return(res)
      } # end of function
     ) # end of setAs
# a 'simpleGraph' object is returned

# typecasting from 'multiGraph' to 'simpleGraph'
setAs("multiGraph","simpleGraph",
      function(from,to){
        warning("Coercing multiGraph to simpleGraph, possibly loosing information...")
        if(!isEmpty(from@incidenceList)){ # build from incidence list
          res<-new(to,incidenceList=from@incidenceList)
          if(!isEmpty(from@incidenceMatrix))
            res@incidenceMatrix<-as(res@incidenceList,"incidenceMatrix")
          if(!isEmpty(from@adjacencyList))
            res@adjacencyList<-as(res@incidenceList,"adjacencyList")
        }else if(!isEmpty(from@incidenceMatrix)){ # build from incidence matrix
          res<-new(to,incidenceMatrix=from@incidenceMatrix)
          if(!isEmpty(from@adjacencyList))
            res@adjacencyList<-as(res@incidenceMatrix,"adjacencyList")
        }else{ # build from adjacency list
          res<-new(to,adjacencyList=from@adjacencyList)
        } # end of if-else if-else
        return(res)
      } # end of function
     ) # end of setAs
# a 'simpleGraph' object is returned

## ----------------------------------------------------------
## Getting graph representations
## ----------------------------------------------------------

# getting the adjacency matrix of a simple graph
setMethod("adjacencyMatrix", "simpleGraph",
          function(object, ...) {
          if (!isEmpty(object@adjacencyMatrix))
            object@adjacencyMatrix
          else if (!isEmpty(object@adjacencyList))
            as(object@adjacencyList,"adjacencyMatrix")
          else if (!isEmpty(object@incidenceMatrix))
            as(object@incidenceMatrix,"adjacencyMatrix")
          else
            as(object@incidenceList,"adjacencyMatrix")
          } # end of function
         ) # end of setMethod
# an 'adjacencyMatrix' object is returned

# getting the adjacency list of a simple graph
setMethod("adjacencyList", "simpleGraph",
          function(object, ...) {
          if (!isEmpty(object@adjacencyList))
            object@adjacencyList
          else if (!isEmpty(object@adjacencyMatrix))
            as(object@adjacencyMatrix,"adjacencyList")
          else if (!isEmpty(object@incidenceMatrix))
            as(object@incidenceMatrix,"adjacencyList")
          else
            as(object@incidenceList,"adjacencyList")
          } # end of function
         ) # end of setMethod
# an 'adjacencyList' object is returned

# getting the incidence matrix of a simple graph
setMethod("incidenceMatrix", "simpleGraph",
          function(object, ...) {
          if (!isEmpty(object@incidenceMatrix))
            object@incidenceMatrix
          else if (!isEmpty(object@adjacencyList))
            as(object@adjacencyList,"incidenceMatrix")
          else if (!isEmpty(object@adjacencyMatrix))
            as(object@adjacencyMatrix,"incidenceMatrix")
          else
            as(object@incidenceList,"incidenceMatrix")
          } # end of function
         ) # end of setMethod
# an 'incidenceMatrix' object is returned

# getting the incidence list of a simple graph
setMethod("incidenceList", "simpleGraph",
          function(object, ...) {
          if (!isEmpty(object@incidenceList))
            object@incidenceList
          else if (!isEmpty(object@incidenceMatrix))
            as(object@incidenceMatrix,"incidenceList")
          else if (!isEmpty(object@adjacencyList))
            as(object@adjacencyList,"incidenceList")
          else
            as(object@adjacencyMatrix,"incidenceList")
          } # end of function
         ) # end of setMethod
# an 'incidenceList' object is returned

# getting the adjacency list of a multi graph
setMethod("adjacencyList", "multiGraph",
          function(object, ...) {
            if (!isEmpty(object@adjacencyList))
              object@adjacencyList
            else if (!isEmpty(object@incidenceMatrix))
              as(object@incidenceMatrix,"adjacencyList")
            else
              as(object@incidenceList,"adjacencyList")
          } # end of function
         ) # end of setMethod
# an 'adjacencyList' object is returned

# getting the incidence matrix of a multi graph
setMethod("incidenceMatrix", "multiGraph",
          function(object, ...) {
            if (!isEmpty(object@incidenceMatrix))
              object@incidenceMatrix
            else if (!isEmpty(object@adjacencyList))
              as(object@adjacencyList,"incidenceMatrix")
            else
              as(object@incidenceList,"incidenceMatrix")
          } # end of function
         ) # end of setMethod
# an 'incidenceMatrix' object is returned

# getting the incidence list of a multi graph
setMethod("incidenceList", "multiGraph",
          function(object, ...) {
            if (!isEmpty(object@incidenceList))
              object@incidenceList
            else if (!isEmpty(object@incidenceMatrix))
              as(object@incidenceMatrix,"incidenceList")
            else
              as(object@adjacencyList,"incidenceList")
          } # end of function
         ) # end of setMethod
# an 'incidenceList' object is returned

# getting the incidence matrix of a general graph
setMethod("incidenceMatrix", "generalGraph",
          function(object, ...) {
            if (!isEmpty(object@incidenceMatrix))
              object@incidenceMatrix
            else
              as(object@incidenceList,"incidenceMatrix")
          } # end of function
         ) # end of setMethod
# an 'incidenceMatrix' object is returned

# getting the incidence list of a general graph
setMethod("incidenceList", "generalGraph",
          function(object, ...) {
            if (!isEmpty(object@incidenceList))
              object@incidenceList
            else
              as(object@incidenceMatrix,"incidenceList")
          } # end of function
         ) # end of setMethod
# an 'incidenceList' object is returned

# getting the incidence list of any graph
setMethod("incidenceList", "anyGraph", function(object, ...) object@incidenceList)
# an 'incidenceList' object is returned

## ----------------------------------------------------------
## Setting graph representations (replacement methods)
## ----------------------------------------------------------

# setting the adjacency matrix of a simple graph
setReplaceMethod("adjacencyMatrix", "simpleGraph",
                 function(x, force=TRUE,value){
                   if (force) {
                     x@adjacencyMatrix <- value
                     x@incidenceList <- new("incidenceList")
                     x@incidenceMatrix <- new("incidenceMatrix")
                     x@adjacencyList <- new("adjacencyList")
                   } else if(areTheSame(adjacencyMatrix(x),value))
                       x@adjacencyMatrix <- value
                     else
                       warning("Not matching adjacency matrix...")
                   return(x)
                 } # end of function
                ) #Êend of setReplaceMethod
# returns a graph object of the same class as the input object

# setting the adjacency list of a multi graph
setReplaceMethod("adjacencyList", "multiGraph",
                 function(x, force = TRUE, value){
                   if (force) {
                     x@adjacencyList <- value
                     x@incidenceList <- new("incidenceList")
                     x@incidenceMatrix <- new("incidenceMatrix")
                   } else if(areTheSame(adjacencyList(x),value))
                       x@adjacencyList <- value
                     else
                       warning("Not matching adjacency list...")
                   return(x)
                 } # end of function
                ) #Êend of setReplaceMethod
# returns a graph object of the same class as the input object

# setting the adjacency list of a simple graph
setReplaceMethod("adjacencyList", "simpleGraph",
                 function(x, force = TRUE, value){
                   if (force) {
                     x@adjacencyList <- value
                     x@incidenceList <- new("incidenceList")
                     x@incidenceMatrix <- new("incidenceMatrix")
                     x@adjacencyMatrix <- new("adjacencyMatrix")
                   } else if(areTheSame(adjacencyList(x),value))
                       x@adjacencyList <- value
                     else
                       warning("Not matching adjacency list...")
                   return(x)
                 } # end of function
                ) #Êend of setReplaceMethod
# returns a graph object of the same class as the input object

# setting the incidence matrix of a general graph
setReplaceMethod("incidenceMatrix", "generalGraph",
                 function(x, force = TRUE, value){
                   if (force) {
                     x@incidenceMatrix <- value
                     x@incidenceList <- new("incidenceList")
                   } else if(areTheSame(incidenceMatrix(x),value))
                       x@incidenceMatrix <- value
                     else
                       warning("Not matching incidence matrix...")
                   return(x)
                 } # end of function
                ) #Êend of setReplaceMethod
# returns a graph object of the same class as the input object

# setting the incidence matrix of a multi graph
setReplaceMethod("incidenceMatrix", "multiGraph",
                 function(x, force = TRUE, value){
                   if (force) {
                     x@incidenceMatrix <- value
                     x@incidenceList <- new("incidenceList")
                     x@adjacencyList <- new("adjacencyList")
                   } else if(areTheSame(incidenceMatrix(x),value))
                       x@incidenceMatrix <- value
                     else
                       warning("Not matching incidence matrix...")
                   return(x)
                 } # end of function
                ) #Êend of setReplaceMethod
# returns a graph object of the same class as the input object

# setting the incidence matrix of a simple graph
setReplaceMethod("incidenceMatrix", "simpleGraph",
                 function(x, force = TRUE, value){
                   if (force) {
                     x@incidenceMatrix <- value
                     x@incidenceList <- new("incidenceList")
                     x@adjacencyList <- new("adjacencyList")
                     x@adjacencyMatrix <- new("adjacencyMatrix")
                   } else if(areTheSame(incidenceMatrix(x),value))
                       x@incidenceMatrix <- value
                     else
                       warning("Not matching incidence matrix...")
                   return(x)
                 } # end of function
                ) #Êend of setReplaceMethod
# returns a graph object of the same class as the input object

# setting the incidence list of any graph
setReplaceMethod("incidenceList", "anyGraph",
                 function(x, force = TRUE, value){
                   if (force) {
                     x@incidenceList <- value
                   } else if(areTheSame(incidenceList(x),value))
                       x@incidenceList <- value
                     else
                       warning("Not matching incidence list...")
                   return(x)
                 } # end of function
                ) #Êend of setReplaceMethod
# returns a graph object of the same class as the input object

# setting the incidence list of a general graph
setReplaceMethod("incidenceList", "generalGraph",
                 function(x, force = TRUE, value){
                   if (force) {
                     x@incidenceList <- value
                     x@incidenceMatrix <- new("incidenceMatrix")
                   } else if(areTheSame(incidenceList(x),value))
                       x@incidenceList <- value
                     else
                       warning("Not matching incidence list...")
                   return(x)
                 } # end of function
                ) #Êend of setReplaceMethod
# returns a graph object of the same class as the input object

# setting the incidence list of a multi graph
setReplaceMethod("incidenceList", "multiGraph",
                 function(x, force = TRUE, value){
                   if (force) {
                     x@incidenceList <- value
                     x@incidenceMatrix <- new("incidenceMatrix")
                     x@adjacencyList <- new("adjacencyList")
                   } else if(areTheSame(incidenceList(x),value))
                       x@incidenceList <- value
                     else
                       warning("Not matching incidence list...")
                   return(x)
                 } # end of function
                ) #Êend of setReplaceMethod
# returns a graph object of the same class as the input object

# setting the incidence list of a simple graph
setReplaceMethod("incidenceList", "simpleGraph",
                 function(x, force = TRUE, value){
                   if (force) {
                     x@incidenceList <- value
                     x@incidenceMatrix <- new("incidenceMatrix")
                     x@adjacencyList <- new("adjacencyList")
                     x@adjacencyMatrix <- new("adjacencyMatrix")
                   } else if(areTheSame(incidenceList(x),value))
                       x@incidenceList <- value
                     else
                       warning("Not matching incidence list...")
                   return(x)
                 } # end of function
                ) #Êend of setReplaceMethod
# returns a graph object of the same class as the input object
