## f2-genericFunctions.R --- 
## Author          : Jens Henrik Badsberg, Claus Dethlefsen, Luca La Rocca
## Created On      : Tue Nov 30 14:23:00 2004
## Last Modified By: Claus Dethlefsen
## Last Modified On: Sun Dec 17 21:38:36 2006
## Update Count    : 26
## Status          : Unknown, Use with caution!
######################################################

# showing an object relative to a given code
if (!isGeneric("showRel")) {
    setGeneric("showRel",function(object,code) standardGeneric("showRel"))
} # end of if

###

# getting the character identifiers of an object
#if (!isGeneric("names")) {
#    setMethod("names",function(x) standardGeneric("names"))
#} # end of if

# setting the character identifiers of an object
#if (!isGeneric("names<-")) {
#    setGeneric("names<-",function(x,value) standardGeneric("names<-"))
#  } # end of if

###

# getting the cardinality of an object
if (!isGeneric("card")) {
    setGeneric("card",function(object,...) standardGeneric("card"))
} # end of if

# 'card' method for class 'vector'
setMethod("card","vector",function(object,...){length(object)})
# returns the length of the vector

###

# checking whether an object is empty
if (!isGeneric("isEmpty")) {
    setGeneric("isEmpty",function(object,...) standardGeneric("isEmpty"))
} # end of if

# isEmpty method for class 'vector'
setMethod("isEmpty","vector",function(object,...){length(object)==0})
# a 'logical' value answering the question is returned

# 'isEmpty' method for class 'NULL'
setMethod("isEmpty","NULL",function(object,...) TRUE)
# a 'NULL' object is always empty

###

# checking whether two objects represent the same mathematical entity
if (!isGeneric("areTheSame")) {
    setGeneric("areTheSame",function(x,y) standardGeneric("areTheSame"))
} # end of if

###

# checking whether an object is present in another object
if (!isGeneric("isPresent")) {
    setGeneric("isPresent",function(el,ou) standardGeneric("isPresent"))
} # end of if

###

# getting the maximum numeric identifier of an object
if (!isGeneric("maxId")) {
    setGeneric("maxId",function(x) standardGeneric("maxId"))
} # end of if

###

# recoding an object from a source code to a destination code
if (!isGeneric("recode")) {
    setGeneric("recode",function(object,src,dst) standardGeneric("recode"))
} # end of if

###

# getting the incidence list
if (!isGeneric("incidenceList")) {
    setGeneric("incidenceList",function(object, ...) standardGeneric("incidenceList"))
  } # end of if

# getting the incidence matrix
if (!isGeneric("incidenceMatrix")) {
    setGeneric("incidenceMatrix",function(object, ...) standardGeneric("incidenceMatrix"))
  } # end of if

# getting the adjacency list
if (!isGeneric("adjacencyList")) {
    setGeneric("adjacencyList",function(object, ...) standardGeneric("adjacencyList"))
  } # end of if

# getting the adjacency matrix
if (!isGeneric("adjacencyMatrix")) {
    setGeneric("adjacencyMatrix",function(object, ...) standardGeneric("adjacencyMatrix"))
  } # end of if

###

# setting the incidence list
if (!isGeneric("incidenceList<-")) {
    setGeneric("incidenceList<-",function(x, force=TRUE, value) standardGeneric("incidenceList<-"))
  } # end of if

# setting the incidence matrix
if (!isGeneric("incidenceMatrix<-")) {
    setGeneric("incidenceMatrix<-",function(x, force=TRUE, value) standardGeneric("incidenceMatrix<-"))
  } # end of if

# setting the adjacency list
if (!isGeneric("adjacencyList<-")) {
    setGeneric("adjacencyList<-",function(x, force=TRUE, value) standardGeneric("adjacencyList<-"))
  } # end of if

# setting the adjacency matrix
if (!isGeneric("adjacencyMatrix<-")) {
    setGeneric("adjacencyMatrix<-",function(x, force=TRUE, value) standardGeneric("adjacencyMatrix<-"))
  } # end of if

###

# 'dynamic.Graph'
if (!isGeneric("dynamic.Graph")) {
  setGeneric("dynamic.Graph",function(object, ...) standardGeneric("dynamic.Graph"))
  } # end of if

# 'display'
if (!isGeneric("display")) {
    setGeneric("display",function(x,...) standardGeneric("display"))
  } # end of if
