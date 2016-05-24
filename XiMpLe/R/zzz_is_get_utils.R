# Copyright 2011-2014 Meik Michalke <meik.michalke@hhu.de>
#
# This file is part of the R package XiMpLe.
#
# XiMpLe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# XiMpLe is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with XiMpLe.  If not, see <http://www.gnu.org/licenses/>.


## the name "zzz_*" is just to ensure roxygen doesn't parse it before XMLNode.R and XMLTree.R

#' @param x An arbitrary \code{R} object.
#' @rdname XMLNode
#' @export
is.XiMpLe.node <- function(x){
  inherits(x, "XiMpLe.node")
}

#' @param x An arbitrary \code{R} object.
#' @rdname XMLTree
#' @export
is.XiMpLe.doc <- function(x){
  inherits(x, "XiMpLe.doc")
}


#' Getter/setter methods for S4 objects of XiMpLe XML classes
#'
#' Used to get/set certain slots from objects of class \code{\link[XiMpLe:XiMpLe.doc-class]{XiMpLe.doc}}
#' and \code{\link[XiMpLe:XiMpLe.node-class]{XiMpLe.node}}.
#' 
#' These are convenience methods to get or set slots from XML objects without using the \code{@@} operator.
#'
#' \itemize{
#'    \item{\code{XMLName()}: }{get/set the XML node name (slot \code{name} of class \code{XiMpLe.node})}
#'    \item{\code{XMLAttrs()}: }{get/set the XML node attributes (slot \code{attrs} of class \code{XiMpLe.node})}
#'    \item{\code{XMLValue()}: }{get/set the XML node value (slot \code{value} of class \code{XiMpLe.node})}
#'    \item{\code{XMLChildren()}: }{get/set the XML child nodes (slot \code{children} of both classes \code{XiMpLe.node}
#'      and  \code{XiMpLe.doc})}
#'    \item{\code{XMLFile()}: }{get/set the XML document file name  (slot \code{file} of class \code{XiMpLe.doc})}
#'    \item{\code{XMLDecl()}: }{get/set the XML document declaration (slot \code{xml} of class \code{XiMpLe.doc})}
#'    \item{\code{XMLDTD()}: }{get/set the XML document doctype definition (slot \code{dtd} of class \code{XiMpLe.doc})}
#' }
#'
#' Another special method can scan a node/document tree object for appearances of nodes with a particular name:
#'
#' \itemize{
#'    \item{\code{XMLScan(obj, name, as.list=FALSE)}: }{get/set the XML nodes by name (recursively searches slot \code{name} of both classes
#'      \code{XiMpLe.node} and  \code{XiMpLe.doc}). If \code{as.list=TRUE} allways returns a list (or NULL), otherwise if exactly one result is found,
#'      it will be returned as as single \code{XiMpLe.node}.}
#' }
#'
#' Finally, there is a method to scan for certain values in XiMpLe objects and just list them. For instance, it can be used to
#' list all instances of a certain attribute type in a document tree:
#'
#' \itemize{
#'    \item{\code{XMLScanDeep(obj, find, search="attributes")}: }{returns all found instances of \code{find} in all slots defined by \code{search}.}
#' }
#' @param obj An object of class \code{XiMpLe.node} or \code{XiMpLe.doc}
#' @seealso
#'    \code{\link[XiMpLe:node]{node}},
#'    \code{\link[XiMpLe:XiMpLe.doc-class]{XiMpLe.doc}},
#'    \code{\link[XiMpLe:XiMpLe.node-class]{XiMpLe.node}}
#' @keywords methods
#' @docType methods
#' @rdname XMLGetters-methods
#' @export
#' @examples
#' xmlTestNode <- XMLNode("foo", XMLNode("testchild"))
#' XMLName(xmlTestNode) # returns "foo"
#' XMLName(xmlTestNode) <- "bar"
#' XMLName(xmlTestNode) # now returns "bar"
#'
#' # search for a child node
#' XMLScan(xmlTestNode, "testchild")
#' # remove nodes of that name
#' XMLScan(xmlTestNode, "testchild") <- NULL
setGeneric("XMLName", function(obj){standardGeneric("XMLName")})

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLName,-methods
#'    XMLName,XiMpLe.node-method
#' @docType methods
#' @include 00_class_01_XiMpLe.node.R
setMethod("XMLName",
  signature=signature(obj="XiMpLe.node"),
  function(obj){
    return(obj@name)
  }
)

#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLName<-", function(obj, value){standardGeneric("XMLName<-")})

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLName<-,-methods
#'    XMLName<-,XiMpLe.node-method
#' @docType methods
#' @include 00_class_01_XiMpLe.node.R
setMethod("XMLName<-",
  signature=signature(obj="XiMpLe.node"),
  function(obj, value){
    obj@name <- value
    return(obj)
  }
)

#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLAttrs", function(obj){standardGeneric("XMLAttrs")})

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLAttrs,-methods
#'    XMLAttrs,XiMpLe.node-method
#' @docType methods
#' @include 00_class_01_XiMpLe.node.R
setMethod("XMLAttrs",
  signature=signature(obj="XiMpLe.node"),
  function(obj){
    return(obj@attributes)
  }
)

#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLAttrs<-", function(obj, value){standardGeneric("XMLAttrs<-")})

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLAttrs<-,-methods
#'    XMLAttrs<-,XiMpLe.node-method
#' @docType methods
#' @include 00_class_01_XiMpLe.node.R
setMethod("XMLAttrs<-",
  signature=signature(obj="XiMpLe.node"),
  function(obj, value){
    obj@attributes <- value
    return(obj)
  }
)

#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLChildren", function(obj){standardGeneric("XMLChildren")})

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLChildren,-methods
#'    XMLChildren,XiMpLe.node-method
#' @docType methods
#' @include 00_class_01_XiMpLe.node.R
setMethod("XMLChildren",
  signature=signature(obj="XiMpLe.node"),
  function(obj){
    return(obj@children)
  }
)

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLChildren,-methods
#'    XMLChildren,XiMpLe.doc-method
#' @docType methods
#' @include 00_class_02_XiMpLe.doc.R
setMethod("XMLChildren",
  signature=signature(obj="XiMpLe.doc"),
  function(obj){
    return(obj@children)
  }
)

#' @param value The new value to set.
#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLChildren<-", function(obj, value){standardGeneric("XMLChildren<-")})

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLChildren<-,-methods
#'    XMLChildren<-,XiMpLe.node-method
#' @docType methods
#' @include 00_class_01_XiMpLe.node.R
setMethod("XMLChildren<-",
  signature=signature(obj="XiMpLe.node"),
  function(obj, value){
    obj@children <- child.list(value)
    return(obj)
  }
)

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLChildren<-,-methods
#'    XMLChildren<-,XiMpLe.doc-method
#' @docType methods
#' @include 00_class_02_XiMpLe.doc.R
setMethod("XMLChildren<-",
  signature=signature(obj="XiMpLe.doc"),
  function(obj, value){
    obj@children <- child.list(value)
    return(obj)
  }
)


#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLValue", function(obj){standardGeneric("XMLValue")})

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLValue,-methods
#'    XMLValue,XiMpLe.node-method
#' @docType methods
#' @include 00_class_01_XiMpLe.node.R
setMethod("XMLValue",
  signature=signature(obj="XiMpLe.node"),
  function(obj){
    return(obj@value)
  }
)

#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLValue<-", function(obj, value){standardGeneric("XMLValue<-")})

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLValue<-,-methods
#'    XMLValue<-,XiMpLe.node-method
#' @docType methods
#' @include 00_class_01_XiMpLe.node.R
setMethod("XMLValue<-",
  signature=signature(obj="XiMpLe.node"),
  function(obj, value){
    obj@value <- value
    return(obj)
  }
)

#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLFile", function(obj){standardGeneric("XMLFile")})

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLFile,-methods
#'    XMLFile,XiMpLe.doc-method
#' @docType methods
#' @include 00_class_02_XiMpLe.doc.R
setMethod("XMLFile",
  signature=signature(obj="XiMpLe.doc"),
  function(obj){
    return(obj@file)
  }
)

#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLFile<-", function(obj, value){standardGeneric("XMLFile<-")})

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLFile<-,-methods
#'    XMLFile<-,XiMpLe.doc-method
#' @docType methods
#' @include 00_class_02_XiMpLe.doc.R
setMethod("XMLFile<-",
  signature=signature(obj="XiMpLe.doc"),
  function(obj, value){
    obj@file <- value
    return(obj)
  }
)

#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLDecl", function(obj){standardGeneric("XMLDecl")})

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLDecl,-methods
#'    XMLDecl,XiMpLe.doc-method
#' @docType methods
#' @include 00_class_02_XiMpLe.doc.R
setMethod("XMLDecl",
  signature=signature(obj="XiMpLe.doc"),
  function(obj){
    return(obj@xml)
  }
)

#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLDecl<-", function(obj, value){standardGeneric("XMLDecl<-")})

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLDecl<-,-methods
#'    XMLDecl<-,XiMpLe.doc-method
#' @docType methods
#' @include 00_class_02_XiMpLe.doc.R
setMethod("XMLDecl<-",
  signature=signature(obj="XiMpLe.doc"),
  function(obj, value){
    obj@xml <- value
    return(obj)
  }
)

#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLDTD", function(obj){standardGeneric("XMLDTD")})

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLDTD,-methods
#'    XMLDTD,XiMpLe.doc-method
#' @docType methods
#' @include 00_class_02_XiMpLe.doc.R
setMethod("XMLDTD",
  signature=signature(obj="XiMpLe.doc"),
  function(obj){
    return(obj@dtd)
  }
)

#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLDTD<-", function(obj, value){standardGeneric("XMLDTD<-")})

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLDTD<-,-methods
#'    XMLDTD<-,XiMpLe.doc-method
#' @docType methods
#' @include 00_class_02_XiMpLe.doc.R
setMethod("XMLDTD<-",
  signature=signature(obj="XiMpLe.doc"),
  function(obj, value){
    obj@dtd <- value
    return(obj)
  }
)

## scan a tree for appearances of nodes
#' @param name Character, name of nodes to scan for.
#' @param as.list Logical, if \code{TRUE} allways returns a list (or NULL), otherwise if exactly one result is found,
#'    it will be returned as as single \code{XiMpLe.node}.
#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLScan", function(obj, name, as.list=FALSE){standardGeneric("XMLScan")})

# internal helper function
find.nodes <- function(nodes, nName){
  res <- list()
  for (thisNode in nodes){
      if(identical(XMLName(thisNode), nName)){
        res <- append(res, thisNode)
      } else if(length(XMLChildren(thisNode)) > 0){
        res <- append(res, find.nodes(XMLChildren(thisNode), nName=nName))
      } else {}
    }
  return(res)
}

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLScan,-methods
#'    XMLScan,XiMpLe.node-method
#' @docType methods
#' @include 00_class_01_XiMpLe.node.R
setMethod("XMLScan",
  signature=signature(obj="XiMpLe.node"),
  function(obj, name, as.list=FALSE){
    node.list <- find.nodes(
      nodes=child.list(obj),
      nName=name)
    if(identical(node.list, list())){
      return(NULL)
    } else if(length(node.list) == 1 && !isTRUE(as.list)){
      return(node.list[[1]])
    } else {
      return(node.list)
    }
  }
)

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLScan,XiMpLe.doc-method
#' @docType methods
#' @include 00_class_02_XiMpLe.doc.R
setMethod("XMLScan",
  signature=signature(obj="XiMpLe.doc"),
  function(obj, name, as.list=FALSE){
    node.list <- find.nodes(
      nodes=XMLChildren(obj),
      nName=name)
    if(identical(node.list, list())){
      return(NULL)
    } else if(length(node.list) == 1 && !isTRUE(as.list)){
      return(node.list[[1]])
    } else {
      return(node.list)
    }
  }
)

#' @rdname XMLGetters-methods
#' @exportMethod XMLScan<-
setGeneric("XMLScan<-", function(obj, name, value){standardGeneric("XMLScan<-")})

# internal helper function
replace.nodes <- function(nodes, nName, replacement){
  nodes <- sapply(nodes, function(thisNode){
      if(identical(XMLName(thisNode), nName)){
        return(replacement)
      } else if(length(XMLChildren(thisNode)) > 0){
        XMLChildren(thisNode) <- replace.nodes(
            XMLChildren(thisNode),
            nName=nName,
            replacement=replacement)
        return(thisNode)
      } else {
        return(thisNode)
      }
    })
  # get rid of NULL in list
  nodes <- Filter(Negate(is.null), nodes)
  return(nodes)
}

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLScan<-,-methods
#'    XMLScan<-,XiMpLe.node-method
#' @docType methods
#' @include 00_class_01_XiMpLe.node.R
setMethod("XMLScan<-",
  signature=signature(obj="XiMpLe.node"),
  function(obj, name, value){
    # prevent the creation of invalid results
    stopifnot(is.XiMpLe.node(value) || is.null(value))
    obj <- replace.nodes(
      nodes=child.list(obj),
      nName=name,
      replacement=value)
    stopifnot(validObject(object=obj, test=TRUE, complete=TRUE))
    if(identical(obj, as.list(value))){
      # it seems the full object was replaced by value
      return(value)
    } else {
      return(obj[[1]])
    }
  }
)

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLScan<-,XiMpLe.doc-method
#' @docType methods
#' @include 00_class_02_XiMpLe.doc.R
setMethod("XMLScan<-",
  signature=signature(obj="XiMpLe.doc"),
  function(obj, name, value){
    # prevent the creation of invalid results
    stopifnot(is.XiMpLe.node(value) || is.null(value))
    XMLChildren(obj) <- replace.nodes(
        nodes=XMLChildren(obj),
        nName=name,
        replacement=value)[[1]]
    stopifnot(validObject(object=obj, test=TRUE, complete=TRUE))
    return(obj)
  }
)

#' @param find Character, name of element to scan for.
#' @param search Character, name of the slot to scan, one of \code{"attributes"},
#'    \code{"name"}, or \code{"value"} for nodes.
#' @rdname XMLGetters-methods
#' @docType methods
#' @export
setGeneric("XMLScanDeep", function(obj, find=NULL, search="attributes"){standardGeneric("XMLScanDeep")})

# internal helper function
recursiveScan <- function(robj, rfind, rsearch, recResult=list(), result, envID="all"){
  if(is.XiMpLe.doc(robj)){
    recResult <- append(recResult, lapply(robj@children, function(this.child){
      recursiveScan(robj=this.child, rfind=rfind, rsearch=rsearch, recResult=recResult, result=result, envID=envID)
    }))
  } else if(is.XiMpLe.node(robj)){
    foundThis <- NULL
    if(identical(rsearch, "attributes")){
      foundThis <- XMLAttrs(robj)[[rfind]]
    } else if(identical(rsearch, "name")){
      objName <- XMLName(robj)
      if(identical(objName, rfind)){
        foundThis <- objName
      } else {}
    } else if(identical(rsearch, "value")){
      objValue <- XMLValue(robj)
      if(identical(objValue, rfind)){
        foundThis <- objValue
      } else {}
    } else {
      stop(simpleError("Only \"attributes\", \"name\" or \"value\" are valid for search!"))
    }
    if(!is.null(foundThis)){
      thisResult <- as.list(result)
      nodeName <- XMLName(robj)
      thisResult[[envID]] <- append(thisResult[[envID]], foundThis)
      names(thisResult[[envID]])[length(thisResult[[envID]])] <- nodeName
      list2env(thisResult, envir=result)
    } else {}
    recResult <- append(recResult, lapply(robj@children, function(this.child){
      recursiveScan(robj=this.child, rfind=rfind, rsearch=rsearch, recResult=recResult, result=result, envID=envID)
    }))
  }
  return(recResult)
}

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLScanDeep,-methods
#'    XMLScanDeep,XiMpLe.node-method
#' @docType methods
#' @include 00_class_01_XiMpLe.node.R
setMethod("XMLScanDeep",
  signature=signature(obj="XiMpLe.node"),
  function(obj, find, search){
    result <- new.env()
    assign(find, c(), envir=result)
    recursiveScan(robj=obj, rfind=find, rsearch=search, recResult=list(), result=result, envID=find)
    return(get(find, envir=result))
  }
)

#' @rdname XMLGetters-methods
#' @aliases
#'    XMLScanDeep,XiMpLe.doc-method
#' @docType methods
#' @include 00_class_02_XiMpLe.doc.R
setMethod("XMLScanDeep",
  signature=signature(obj="XiMpLe.doc"),
  function(obj, find, search){
    result <- new.env()
    assign(find, c(), envir=result)
    recursiveScan(robj=obj, rfind=find, rsearch=search, recResult=list(), result=result, envID=find)
    return(get(find, envir=result))
  }
)
