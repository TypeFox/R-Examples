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


#' Extract/manipulate a node or parts of it from an XML tree
#'
#' This method can be used to get parts of a parsed XML tree object, or to fill it with new values.
#'
#' @include 00_class_01_XiMpLe.node.R
#' @include 00_class_02_XiMpLe.doc.R
#' @include XiMpLe-internal.R
#' @import methods
#' @examples
#' \dontrun{
#' node(my.xml.tree, node=list("html","body"), what="attributes")
#' node(my.xml.tree, node=list("html","head","title"), what="value") <- "foobar"
#' }
#' @docType methods
#' @rdname node
#' @export
setGeneric("node", function(obj, node=list(), what=NULL, cond.attr=NULL, cond.value=NULL, element=NULL){standardGeneric("node")})

# define class union to make life easier
#' @export
setClassUnion("XiMpLe.XML", members=c("XiMpLe.node", "XiMpLe.doc"))

#' @rdname node
#' @aliases
#'    node,-methods
#'    node,XiMpLe.doc-method
#'    node,XiMpLe.node-method
#'    node,XiMpLe.XML-method
#' @param obj An object of class \code{\link[XiMpLe:XiMpLe.doc-class]{XiMpLe.doc}} or
#'    \code{\link[XiMpLe:XiMpLe.node-class]{XiMpLe.node}}.
#' @param node A list of node names (or their numeric values), where each element is
#'    the child of its previous element. duplicate matches will be returned as a list.
#' @param what A character string, must be a valid slot name of class \code{\link[XiMpLe:XiMpLe.node-class]{XiMpLe.node}}, like
#'    \code{"attributes"} or \code{"value"}. If not \code{NULL}, only that part of a node will be returned.
#'    There's also two special properties for this option: \code{what="@@path"} will not return the
#'    node or it's contents, but a character string with the "path" to it in the object; \code{what="obj@@path"}
#'    is the same but won't have \code{obj} substituted with the object's name.
#' @param cond.attr A named character string, to further filter the returned results.
#'    If not \code{NULL}, only nodes with fully matching attributes will be considered.
#' @param cond.value A character string, similar to \code{cond.attr}, but is matched
#'    against the value between a pair of tags.
#' @param element A character string naming one list element of the node slot. If \code{NULL}, all
#'    elements will be returned.
setMethod("node",
  signature(obj="XiMpLe.XML"),
  function(obj, node=list(), what=NULL, cond.attr=NULL, cond.value=NULL, element=NULL){

    # check top level if this is a node, not a tree
    if(is.XiMpLe.node(obj)){
      got.this <- identical(slot(obj, "name"), node[[1]])
      if(!isTRUE(got.this)){
        # apparently, this node doesn't exist
        stop(simpleError(paste0("Can't find node ", node[[1]], " in ", sQuote(deparse(substitute(obj))), "!")))
      } else {
        # remove first element in list node
        node[[1]] <- NULL
      }
    } else {}
    result.node.path <- "obj"
    for (this.node in node){
      for (this.path in result.node.path){
        this.node.part <- eval(parse(text=this.path))
        got.this <- lapply(slot(this.node.part, "children"), function(this.child){slot(this.child, "name")}) %in% this.node
        if(!any(got.this)){
          # apparently, this node doesn't exist
          stop(simpleError(paste0("Can't find node ", sQuote(this.node), " in ", sQuote(deparse(substitute(obj))), "!")))
        } else {
          result.node.path <- unique(paste0(result.node.path, paste0("@children[[",which(got.this),"]]")))
        }
      }
    }

    # filter by attributes
    if(!is.null(cond.attr)){
      filter <- names(cond.attr)
      filtered.paths <- c()
      for (this.path in result.node.path){
        this.node.part <- eval(parse(text=this.path))
        this.attrs <- slot(this.node.part, "attributes")
        attr.found <- filter %in% names(this.attrs)
        # were the wanted attributes found at all?
        if(all(attr.found)){
          # check if they're all equal
          found.this <- sapply(filter, function(this.attr){
              results <- unlist(cond.attr[this.attr]) == unlist(this.attrs[this.attr])
              return(results)
            })
          if(all(found.this)){
            filtered.paths <- unique(c(filtered.paths, this.path))
          } else {}
        } else {}
      }
      result.node.path <- filtered.paths
    } else {}
    # create a list with matching node objects
    result.cond <- sapply(result.node.path, function(this.path){eval(parse(text=this.path))})
    names(result.cond) <- NULL

    if(!is.null(cond.value)){
      stopifnot(length(cond.value) == 1)
      filtered.paths <- c()
      for (this.path in result.node.path){
        this.node.part <- eval(parse(text=this.path))
        this.value <- slot(this.node.part, "value")
        if(identical(this.value, cond.value)){
          filtered.paths <- unique(c(filtered.paths, this.path))
        } else {}
      }
      result.node.path <- filtered.paths
    } else {}

    if(!is.null(what)){
      stopifnot(length(what) == 1)
      if(!what %in% c(slotNames(new("XiMpLe.node")), "@path", "obj@path")){
        stop(simpleError(paste0("Invalid slot for class XiMpLe.node:", paste(sQuote(what), collapse=", "), "!")))
      } else {}
      if(identical(what, "@path")){
        ## return subtituted path info
        result.node.path <- gsub("^obj", paste(deparse(substitute(obj))), result.node.path)
        return(result.node.path)
      } else if(identical(what, "obj@path")){
        ## return path info
        return(result.node.path)
      } else {}
      result <- unlist(lapply(result.node.path, function(this.path){
          this.node <- eval(parse(text=this.path))
          results <- slot(this.node, what)
          # special case: text values can either be directly in the value slot of a node,
          # or in a pseudo tag as a child node, so we check both
          if(identical(what, "value")){
            for (this.child in slot(this.node, "children")){
                if(identical(slot(this.child, "name"), "") & isTRUE(nchar(slot(this.child, "value")) > 0))
                  results <- paste(slot(this.child, "value"), results, sep=" ")
              }
          } else {}
          if(!is.null(element)){
            results <- results[element]
          } else {}
          return(results)
        }))
      # turn from vector to list for attributes, so they can be reached with "$"
      if(identical(what, "attributes")){
        result <- as.list(result)
      } else {}
    } else {
      result <- unlist(lapply(result.node.path, function(this.path){
          return(eval(parse(text=this.path)))
        }))
    }

    # no need for a list if it's inly one node
    if(length(result) == 1){
      if(is.XiMpLe.node(result[[1]]) | !is.null(element)){
        result <- result[[1]]
      } else {}
    } else {}

    return(result)
  }
)

#' @param value The value to set.
#' @export
#' @rdname node
setGeneric("node<-", function(obj, node=list(), what=NULL, cond.attr=NULL, cond.value=NULL, element=NULL, value){standardGeneric("node<-")})

#' @rdname node
#' @aliases
#'    node<-,-methods
#'    node<-,XiMpLe.doc-method
#'    node<-,XiMpLe.node-method
#'    node<-,XiMpLe.XML-method
setMethod("node<-",
  signature(obj="XiMpLe.XML"),
  function(obj, node=list(), what=NULL, cond.attr=NULL, cond.value=NULL, element=NULL, value){

  # get path to node in object
  obj.paths <- node(obj, node=node, what="obj@path", cond.attr=cond.attr, cond.value=cond.value)
  if(is.null(obj.paths)){
    # seems we need to create this node
    stop(simpleError("Node not found."))
  } else {}
  for (this.node in obj.paths){
    if(!is.null(what)){
      # special case: text values can either be directly in the value slot of a node,
      # or in a pseudo tag as a child node, so we check both and remove all
      if(identical(what, "value")){
        eval(parse(text=paste0(this.node, "@value <- character()")))
        all.node.children <- slot(eval(parse(text=this.node)), "children")
        child.is.value <- sapply(all.node.children, function(this.child){
            if(identical(slot(this.child, "name"), "") & isTRUE(nchar(slot(this.child, "value")) > 0)){
              return(TRUE)
            } else {
              return(FALSE)
            }
          })
        # if we have a mix of pseudo and actual tags, we probably messed up the markup
        if(length(all.node.children) != length(child.is.value)){
          warning("a child node contained text values and other nodes, we probably messed up the markup!")
        } else {}
        remove.nodes <- paste0(this.node, "@children[child.is.value] <- NULL")
        eval(parse(text=remove.nodes))

        # paste new value into a single pseudo node
        pseudo.node <- paste0(this.node, "@children <- append(", this.node, "@children, ",
          "new(\"XiMpLe.node\", name=\"\", value=\"", value, "\"), after=0)")
        eval(parse(text=pseudo.node))

        # now return the object
        return(obj)
      } else {
        this.node <- paste0(this.node, "@", what)
      }

      if(!is.null(element)){
        this.node <- paste0(this.node, "[[\"",element,"\"]]")
      } else {}
    } else {}

    eval(parse(text=paste0(this.node, " <- ", deparse(value))))
  }
  
    return(obj)
  }
)
