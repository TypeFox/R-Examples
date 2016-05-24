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


#' Paste methods for XiMpLe XML objects
#' 
#' These methods can be used to paste objects if class \code{\link[XiMpLe:XiMpLe.node-class]{XiMpLe.node}}
#' or \code{\link[XiMpLe:XiMpLe.doc-class]{XiMpLe.doc}}.
#'
#' @note The functions pasteXMLNode() and pasteXMLTree() have been replaced by the pasteXML methods.
#'    They should no longer be used.
#'
#' @param obj An object of class \code{XiMpLe.node} or \code{XiMpLe.doc}.
#' @param ... Additional options for the generic method, see options for a specific method, respectively.
#' @aliases
#'    pasteXML,-methods
#'    pasteXML,XiMpLe.doc-method
#'    pasteXMLNode
#'    pasteXMLTree
#' @seealso \code{\link[XiMpLe:XiMpLe.node-class]{XiMpLe.node}}, 
#'    \code{\link[XiMpLe:XiMpLe.doc-class]{XiMpLe.doc}}
#' @keywords methods
#' @import methods
#' @rdname pasteXML-methods
#' @include 00_class_01_XiMpLe.node.R
#' @include 00_class_02_XiMpLe.doc.R
#' @docType methods
#' @export
setGeneric("pasteXML", function(obj, ...){
  standardGeneric("pasteXML")
})

#' @param level Indentation level.
#' @param shine Integer, controlling if the output should be formatted for better readability. Possible values:
#'    \describe{
#'      \item{0}{No formatting.}
#'      \item{1}{Nodes will be indented.}
#'      \item{2}{Nodes will be indented and each attribute gets a new line.}
#'    }
#' @param indent.by A charachter string defining how indentation should be done. Defaults to tab.
#' @param tidy Logical, if \code{TRUE} the special characters "<" and ">" will be replaced with the entities
#'    "&lt;" and "gt;" in attributes and text values.
#' @rdname pasteXML-methods
#' @aliases
#'    pasteXML,XiMpLe.node-method
setMethod("pasteXML",
  signature=signature(obj="XiMpLe.node"),
  function(obj, level=1, shine=1, indent.by="\t", tidy=TRUE){

    new.indent <- ifelse(shine > 0, indent(level+1, by=indent.by), "")
    new.node   <- ifelse(shine > 0, "\n", "")

    # get the slot contents
    node.name <- slot(obj, "name")
    node.attr <- slot(obj, "attributes")
    node.chld <- slot(obj, "children")
    node.val  <- slot(obj, "value")

    if(!length(node.attr) > 0){
      node.attr <- NULL
    } else {}

    if(length(node.chld) > 0){
      node.chld <- paste0(unlist(sapply(node.chld, function(this.node){
        if(slot(this.node, "name") == ""){
          this.node.pasted <- paste0(new.indent, pasteXML(this.node, level=level, shine=shine, indent.by=indent.by, tidy=tidy))
        } else {
          this.node.pasted <- pasteXML(this.node, level=(level + 1), shine=shine, indent.by=indent.by, tidy=tidy)
        }
        return(this.node.pasted)})), collapse="")
      node.empty <- FALSE
    } else {
      node.chld <- NULL
      node.empty <- TRUE
    }

    # take care of text value
    if(length(node.val) > 0){
      node.empty <- FALSE
      if(nchar(node.val) > 0){
        if(isTRUE(tidy)){
          node.val <- sapply(node.val, xml.tidy)
        } else {}
        node.chld <- paste0(node.chld, paste(node.val, new.node, collapse=" "))
      } else {}
    } else {}

    pasted.node <- pasteXMLTag(node.name, attr=node.attr, child=node.chld, empty=node.empty, level=level, allow.empty=TRUE, rename=NULL, shine=shine, indent.by=indent.by, tidy=tidy)
    
    return(pasted.node)
  }
)

#' @rdname pasteXML-methods
setMethod("pasteXML",
  signature=signature(obj="XiMpLe.doc"),
  function(obj, shine=1, indent.by="\t", tidy=TRUE){

    filename <- slot(obj, "file")
    tree.xml <- slot(obj, "xml")
    tree.doctype <- slot(obj, "dtd")
    tree.nodes <- slot(obj, "children")

    if(any(nchar(unlist(tree.xml)) > 0)) {
      doc.xml <- pasteXMLTag("?xml", attr=tree.xml, child=NULL, empty=TRUE, level=1, allow.empty=FALSE, rename=NULL, shine=min(1, shine), indent.by=indent.by, tidy=tidy)
      doc.xml <- gsub("/>", "\\?>", doc.xml)
    } else {
      doc.xml <- ""
    }

    if(any(nchar(unlist(tree.doctype)) > 0)) {
      new.node   <- ifelse(shine > 0, "\n", "")
      doc.doctype <- paste("<!DOCTYPE", tree.doctype[["doctype"]], tree.doctype[["decl"]], sep=" ")
      for (elmt in c("id", "refer")){
        if(length(tree.doctype[[elmt]]) > 0) {
          if(nchar(tree.doctype[[elmt]]) > 0){
            doc.doctype <- paste0(doc.doctype, " \"",tree.doctype[[elmt]], "\"")
          } else {}
        } else {}
      }
      doc.doctype <- paste0(doc.doctype, ">", new.node)
    } else {
      doc.doctype <- ""
    }

    if(length(tree.nodes) > 0) {
      doc.nodes <- paste0(unlist(sapply(tree.nodes, function(this.node){
        return(pasteXML(this.node, level=1, shine=shine, indent.by=indent.by, tidy=tidy))})), collapse="")
    } else {
      doc.nodes <- ""
    }

    doc.all <- paste0(doc.xml, doc.doctype, doc.nodes, collapse="")

    return(doc.all)
  }
)

# for compatibility reasons, deploy wrapper functions
#' @export
pasteXMLNode <- function(node, level=1, shine=1, indent.by="\t", tidy=TRUE){
  pasteXML(node, level=level, shine=shine, indent.by=indent.by, tidy=tidy)
}

#' @export
pasteXMLTree <- function(obj, shine=1, indent.by="\t", tidy=TRUE){
  pasteXML(obj, shine=shine, indent.by=indent.by, tidy=tidy)
}
