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


# Class XiMpLe.node
#
# This class is used to create DOM trees of XML documents, like objects that are returned
# by \code{\link[XiMpLe:parseXMLTree]{parseXMLTree}}.
# 
# There are certain special values predefined for the \code{name} slot to easily create special XML elements:
# \describe{
#     \item{\code{name=""}}{If the name is an empty character string, a pseudo node is created,
#       \code{\link[XiMpLe:pasteXMLNode]{pasteXMLNode}} will paste its \code{value} as plain text.}
#     \item{\code{name="!--"}}{Creates a comment tag, i.e., this will comment out all its \code{children}.}
#     \item{\code{name="![CDATA["}}{Creates a CDATA section and places all its \code{children} in it.}
#     \item{\code{name="*![CDATA["}}{Creates a CDATA section and places all its \code{children} in it, where the CDATA markers are
#       commented out by \code{/* */}, as is used for JavaScript in XHTML.}
# }
#
# @slot name Name of the node (i.e., the XML tag identifier). For special names see details.
# @slot attributes A list of named character values, representing the attributes of this node.
# @slot children A list of further objects of class XiMpLe.node, representing child nodes of this node.
# @slot value Plain text to be used as the enclosed value of this node. Set to \code{value=""} if you
#    want a childless node to be forced into an non-empty pair of start and end tags by \code{\link[XiMpLe:pasteXMLNode]{pasteXMLNode}}.
# @name XiMpLe.node,-class
# @aliases XiMpLe.node-class XiMpLe.node,-class
#' @import methods
# @keywords classes
# @rdname XiMpLe.node-class
#' @export

setClass("XiMpLe.node",
  representation=representation(
    name="character",
    attributes="list",
    children="list",
    value="character"
  ),
  prototype(
    name=character(),
    attributes=list(),
    children=list(),
    value=character()
  )
)

setValidity("XiMpLe.node", function(object){
    obj.name <- object@name
    obj.attributes <- object@attributes
    obj.children <- object@children
    obj.value <- object@value

    if(isTRUE(!nchar(obj.name) > 0) & isTRUE(!nchar(obj.value) > 0)){
      print(str(object))
      stop(simpleError("Invalid object: A node must at least have a name or a value!"))
    } else {}

    obj.attributes.names <- names(obj.attributes)
    # if there are attributes, check that they all have names
    if(length(obj.attributes) > 0){
      if(length(obj.attributes) != length(obj.attributes.names)){
        stop(simpleError("Invalid object: All attributes must have names!"))
      } else {}
    } else {}

    # check content of children
    if(length(obj.children) > 0){
      child.nodes <- sapply(obj.children, function(this.child){is.XiMpLe.node(this.child)})
      if(!all(child.nodes)){
        stop(simpleError("Invalid object: All list elements of children must be of class XiMpLe.node!"))
      } else {}
    } else {}
  return(TRUE)
})
