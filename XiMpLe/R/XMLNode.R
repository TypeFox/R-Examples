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


#' Constructor function for XiMpLe.node objects
#' 
#' Can be used to create XML nodes.
#' 
#' To generate a CDATA node, set \code{name="![CDATA["}, to create a comment, set \code{name="!--"}.
#' 
#' @param name Character string, the tag name.
#' @param ... Optional children for the tag. Must be either objects of class XiMpLe.node or character strings,
#'    which are treated as simple text values. If this is empty, the tag will be treated as an empty tag. To
#'    force a closing tag, supply an empty string, i.e. \code{""}.
#' @param attrs An optional named list of attributes.
#' @param namespace Currently ignored.
#' @param namespaceDefinitions Currently ignored.
#' @param .children Alternative way of specifying children, if you have them already as a list.
#' @return An object of class \code{\link[XiMpLe:XiMpLe.node-class]{XiMpLe.node}}.
#' @seealso
#'    \code{\link[XiMpLe:XMLTree]{XMLTree}},
#'    \code{\link[XiMpLe:pasteXML]{pasteXML}}
#' @export
#' @rdname XMLNode
#' @examples
#' sample.XML.node <- XMLNode("a",
#'   attrs=list(href="http://example.com", target="_blank"),
#'   .children="klick here!")

XMLNode <- function(name, ..., attrs=NULL, namespace="", namespaceDefinitions=NULL, .children=list(...)){

  all.children <- list()

  # text node?
  if(identical(name, "") &
      (all(unlist(lapply(.children, is.character)))) |
      all(unlist(lapply(.children, is.numeric)))){
    value <- paste(..., sep=" ")
  } else if(identical(.children, list(""))){
    value <- ""
  } else {
    # remove NULLs
    .children <- .children[unlist(lapply(.children, length) != 0)]
    # check for text values
    all.children <- sapply(child.list(.children), function(this.child){
      if(is.character(this.child) | is.numeric(this.child)){
        this.child <- new("XiMpLe.node",
            name="",
            value=as.character(this.child)
          )
      } else {}
      return(this.child)
    })
    value <- character()
  }

  newNode <- new("XiMpLe.node",
    name=name,
    attributes=as.list(attrs),
    children=all.children,
    value=value)

  return(newNode)
}
