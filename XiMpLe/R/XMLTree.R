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


#' Constructor function for XiMpLe.doc objects
#'
#' Can be used to create full XML trees.
#'
#' @param ... Optional children for the XML tree. Must be either objects of class
#'    \code{\link[XiMpLe:XiMpLe.node-class]{XiMpLe.node}} or character strings,
#'    which are treated as simple text values.
#' @param xml A named list, XML declaration of the XML tree. Currently just pasted, no checking is done.
#' @param dtd A named list, doctype definition of the XML tree. Valid elements are \code{doctype} (root element), \code{decl}
#' ("PUBLIC" or "SYSTEM"), \code{id} (the identifier) and \code{refer} (URI to .dtd).
#'    Currently just pasted, no checking is done.
#' @param .children Alternative way of specifying children, if you have them already as a list.
#' @return An object of class \code{\link[XiMpLe:XiMpLe.doc-class]{XiMpLe.doc}}
#' @seealso
#'    \code{\link[XiMpLe:XMLNode]{XMLNode}},
#'    \code{\link[XiMpLe:pasteXML]{pasteXML}}
#' @export
#' @rdname XMLTree
#' @examples
#' sample.XML.a <- XMLNode("a",
#'   attrs=list(href="http://example.com", target="_blank"),
#'   .children="klick here!")
#' sample.XML.body <- XMLNode("body", .children=list(sample.XML.a))
#' sample.XML.html <- XMLNode("html", .children=list(XMLNode("head", ""),
#'   sample.XML.body))
#' sample.XML.tree <- XMLTree(sample.XML.html,
#'   xml=list(version="1.0", encoding="UTF-8"),
#'   dtd=list(doctype="html", decl="PUBLIC",
#'     id="-//W3C//DTD XHTML 1.0 Transitional//EN",
#'     refer="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"))

XMLTree <- function(..., xml=NULL, dtd=NULL, .children=list(...)){

  # remove NULLs
  .children <- .children[unlist(lapply(.children, length) != 0)]

  # check for text values
  all.children <- sapply(child.list(.children), function(this.child){
    if(is.character(this.child)){
      this.child <- new("XiMpLe.node",
          name="",
          value=this.child
        )
    } else {}
    return(this.child)
  })

  if(is.null(xml)){
    xml <- list()
  } else {}
  if(is.null(dtd)){
    dtd <- list()
  } else {}
  
  newTree <- new("XiMpLe.doc",
    xml=xml,
    dtd=dtd,
    children=all.children
  )

  return(newTree)
}