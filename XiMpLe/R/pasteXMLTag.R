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


#' Write an XML tag
#'
#' Creates a whole XML tag with attributes and, if it is a pair of start and end tags,
#' also one object as child. This can be used recursively to create whole XML tree structures
#' with this one function.
#' 
#' @note However, you will probably not want to use this function at all, as it is much more
#'   comfortable to create XML nodes or even nested trees with \code{\link[XiMpLe:XMLNode]{XMLNode}} and
#'   \code{\link[XiMpLe:XMLTree]{XMLTree}}, and then feed the result to \code{\link[XiMpLe:pasteXML]{pasteXML}}.
#'
#' @param tag Character string, name of the XML tag.
#' @param attr A list of attributes for the tag.
#' @param child If \code{empty=FALSE}, a character string to be pasted as a child node between start and end tag.
#' @param empty Logical, <true /> or <false></false>
#' @param level Indentation level.
#' @param allow.empty Logical, if \code{FALSE}, tags without attributes will not be returned.
#' @param rename An optional named list if the attributes in XML need to be renamed from their list names in \code{attr}.
#'    This list must in turn have a list element named after \code{tag}, containing named character elements, where the
#'    names represent the element names in \code{attr} and their values the names the XML attribute should get.
#' @param shine Integer, controlling if the output should be formatted for better readability. Possible values:
#'    \describe{
#'      \item{0}{No formatting.}
#'      \item{1}{Nodes will be indented.}
#'      \item{2}{Nodes will be indented and each attribute gets a new line.}
#'    }
#' @param indent.by A charachter string defining how indentation should be done. Defaults to tab.
#' @param tidy Logical, if \code{TRUE} the special characters "<", ">" and "&" will be replaced with the entities
#'    "&lt;", "&gt;" and "&amp;" in attribute values. For comment or CDATA tags, if the text includes newline characters
#'    they will also be indented.
#' @return A character string.
#' @seealso \code{\link[XiMpLe:XMLNode]{XMLNode}},
#'    \code{\link[XiMpLe:XMLTree]{XMLTree}},
#'    \code{\link[XiMpLe:pasteXML]{pasteXML}}
#' @export
#' @examples
#' sample.XML.tag <- pasteXMLTag("a",
#'   attr=list(href="http://example.com", target="_blank"),
#'   child="klick here!",
#'   empty=FALSE)

pasteXMLTag <- function(tag, attr=NULL, child=NULL, empty=TRUE, level=1, allow.empty=FALSE, rename=NULL, shine=2, indent.by="\t", tidy=TRUE){
  # what attributes do we have?
  all.attributes <- pasteXMLAttr(attr, tag=tag, level=level, rename=rename, shine=shine, indent.by=indent.by, tidy=tidy)
  # probaly don't produce empty tags
  if(!isTRUE(allow.empty) & is.null(all.attributes)){
    return("")
  } else {}

  new.node   <- ifelse(shine > 0, "\n", "")
  new.indent <- ifelse(shine > 0, indent(level, by=indent.by), "")
  new.attr   <- ifelse(shine > 1, "\n", "")
  new.attr.indent <- ifelse(shine > 1, indent(level, by=indent.by), "")
  attr.space <- ifelse(nchar(all.attributes) > 0, " ", "")
  new.cmmt   <- ifelse(shine > 0, "\n", " ")
  new.cmmt.indent <- ifelse(shine > 1, indent(level + 1, by=indent.by), "")
  comment.indent <- ifelse(shine > 0, indent(level + 1, by=indent.by), "")

  # three special cases: value pseudotags, comments and CDATA
  if(isTRUE(nchar(tag) == 0) | length(tag) == 0){
    if(isTRUE(tidy)){
      child <- trim(child)
      child <- gsub("\n", new.cmmt, trim(setMinIndent(child, level=2, indent.by=indent.by)))
    }
    full.tag <- paste0(child, " ")
  } else {
    switch(tag,
     "!--"={
       # clean up value if needed
       if(!is.null(child)){
         if(isTRUE(tidy)){
           child <- trim(child)
           child <- gsub("\n", new.cmmt, trim(setMinIndent(child, level=level, indent.by=indent.by)))
         }
       } else {}
       full.tag <- paste0(new.indent, "<!-- ", new.attr, new.cmmt.indent,
         child, " ", new.attr, new.attr.indent,
         "-->", new.node)
       },
     "![CDATA["={
       # clean up value if needed
       if(!is.null(child)){
          child <- trim(child)
          if(isTRUE(tidy)){
            child <- gsub("\n", new.cmmt, trim(setMinIndent(child, level=level, indent.by=indent.by)))
          }
        } else {}
        full.tag <- paste0(new.indent, "<![CDATA[ ", new.cmmt, comment.indent,
          child, " ", new.cmmt, new.indent,
          "]]>", new.node)},
      "*![CDATA["={
        # clean up value if needed
        if(!is.null(child)){
          child <- trim(child)
          if(isTRUE(tidy)){
            child <- gsub("\n", new.cmmt, trim(setMinIndent(child, level=level, indent.by=indent.by)))
          }
        } else {}
        # 
        full.tag <- paste0(new.indent, "/* <![CDATA[ */ ", new.cmmt, comment.indent,
          child, " ", new.cmmt, new.indent,
          "/* ]]> */", new.node)},
      # last but not least, the default value
      {
        # only put attributes in new lines if there's more than one
        new.attr <- ifelse((length(attr) > 1), new.attr, "")
        new.attr.indent <- ifelse((length(attr) > 1), new.attr.indent, "")
        new.cmmt.indent <- ifelse((length(attr) > 1), new.cmmt.indent, "")
        val.indent <- ifelse(shine > 0, indent(level + 1, by=indent.by), "")
        # empty decides whether this is a empty tag or a pair of start and end tags
        if(isTRUE(empty)){
          full.tag <- paste0(new.indent, "<", tag, attr.space, new.attr, new.cmmt.indent, all.attributes, new.attr, new.attr.indent, " />", new.node)
        } else {
          full.tag <- paste0(
            new.indent, "<", tag, attr.space, new.attr, new.cmmt.indent, all.attributes, new.attr, new.attr.indent, ">", new.node,
            if(!is.null(child)){paste0(val.indent, trim(child), new.node)},
            new.indent, "</", tag, ">", new.node)
        }
      }
    )
  }
  return(full.tag)
}
