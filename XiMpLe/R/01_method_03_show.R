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


#' Show method for S4 objects of XiMpLe XML classes
#' 
#' Used to display objects of class \code{\link[XiMpLe:XiMpLe.doc-class]{XiMpLe.doc}}
#' and \code{\link[XiMpLe:XiMpLe.node-class]{XiMpLe.node}}
#'
#' @param object An object of class \code{XiMpLe.doc} or \code{XiMpLe.node}
#' @aliases
#'    show,-methods
#'    show,XiMpLe.doc-method
#'    show,XiMpLe.node-method
#'    show,XiMpLe.XML-method
#' @seealso  
#'    \code{\link[XiMpLe:XiMpLe.doc-class]{XiMpLe.doc}}
#'    \code{\link[XiMpLe:XiMpLe.node-class]{XiMpLe.node}}
#' @keywords methods
#' @docType methods
#' @exportMethod show
#' @rdname show-methods
#' @include 00_class_01_XiMpLe.node.R
#' @include 00_class_02_XiMpLe.doc.R
setMethod("show", signature(object="XiMpLe.XML"), function(object){
  cat(pasteXML(object))
})
