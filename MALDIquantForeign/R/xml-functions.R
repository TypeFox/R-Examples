## Copyright 2013 Sebastian Gibb
## <mail@sebastiangibb.de>
##
## This file is part of MALDIquantForeign for R and related languages.
##
## MALDIquantForeign is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## MALDIquantForeign is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with MALDIquantForeign. If not, see <http://www.gnu.org/licenses/>

.writeXmlHeader <- function(file, encoding) {
  cat("<?xml version=\"1.0\" encoding=\"", encoding, "\"?>\n", file=file,
      sep="")
 }

.writeCloseXmlTag <- function(tag, intend=0, file) {
  cat(paste0(rep(" ", times=intend)), "</", tag, ">\n",
      file=file, sep="", append=TRUE)
}

.writeXmlTag <- function(tag, text=NULL, attrs=NULL, intend=0, close=TRUE,
                         file) {
  intend <- paste0(rep(" ", times=intend))

  if (length(attrs)) {
    attrs <- paste0(" ", names(attrs), "=\"", attrs, "\"")
  }
  if (length(text)) {
    text <- paste0(">", text)

    if (close) {
      text <- paste0(text, "</", tag, ">")
    }
  } else {
    text <- ">"

    if (close) {
      text <- "/>"
    }
  }

  cat(intend, "<", tag, attrs, text, "\n", file=file, sep="", append=TRUE)
}
