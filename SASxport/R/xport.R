### This file is part of the 'foreign' package for R.

###
###             Read SAS xport format libraries
###
### Copyright 1999-1999 Douglas M. Bates <bates$stat.wisc.edu>,
###                     Saikat DebRoy <saikat$stat.wisc.edu>
###
### This file is part of the `foreign' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, a copy is available at
### http://www.r-project.org/Licenses/

lookup.xport.inner <- function(file) {
  .Call('xport_info', file, PACKAGE = "SASxport")
}

read.xport.inner <- function(file, stringsAsFactors=FALSE) {
    data.info <- lookup.xport.inner(file)
    ans <- .Call('xport_read', file, data.info, PACKAGE = "SASxport")
    if (length(ans) == 1L)
      as.data.frame(ans[[1L]], stringsAsFactors=stringsAsFactors)
    else
      lapply(ans, as.data.frame, stringsAsFactors=stringsAsFactors)
}
