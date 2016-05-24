{###############################################################################
# T.R
# This file is part of the R package dostats
# 
# Copyright 2012 Andrew Redd
# Date: 5/30/2012
# 
# DESCRIPTION
# ===========
# Helper function for specifying text without quotations.
# 
# LICENSE
# ========
# dostats is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
# 
# dostats is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with 
# dostats. If not, see http://www.gnu.org/licenses/.
# 
}###############################################################################

#' create a text vector
#' @rdname T
#' @param ... names, quoted or not, no substitution made
#' @export
#' @examples
#' .T(min, mean, 'median')
.T <- function(...){
    c <- as.character(x <- substitute(c(...)))[-1]
    names(c) <- names(as.list(x))[-1]
    return(c)
}



