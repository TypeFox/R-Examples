#  Copyright (C) 2012 Heather Turner
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

Nonlin <- function(functionCall){
    .Defunct(msg = paste("'Nonlin' is defunct.",
             "\nUse functions of class \"nonlin\" instead.",
             "\nSee ?nonlin.function for more details."))
}
class(Nonlin) <- "nonlin"

getModelFrame <- function() {
    .Defunct(msg = paste("'getModelFrame' is deprecated as it was designed to ",
             "work with the old plug-in architecture for nonlinear terms."))
}

qrSolve <- function(A, b, rank = NULL, ...) {
    .Defunct(msg = paste("'qrSolve' is deprecated as it is no longer used ",
             "by gnm."))
}
