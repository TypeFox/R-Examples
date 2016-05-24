#  sybilStack.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybil.
#
#  Sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
# initialize stack

stinit <- function(stname) {

    .SYBILenv$sybilStack[[stname]] <- vector(mode = "list", length = 0)

}


#------------------------------------------------------------------------------#
# remove stack

stclear <- function(stname) {

#    if (!(stname %in% names(.SYBILenv$sybilStack))) {
#        stop("stack ", sQuote(stname), " does not exist")
#    }
    stopifnot(stexists(stname))
    .SYBILenv$sybilStack[stname] <- NULL

}


#------------------------------------------------------------------------------#
# add element at last position

stpush <- function(stname, value) {

    stopifnot(stexists(stname))

    .SYBILenv$sybilStack[[stname]] <- append(.SYBILenv$sybilStack[[stname]],
                                 value,
                                 after = length(.SYBILenv$sybilStack[[stname]]))
}


#------------------------------------------------------------------------------#
# remove last element and return invisibly

stpop <- function(stname) {

    stopifnot(stexists(stname))

    lg <- length(.SYBILenv$sybilStack[[stname]])
    value <- .SYBILenv$sybilStack[[stname]][lg]
    length(.SYBILenv$sybilStack[[stname]]) <- lg - 1
    return(invisible(value))

}


#------------------------------------------------------------------------------#
# add element at first position

stunshift <- function(stname, value) {

    stopifnot(stexists(stname))

    .SYBILenv$sybilStack[[stname]] <- append(.SYBILenv$sybilStack[[stname]],
                                             value, after = 0)

}

#------------------------------------------------------------------------------#
# remove first element and return invisibly

stshift <- function(stname) {

    stopifnot(stexists(stname))

    value <- .SYBILenv$sybilStack[[stname]][1]
    .SYBILenv$sybilStack[[stname]] <- .SYBILenv$sybilStack[[stname]][-1]
    return(invisible(value))

}


#------------------------------------------------------------------------------#
# return last element, not removed

stseek <- function(stname) {

    stopifnot(stexists(stname))

    lg <- length(.SYBILenv$sybilStack[[stname]])
    value <- .SYBILenv$sybilStack[[stname]][lg]
    return(value)

}


#------------------------------------------------------------------------------#
# return first element, not removed

stfirst <- function(stname) {

    stopifnot(stexists(stname))

    value <- .SYBILenv$sybilStack[[stname]][1]
    return(value)

}


#------------------------------------------------------------------------------#
# return complete stack as list

stlist <- function(stname) {

    stopifnot(stexists(stname))

    return(.SYBILenv$sybilStack[[stname]])

}


#------------------------------------------------------------------------------#
# return number of elements in stack

stlength <- function(stname) {

    stopifnot(stexists(stname))

    return(length(.SYBILenv$sybilStack[[stname]]))

}


#------------------------------------------------------------------------------#
# return TRUE is stack stname exists, FALSE otherwise

stexists <- function(stname) {

    return(stname %in% names(.SYBILenv$sybilStack))

}

