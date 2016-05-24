# /usr/bin/r
#
# Copyright 2015-2015 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of madness.
#
# madness is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# madness is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with madness.  If not, see <http://www.gnu.org/licenses/>.
#
# Created: 2015.11.30
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @include AllClass.r
#' @include utils.r
NULL

# coerce#FOLDUP
# 2FIX: add documentation and export

# http://stackoverflow.com/a/7766762/164611
#' @title Coerce madness to something else
#'
#' @description 
#'
#' Coerce as something else
#' 
#' @name as
#' @rdname as
#' @template etc
NULL

setAs(from='madness', to='array', def=function(from) val(from) )
setAs(from='madness', to='matrix', def=function(from) as.matrix(val(from)) )
setAs(from='madness', to='numeric', def=function(from) as.numeric(val(from)) )
setAs(from='madness', to='logical', def=function(from) as.logical(val(from)) )
setAs(from='madness', to='integer', def=function(from) as.integer(val(from)) )
setAs(from='madness', to='complex', def=function(from) as.complex(val(from)) )

# what. a fucking mess.
# http://stackoverflow.com/a/11285863/164611

#' @rdname as
#' @aliases as.array
#' @exportMethod as.array
#' @param x a \code{madness} object
#' @param ... further arguments passed to or from other methods.
setGeneric('as.array', signature="x", function(x,...) standardGeneric('as.array'))
#' @rdname as
#' @aliases as.array,madness-method
setMethod('as.array', 'madness', function(x,...) { as(x,'array') })

#' @rdname as
#' @aliases as.matrix
#' @exportMethod as.matrix
setGeneric('as.matrix', signature="x", function(x,...) standardGeneric('as.matrix'))
#' @rdname as
#' @aliases as.matrix,madness-method
setMethod('as.matrix', 'madness', function(x,...) { as(x,'matrix') })

#' @rdname as
#' @aliases as.numeric
#' @exportMethod as.numeric
#setGeneric('as.numeric', signature="x", function(x,...) standardGeneric('as.numeric'))
#' @rdname as
#' @aliases as.numeric,madness-method
setMethod('as.numeric', 'madness', function(x,...) { as(x,'numeric') })

#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
