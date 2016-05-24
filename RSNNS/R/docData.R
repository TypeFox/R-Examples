#############################################################################
#
#   This file is part of the R package "RSNNS".
#
#   Author: Christoph Bergmeir
#   Supervisor: José M. Benítez
#   Copyright (c) DiCITS Lab, Sci2s group, DECSAI, University of Granada.
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Library General Public
#   License as published by the Free Software Foundation; either
#   version 2 of the License, or (at your option) any later version.
# 
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Library General Public License for more details.
# 
#   You should have received a copy of the GNU Library General Public License
#   along with this library; see the file COPYING.LIB.  If not, write to
#   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
#   Boston, MA 02110-1301, USA.
#
#############################################################################


#' This is data from the original SNNS examples directory ported to R and stored as one list.
#' The function \code{\link{readPatFile}} was used to parse all pattern files (.pat) from the 
#' original SNNS examples directory. Due to limitations of that function, pattern files
#' containing patterns with variable size were omitted. 
#'
#' @title Example data of the package
#' @name snnsData
#' @docType data
#' @keywords data
#' @examples
#' data(snnsData)
#' names(snnsData)
NULL