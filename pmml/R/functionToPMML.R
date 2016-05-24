# PMML: Predictive Model Markup Language
#
# Copyright (c) 2009-2015, some parts by Togaware Pty Ltd and other by Zementis, Inc. 
#
# This file is part of the PMML package for R.
#
# The PMML package is free software: you can redistribute it and/or 
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 2 of 
# the License, or (at your option) any later version.
#
# The PMML package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. Please see the
# GNU General Public License for details (http://www.gnu.org/licenses/).
######################################################################################


#' Convert an R expression to PMML.
#' 
#' @param expr an R expression enclosed in quotes

#' @details 
#' As long as the expression passed to the function is a valid R expression (e.g., no unbalanced parenthesis), 
#' it can contain arbitrary function names not defined in R. Variables in the expression passed 
#' to `FunctionXform` are always assumed to be fields, and not substituted. That is, even if `x` has a value in the 
#' R environment, the resulting expression will still use `x`.
#' 
#' An expression such as `foo(x)` is treated as a function `foo` with argument `x`. Consequently, passing in an 
#' R vector `c(1,2,3)` to `functionToPMML()` will produce PMML where `c` is a function and `1,2,3` are the arguments.
#'
#' @return PMML version of the input expression
#' 
#' @author Dmitriy Bolotov
#' 
#' @examples
#' 
#' # Operator precedence and parenthesis
#' functionToPMML("1 + 3/5 - (4 * 2)")
#' 
#' # Nested arbitrary functions
#' functionToPMML("foo(bar(x)) - bar(foo(y-z))")
#' 
#' # If-else expression
#' functionToPMML("if (x==3) { 3 } else { 0 }")
#' 
#' # Function with string argument types
#' functionToPMML("colors('red','green','blue')")

# PMML (Predictive Model Markup Language) Transformations 
#
# Copyright (c) 2015 Zementis, Inc.
#
# This file is part of the pmmlTransformations package 
#
# The pmmlTransformations package is free: you can redistribute it and/or 
# modify it under the terms of the GNU General Public License as published 
# by the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# The pmmlTransformations package is distributed in the hope that it will 
# be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. Please see the
# GNU General Public License for details (http://www.gnu.org/licenses/).
############################################################################


functionToPMML <- function(expr) {
  
  xformed <- .pmmlU(expr)
  
  return(xformed)
}
