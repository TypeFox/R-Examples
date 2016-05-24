{############################################################################### 
# spacing.patterns.R
# This file is part of the R package lint.
# 
# Copyright 2012 Andrew Redd
# Date: 6/16/2012
# 
# DESCRIPTION
# ===========
# predefined spacing patterns.
# 
# LICENSE
# ========
# lint is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
# 
# lint is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see http://www.gnu.org/licenses/.
# 
}###############################################################################

#' @rdname stylechecks
#' @export
#' @include styles.assignment.R
#' @include styles.performance.R
#' @include styles.spacing.R
lint.style <- list(
    spacing.twobeforecomments
  , spacing.spacearoundinfix
  , spacing.spacearoundequals
  , spacing.indentation.notabs
  , spacing.linelength.80
  , styles.assignment.noeq
  , styles.assignment.norightassign
  , styles.assignment.notinfcall
)




