# PMML: Predictive Modelling Markup Language
#
# Part of the Rattle package for Data Mining
#
# Dummy pmmltoc
#
# Time-stamp: <2010-12-03 06:00:43 Graham Williams>
#
# Copyright (c) 2009 Togaware Pty Ltd
#
# This files is part of the Rattle suite for Data Mining in R.
#
# Rattle is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Rattle is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rattle. If not, see <http://www.gnu.org/licenses/>.

########################################################################
# MAIN ENTRY POINT - SIGNATURE

pmmltoc <- function(p,
                    name=NULL,
                    includePMML=TRUE,
                    includeMetaData=TRUE,
                    exportClass=TRUE)
{
  # 101203 Use stop() rather than warning() since it should be brought
  # directly to the user's attention that the export to C is not
  # functional.
  
  stop("The PMMLtoC functionality is not implemented.")
}
