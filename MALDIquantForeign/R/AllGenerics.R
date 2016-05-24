## Copyright 2012-2013 Sebastian Gibb
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

if (is.null(getGeneric(".composeFilename")))
  setGeneric(".composeFilename", function(x, ...) standardGeneric(".composeFilename"))
if (is.null(getGeneric(".exportTab")))
  setGeneric(".exportTab", function(x, ...) standardGeneric(".exportTab"))
if (is.null(getGeneric(".exportCsv")))
  setGeneric(".exportCsv", function(x, ...) standardGeneric(".exportCsv"))
if (is.null(getGeneric(".exportImzMl")))
  setGeneric(".exportImzMl", function(x, ...) standardGeneric(".exportImzMl"))
if (is.null(getGeneric(".exportMsd")))
  setGeneric(".exportMsd", function(x, ...) standardGeneric(".exportMsd"))
if (is.null(getGeneric(".exportMzMl")))
  setGeneric(".exportMzMl", function(x, ...) standardGeneric(".exportMzMl"))
if (is.null(getGeneric("export")))
  setGeneric("export", function(x, ...) standardGeneric("export"))
if (is.null(getGeneric("exportTab")))
  setGeneric("exportTab", function(x, ...) standardGeneric("exportTab"))
if (is.null(getGeneric("exportCsv")))
  setGeneric("exportCsv", function(x, ...) standardGeneric("exportCsv"))
if (is.null(getGeneric("exportImzMl")))
  setGeneric("exportImzMl", function(x, ...) standardGeneric("exportImzMl"))
if (is.null(getGeneric("exportMsd")))
  setGeneric("exportMsd", function(x, ...) standardGeneric("exportMsd"))
if (is.null(getGeneric("exportMzMl")))
  setGeneric("exportMzMl", function(x, ...) standardGeneric("exportMzMl"))
