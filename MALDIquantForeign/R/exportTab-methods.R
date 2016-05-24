## Copyright 2012 Sebastian Gibb
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

setMethod(f=".exportTab",
  signature=signature(x="AbstractMassObject"),
  definition=function(x, file="", row.names=FALSE, col.names=FALSE, ...) {
  return(write.table(as.matrix(x), file=file,
                     row.names=row.names, col.names=col.names, ... ))
})

setMethod(f=".exportCsv",
  signature=signature(x="AbstractMassObject"),
  definition=function(x, file="", sep=",", row.names=FALSE, col.names=TRUE,
                      ...) {
  return(.exportTab(x, file=file, sep=sep, row.names=row.names,
                    col.names=col.names, ...))
})
