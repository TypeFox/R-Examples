## Copyright 2013-2014 Sebastian Gibb
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

.importCiphergenXml <- function(file, centroided=FALSE, massRange=c(0, Inf),
                                minIntensity=0, verbose=FALSE) {

  .msg(verbose, "Reading spectrum from ", sQuote(file), " ...")

  if (!file.exists(file)) {
    stop("File ", sQuote(file), " doesn't exists!")
  }

  ## read file
  s <- .parseCiphergenXml(file=file)

  list(.createMassObject(mass=s$spectrum$mass, intensity=s$spectrum$intensity,
                         metaData=s$metaData, centroided=centroided,
                         massRange=massRange, minIntensity=minIntensity,
                         verbose=verbose))
}
