## Copyright 2015 Sebastian Gibb
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

.writeIbd <- function(x, filename, uuid, processed=TRUE) {
  ibd <- file(filename, open="wb")
  on.exit(close(ibd))

  uuid <- gsub("-", "", tolower(uuid))
  uuid <- substring(uuid, seq(1L, 32L, by=2L), seq(2L, 32L, by=2L))
  uuid <- strtoi(uuid, base=16L)
  writeBin(uuid, con=ibd, size=1L, endian="little")

  for (i in seq(along=x)) {
    if (processed || i == 1L) {
      writeBin(as.double(x[[i]]@mass), con=ibd, size=8L, endian="little")
    }
    writeBin(as.double(x[[i]]@intensity), con=ibd, size=8L, endian="little")
  }
}

.ibdOffsets <- function(x, processed=TRUE, encodedLengthSize=8L) {
  ## start at 16 (16 bytes for UUID)
  n <- rep(unlist(lapply(x, length)), each=2L)
  encodedLength <- n * encodedLengthSize

  if (processed) {
    offsets <- cumsum(c(16L, encodedLength[-length(n)]))
  } else {
    sel <- seq(from=2L, to=length(n), by=2L)
    offsets <- rep.int(16L, length(n))
    offsets[sel] <- 16L + cumsum(encodedLength[sel])
  }

  matrix(c(offsets, n, encodedLength), nrow=length(n),
         dimnames=list(rep(c("mass", "intensity"), times=length(x)),
                       c("offset", "length", "encodedLength")))
}

.addIbdOffsets <- function(x, processed=TRUE, encodedLengthSize=8L) {
  offsets <- .ibdOffsets(x, processed=processed,
                         encodedLengthSize=encodedLengthSize)
  i <- split(1:nrow(offsets), rep(1:length(x), each=2L))
  for (j in seq(along=x)) {
    x[[j]]@metaData$imaging$offsets <- offsets[i[[j]],]
  }
  x
}
