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

## generate uuid needed for imzML idb files
## based on dplR::uuid.gen by Mikko Korpela
## should be compatible to:
## http://en.wikipedia.org/wiki/Universally_unique_identifier#Version_4_.28random.29
## @param init used to create the md5 string, just for testing; don't change the
## default in other use cases.
## @return uuid string
## @noRd
.uuid <- function(init=paste(c(Sys.info(), Sys.getpid(), unlist(R.version), 
                               Sys.getpid(), 
                               format(Sys.time(), "%Y%m%d%H%M%OS6 %Z"), 
                               runif(5L)), collapse = "")) {
  md5 <- digest::digest(init, algo="md5", serialize=FALSE)
  md5 <- strsplit(md5, "", fixed=TRUE)[[1L]]
  md5[13] <- "4"
  md5[17] <- c("8", "9", "a", "b")[strtoi(md5[17L], base=16L)%%4L + 1L]
  .paste0 <- function(...)paste0(..., collapse="")
  paste(.paste0(md5[1L:8L]), 
        .paste0(md5[9L:12L]), 
        .paste0(md5[13L:16L]), 
        .paste0(md5[17L:20L]),
        .paste0(md5[21L:32L]), sep="-")
}

## test uuid version 4
## @param x character vector to test
## @return logical
## @noRd
.isUuidV4 <- function(x) {
  grepl(paste0("^[0-9a-f]{8}-?",
               "[0-9a-f]{4}-?",
               "4[0-9a-f]{3}-?",
               "[89ab][0-9a-f]{3}-?",
               "[0-9a-f]{12}"), tolower(x))
}
