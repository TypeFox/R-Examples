### This file is part of 'PGRdup' package for R.

### Copyright (C) 2014, ICAR-NBPGR.
#
# PGRdup is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# PGRdup is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#' Prints summary of \code{ProbDup} object.
#' 
#' \code{print.ProbDup} prints to console the summary of an object of class 
#' \code{ProbDup} including the method used ("a", "b" or "c"), the database
#' fields(columns) considered, the number of probable duplicate sets of each
#' kind along with the corresponding number of records.
#' 
#' @param x An object of class \code{ProbDup}.
#' @param ... Unused
#' @seealso \code{\link[PGRdup]{ProbDup}}
#'   
#' @export
print.ProbDup <- function(x,...) {
  attr <- attributes(x)
  x <- x[!sapply(x, is.null)]
  attributes(x) <- append(attributes(x), attr[2:4])
  m <- data.frame(`No. of Sets` = sapply(x, function(x) dim(x)[1]),
                  `No. of Records` = sapply(x,
                                            function(x) length(unique(unlist(strsplit(x$ID, split = ", "))))))
  Total <- as.character(colSums(m))
  Total[2] <- paste(Total[2], "(Distinct:",
                    length(unique(unlist(lapply(x,
                                                function(x) unlist(strsplit(x$ID,split = ", ")))))),
                    ")", sep = "")
  m <- rbind(m, Total = Total)
  cat(paste("Method : ", attributes(x)$method, "\n", sep = ""))
  cat(paste("\n", "KWIC1 fields : ", sep = ""))
  cat(paste(attributes(x)$fields$k1, sep = ""))
  if (attributes(x)$method != "a") {
    cat(paste("\n", "\n", "KWIC2 fields : ", sep = ""))
    cat(paste(attributes(x)$fields$k2, sep = ""))
  }
  cat("\n", "\n")
  print(m)
}
