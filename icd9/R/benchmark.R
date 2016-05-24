# Copyright (C) 2014 - 2015  Jack O. Wasey
#
# This file is part of icd9.
#
# icd9 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# icd9 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with icd9. If not, see <http:#www.gnu.org/licenses/>.

# nocov start

randomPatients <- function(...)
  randomOrderedPatients(...)

randomOrderedPatients <- function(...) {
  x <- randomUnorderedPatients(...)
  x[order(x$visitId), ]
}

randomUnorderedPatients <- function(num_patients = 50000, dz_per_patient = 20,
                                    n = num_patients, np = dz_per_patient) {
  set.seed(1441)
  pts <- round(n / np)
  data.frame(
    visitId = sample(seq(1, pts), replace = TRUE, size = n),
    icd9 = c(randomShortIcd9(round(n / 2)), randomShortAhrq(n - round(n / 2))),
    poa = as.factor(
      sample(x = c("Y", "N", "n", "n", "y", "X", "E", "", NA),
             replace = TRUE, size = n)),
    stringsAsFactors = FALSE
  )
}

#' genereate random short icd9 codes
#' @keywords internal
#' @importFrom stats runif
randomShortIcd9 <- function(n = 50000)
  as.character(floor(stats::runif(min = 1, max = 99999, n = n)))

randomShortAhrq <- function(n = 50000)
  sample(unname(unlist(icd9::ahrqComorbid)), size = n, replace = TRUE)

randomDecimalIcd9 <- function(n = 50000)
  paste(
    round(stats::runif(min = 1, max = 999, n = n)),
    sample(icd9ExpandMinor(), replace = TRUE, size = n),
    sep = "."
  )

runOpenMPVecInt <- function(n = 4, np = 2, threads = 6, chunkSize = 32) {
  pts <- randomPatients(n, np = np)
  icd9ComorbidShortCpp(pts, icd9::ahrqComorbid, threads = threads, chunkSize = chunkSize)
}
# nocov end
