#  reassignFwBwMatch.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybil.
#
#  Sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.


################################################
# Function: .reassignFwBwmatch
#
#
# The function .reassignFwBwmatch() is inspired by the function
# reassingFwBwMatch() contained in the COBRA Toolbox.
# The algorithm is the same.


.reassignFwBwMatch <- function(matchrev, keep) {

  ind               <- keep * 1
  ind[keep == TRUE] <- 1:sum(ind)

  num_match         <- sum(keep == TRUE)
  match             <- integer(num_match)
  j                 <- 0

  for (i in 1:length(matchrev)) {
      if (keep[i] == TRUE) {
          j <- j + 1
          if (matchrev[i] > 0) {
              if (keep[matchrev[i]] == TRUE) {
                  match[j] <- ind[matchrev[i]]
              }
          }
      }
  }

  return(as.integer(match))
}
