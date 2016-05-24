#  progress.R
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
# Functions for a progress bar and similar stuff
#
# .progressDots: Prints one dot per x loops to
#                the standard output.
#     Arguments: x     the number of loops per dot
#                i     the current value of the
#                      loop variable
#                total the maximum nuber of loops
#         Usage: For example in a for loop:
#                for (i in 1:n) {
#                    .progressDots(x = 10, i, total = n)
#                    # something time consuming
#                }

# .progressBar:  Prints a progress bar to the
#                standard output.
#    Arguments:  i          the current value of the
#                           loop variable
#                total      the maximum nuber of loops
#                symb_drawn the number of printed
#                           symbols for the bar
#        Usage:  for example in a for loop:
#                progr <- .progressBar()
#                for (i in 1:n) {
#                    progr <- .progressBar(i, total = n, symb_drawn = progr)
#                    # something time consuming
#                }

.progressDots <- function(x, i, total, symb = ".", num = 40, nl = "\n") {

  # x:     one dot per x loops
  # i:     loop variable
  # total: maximum value of i
  # symb:  printed symbol
  # num:   maximum number of dots per line
  # nl:    newline character

  if (i %% x == 0) {              # print symb after x loops
      cat(symb)
      if (i %% 5 == 0) {
          cat(" ")
      }
      if (i %% (x * num) == 0) {  # print a newline after num symbs
          cat(nl)
      }
  }
  if (i == total) {               # print a newline after the last symb
      cat(nl)
  }

}


.progressBar <- function (i, total, symb_drawn = 0, symb = "=", nl = "\n") {

  # i:          loop variable
  # total:      maximum value of i
  # symb_drawn: number of printet symbs
  # symb:       printed symbol
  # nl:         newline character

  symb_drawnNEW <- symb_drawn

  if (missing(i) && missing(total)) {
      init = TRUE
  }
  else {
      init = FALSE
  }

  if (isTRUE(init)) {
      # print a timeline
      cat("|            :            |            :            | 100 %", nl, "|", sep = "")
  }
  else {
      # print(paste("imax " , total, "nSymbols: ", symb_drawn, "i ", i))
      if (symb_drawn < 51) {

          #progress <- round((i/total) * 51, digits = 0) - symb_drawn
          progress <- floor((i/total) * 51) - symb_drawn
          #print(progress)
          symb_drawnNEW <- progress + symb_drawn

          cat(rep(symb, progress))
          #while (progress > 0) {
          #    cat(symb);
          #    progress <- progress - 1
          #}

          #return(symb_drawn)
      }
      #else {
          if (i >= total) {
              cat("| :-) \n")
          }
      #}
  }

  return(symb_drawnNEW)
}
