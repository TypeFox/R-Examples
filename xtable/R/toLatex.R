### xtable package
###
### Produce LaTeX and HTML tables from R objects.
###
### Copyright 2000-2013 David B. Dahl <dahl@stat.byu.edu>
###
### Maintained by David Scott <d.scott@auckland.ac.nz>
###
### This file is part of the `xtable' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

### The generic for toLatex() is declared in the base package "utils"

toLatex.xtable <- function(object, ...){
  ## Initially just capturing the output of print.xtable().  At some
  ## point this could be refactored to have print.xtable() call
  ## toLatex() instead. - CR, 30/01/2012
  dotArgs <- list(...)
  dotArgs$x <- object
  dotArgs$type <- "latex"
  dotArgs$print.results <- FALSE
  z <- do.call("print.xtable", dotArgs)

  z <- strsplit(z, split="\n")[[1]]
  class(z) <- "Latex"
  z
}
