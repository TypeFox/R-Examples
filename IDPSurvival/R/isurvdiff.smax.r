## File IDPSurvival/R/isurvdiff.smax.r
##
## IDPSurvival package for R (http://www.R-project.org)
##
## Copyright (C) 2014 Dalle Molle Institute for Artificial Intelligence.
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

isurvdiff.smax <- function(formula,...,verbose=FALSE,accuracy=0.05,smax=12) {
  s <- 0
  if(verbose) print(paste('Trying s=',s))
  out <- isurvdiff(formula, ...,s=s, display=FALSE)
  if (out$h==2) {
      s=-1
      return(list("s"=s,"test0"=out))
  }
  ds <- 3
  while(out$h!=2 && s<smax) {
      s <- s + ds
      if(verbose) print(paste('Trying s=',s))
      test0 <- out
      out <- isurvdiff(formula, ...,s=s, display=FALSE)
  }
  if(out$h!=2) {
      return(list("s"=s,"test0"=out))
  }
  mins <- s-ds
  maxs <- s
  while(mins < maxs - accuracy) {
      s <- (mins+maxs)/2
      if(verbose) print(paste('Trying s=',s))
      out <- isurvdiff(formula, ...,s=s,display=FALSE)
      if (out$h==2) maxs <- s
      else {
          mins <- s
          test0 <- out
      }
  }
  return(list("s"=mins, "test0"=test0))
}
