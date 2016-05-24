### entropy.R  (2013-07-16)
###
###    Estimating entropy from observed counts
###
### Copyright 2008-13 Korbinian Strimmer
###
###
### This file is part of the `entropy' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
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


entropy = function(y, lambda.freqs, method=c("ML", "MM", "Jeffreys", "Laplace", 
                   "SG", "minimax", "CS", "NSB", "shrink"),
                   unit=c("log", "log2", "log10"), verbose=TRUE, ...)
{
  method = match.arg(method)

  if (method == "ML")       H = entropy.empirical(y, unit=unit)
  if (method == "MM")       H = entropy.MillerMadow(y, unit=unit)
  if (method == "NSB")      H = entropy.NSB(y, unit=unit, ...)
  if (method == "CS")       H = entropy.ChaoShen(y, unit=unit)

  if (method == "Jeffreys") H = entropy.Dirichlet(y, a=1/2, unit=unit)
  if (method == "Laplace")  H = entropy.Dirichlet(y, a=1, unit=unit)
  if (method == "SG")       H = entropy.Dirichlet(y, a=1/length(y), unit=unit)
  if (method == "minimax")  H = entropy.Dirichlet(y, a=sqrt(sum(y))/length(y), unit=unit)

  if (method == "shrink")   H = entropy.shrink(y, lambda.freqs=lambda.freqs,
      unit=unit, verbose=verbose)

  return(H)
}

freqs = function(y, lambda.freqs, method=c("ML", "MM", "Jeffreys", "Laplace", 
                   "SG", "minimax", "CS", "NSB", "shrink"), verbose=TRUE)
{
  method = match.arg(method)

  if (method == "ML")       H = freqs.empirical(y)
  if (method == "MM")       H = rep(NA, length(y))
  if (method == "NSB")      H = rep(NA, length(y))
  if (method == "CS")       H = rep(NA, length(y))

  if (method == "Jeffreys") H = freqs.Dirichlet(y, a=1/2)
  if (method == "Laplace")  H = freqs.Dirichlet(y, a=1)
  if (method == "SG")       H = freqs.Dirichlet(y, a=1/length(y))
  if (method == "minimax")  H = freqs.Dirichlet(y, a=sqrt(sum(y))/length(y))

  if (method == "shrink")   H = freqs.shrink(y, lambda.freqs=lambda.freqs,
                               verbose=verbose)

  return(H)
}


