### entropy.ChaoShen.R  (2008-08-20)
###
###    Chao-Shen entropy estimator (2003)
###
### Copyright 2008 Korbinian Strimmer
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


# compute Chao-Shen (2003) entropy estimator

# y is a vector of counts (may include zeros)
entropy.ChaoShen = function(y, unit=c("log", "log2", "log10"))
{
  unit = match.arg(unit)

  yx = y[y > 0]           # remove bins with zero counts
  n = sum(yx)             # total number of counts
  p = yx/n                # empirical frequencies

  f1 = sum(yx == 1)       # number of singletons
  if (f1 == n) f1 = n-1   # avoid C=0

  C = 1 - f1/n            # estimated coverage
  pa = C*p                # coverage adjusted empirical frequencies
  la = (1-(1-pa)^n)       # probability to see a bin (species) in the sample

  H = -sum(pa*log(pa)/la) # Chao-Shen (2003) entropy estimator

  if (unit == "log2")       H = H/log(2)  # change from log to log2 scale
  if (unit == "log10")      H = H/log(10) # change from log to log10 scale

  return(H)
}

