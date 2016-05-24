### entropy.empirical.R  (2013-07-16)
###
###    Empirical estimators of entropy, mutual information and related
###    quantities.
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


# compute empirical entropy 
# y is a vector of counts (may include zeros)
entropy.empirical = function(y, unit=c("log", "log2", "log10"))
{
  return( entropy.plugin(freqs.empirical(y), unit=unit) )
}

# empirical frequencies
freqs.empirical = function(y)
{
  return( y/sum(y) )
}


# empirical mutual information
mi.empirical = function(y2d, unit=c("log", "log2", "log10"))
{
  return( mi.plugin(freqs.empirical(y2d), unit=unit) )
}

# empirical chi-squared of independence
chi2indep.empirical = function(y2d, unit=c("log", "log2", "log10"))
{
  return( chi2indep.plugin(freqs.empirical(y2d), unit=unit) )
}


# empirical chi-squared statistic
chi2.empirical = function(y1, y2, unit=c("log", "log2", "log10"))
{
  return( chi2.plugin(freqs.empirical(y1), freqs.empirical(y2), unit=unit) )
}

# empirical KL divergence
KL.empirical = function(y1, y2, unit=c("log", "log2", "log10"))
{
  return( KL.plugin(freqs.empirical(y1), freqs.empirical(y2), unit=unit) )
}


