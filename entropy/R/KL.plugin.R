### KL.plugin.R  (2013-07-16)
###
###    Plug-in estimator of the Kullback-Leibler divergence 
###       and of the Chi-Squared Statistic
###
### Copyright 2013 Korbinian Strimmer
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


KL.plugin = function(freqs1, freqs2, unit=c("log", "log2", "log10") )
{
  unit = match.arg(unit)
  freqs1 = freqs1/sum(freqs1) # just to make sure ...
  freqs2 = freqs2/sum(freqs2) # just to make sure ...

  if( any( !(freqs2 > 0) ) ) warning("Vanishing value(s) in argument freqs2!")

  LR = ifelse(freqs1 > 0, log(freqs1/freqs2), 0)
  KL = sum(freqs1*LR)

  if (unit == "log2")  KL = KL/log(2)  # change from log to log2 scale
  if (unit == "log10") KL = KL/log(10) # change from log to log10 scale

  return(KL)
}


chi2.plugin = function(freqs1, freqs2, unit=c("log", "log2", "log10") )
{
  unit = match.arg(unit)
  freqs1 = freqs1/sum(freqs1) # just to make sure ...
  freqs2 = freqs2/sum(freqs2) # just to make sure ...

  if( any( !(freqs2 > 0) ) ) warning("Vanishing value(s) in argument freqs2!")

  chi2 = sum( (freqs1-freqs2)^2/freqs2 )

  if (unit == "log2")  chi2 = chi2/log(2)  # change from log to log2 scale
  if (unit == "log10") chi2 = chi2/log(10) # change from log to log10 scale

  return(chi2)
}

