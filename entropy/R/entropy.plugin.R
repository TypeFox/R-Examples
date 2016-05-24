### entropy.plugin.R  (2008-09-28)
###
###    Plug-in entropy estimator
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


entropy.plugin = function(freqs, unit=c("log", "log2", "log10"))
{
   unit = match.arg(unit)

   freqs = freqs/sum(freqs) # just to make sure ...

   H = -sum( ifelse(freqs > 0, freqs*log(freqs), 0) )

   if (unit == "log2")  H = H/log(2)  # change from log to log2 scale
   if (unit == "log10") H = H/log(10) # change from log to log10 scale

   return(H)
}

