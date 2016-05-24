### ust.R  (2014-10)
###
###    Univariate Soft Thresholding
###
### Copyright 2014-10 Ghislain DURIF
###
### Adapted from R package "spls"
### Reference: Chun H and Keles S (2010)
### "Sparse partial least squares for simultaneous dimension reduction and variable selection",
### Journal of the Royal Statistical Society - Series B, Vol. 72, pp. 3--25.
###
### This file is part of the `plsgenomics' library for R and related languages.
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



ust <- function(w, lambda.l1) {
	
	## output: thresholded column vector
	w.th <- matrix(0, nrow=length(w), ncol=1)
	
	## soft thresholding
	val.w <- abs(w) - lambda.l1 * max( abs(w) )
	w.th[ val.w >=0 ] <- val.w[ val.w>=0 ] * (sign(w))[ val.w>=0 ]
	
	## return
	return(w.th)
	
}