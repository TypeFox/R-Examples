### generate.split.R  (2007-03-26)
###
###     Generating random cvtings training/test
###
### Copyright 2007-03 Anne-Laure Boulesteix 
###
### 
###
###
### This file is part of the `WilcoxCV' library for R and related languages.
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


generate.split<-function(niter,n,ntest)
{

split<-matrix(0,niter,ntest)
for (i in 1:niter)
 {
 split[i,]<-sample(n,ntest,replace=FALSE)
 }

return(split)
}