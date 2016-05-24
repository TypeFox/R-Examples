### variable.selection.R  (2005-04-06)
###
###     Variable selection using the PLS weights
###
### Copyright 2004-04 Anne-Laure Boulesteix and Korbinian Strimmer
###
### Part of the code was adopted from the pls.pcr package by Ron Wehrens
###
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

variable.selection<-function(X,Y,nvar=NULL)
{
Xscaled<-scale(X,center=TRUE,scale=TRUE)
Y<-as.numeric(Y)
Y<-Y-mean(Y)
a<-t(Xscaled)%*%Y
if (is.null(nvar))
 {
 nvar<-ncol(X)
 }
return(order(-abs(a))[1:nvar])
}
