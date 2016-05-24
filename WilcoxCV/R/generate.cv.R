### generate.cv.R  (2007-03-26)
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


generate.cv<-function(n,m)
{

size<-n/m
cv<-matrix(0,m,ceiling(size))

if (size<5&size!=1)
 stop("m is too large")

if (size==1)
 return(matrix(1:n,nrow=n))

size.int<-floor(size)
size.vector<-rep(size.int,m)
 
if (size.int!=size)
 {
 size.vector[1:((size-size.int)*m)]<-size.vector[1:((size-size.int)*m)]+1
 }
group.index<-c()
for (j in 1:m)
 {
 group.index<-c(group.index,rep(j,size.vector[j]))
 }

group.index<-group.index[sample(n,n,replace=FALSE)]

for (j in 1:m)
 {
 whichj<-which(group.index==j)
 print(length(whichj))
 cv[j,1:length(whichj)]<-whichj
 }

return(cv)
}