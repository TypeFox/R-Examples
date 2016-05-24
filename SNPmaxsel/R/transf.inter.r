###    Transforms a pair of SNPs into a single variable with nine categories
###
### Copyright 2006-10 Anne-Laure Boulesteix 
###
### 
###
###
### This file is part of the `exactmaxstat' library for R and related languages.
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


transf.inter<-function(x1,x2)
{
x<-numeric(length(x1))

if (!(is.vector(x1)&is.vector(x2)))
 stop("x1 and x2 must be vectors")

if (!(is.numeric(x1)&is.numeric(x2)))
 stop("x1 and x2 must be numeric vectors")

if (length(x1)!=length(x2))
 stop("x1 and x2 must have the same length")

x1x2<-c(x1,x2)
x1x2<-x1x2[!is.na(x1x2)]
if (max(x1x2)>3|min(x1x2)<1)
 stop("x1 and x2 must have values 1,2,3")

 
x[x1==1&x2==1]<-1
x[x1==2&x2==1]<-2
x[x1==3&x2==1]<-3
x[x1==1&x2==2]<-4
x[x1==2&x2==2]<-5
x[x1==3&x2==2]<-6
x[x1==1&x2==3]<-7
x[x1==2&x2==3]<-8
x[x1==3&x2==3]<-9
x[x==0]<-NA

return(x)
}
