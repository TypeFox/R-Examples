### wilcox.split.R  (2007-03-26)
###
###     Wilcoxon rank sum statistic in subsampling
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


wilcox.split<-function(x,y,split,algo="new")
{


n<-length(y)
if (!all(is.element(y,c(0,1))))
 stop("y must be coded as 0,1")

if (length(x)!=length(y))
 stop("y and x must have the same length")


n<-length(x)
niter<-nrow(split)

if (algo=="new")
 {
 O<-order(x)
 R<-numeric(n)
 R[O]<-1:n
 R0<-R*(1-y)
 W0<-sum(R0)
 mat<-lower.tri(diag(n))
 Rmat<-mat[R,R] 
 Rmat[y==1,]<-0
 W.split<-numeric(niter)
 
 for (j in 1:niter)
  {
  W.split[j]<-W0-sum(R0[split[j,]])-sum(Rmat[-split[j,],split[j,]])
  }
 return(W.split)
 }

if (algo=="naive")  
 {
 W.split<-numeric(niter)

 for (j in 1:niter)
  {
  W.split[j]<-wilcox.test2(x[-split[j,]],y[-split[j,]],type="stat")
  }
 return(W.split)
 }

}
################

wilcox.split.internal<-function(x,y,split,algo="new",n,niter)
{

if (algo=="new")
 {
 O<-order(x)
 R<-numeric(n)
 R[O]<-1:n
 R0<-R*(1-y)
 W0<-sum(R0)
 mat<-lower.tri(diag(n))
 Rmat<-mat[R,R] 
 Rmat[y==1,]<-0
 W.split<-numeric(niter)
 
 for (j in 1:niter)
  {
  W.split[j]<-W0-sum(R0[split[j,]],na.rm=TRUE)-sum(Rmat[-split[j,],split[j,]],na.rm=TRUE)
  }
  return(W.split)
 }
 

if (algo=="naive.stat")
 {
 W.split<-numeric(niter)

 for (j in 1:niter)
  {
  W.split[j]<-wilcox.test2(x[-split[j,]],y[-split[j,]],type="stat")
  }
 return(W.split)
 }


if (algo=="naive.pvalue")
 {
 W.split<-numeric(niter)

 for (j in 1:niter)
  {
  W.split[j]<-wilcox.test2(x[-split[j,]],y[-split[j,]],type="pvalue")
  }
 }

return(W.split)
}

################
wilcox.test2<-function(x,y,type)
{

if (type=="stat")
 {
 n0<-sum(y==0)
 return(wilcox.test(x[y==0],x[y==1],correct=FALSE,exact=FALSE)$statistic+n0*(n0+1)/2)
 }

if (type=="pvalue")
 {
 return(wilcox.test(x[y==0],x[y==1],correct=FALSE,exact=FALSE)$p.value)
 }
}
#################

