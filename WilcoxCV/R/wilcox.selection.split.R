### wilcox.selection.split.R  (2007-03-26)
###
###     Wilcoxon rank sum statistic in leave-one-out cross-validation
###
### Copyright 2007-03 Anne-Laure Boulesteix 
###
### 
###
###
### This file is part of the `Wilcoxsamp' library for R and related languages.
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


wilcox.selection.split<-function(x,y,split,algo="new",pvalue=FALSE)
{

if (!all(is.element(y,c(0,1,NA))))
 stop("y must be coded as 0,1")

if (is.factor(y))
 {
 y<-as.numeric(y)-1
 }

n<-length(y)

if (nrow(x)!=length(y))
 stop("The length of y must equal the number of rows of x")

 n0<-sum(y==0)
 n1<-n-n0

 ntest<-ncol(split)
 niter<-nrow(split)
 ntrain<-n-ntest

if (algo=="new")
 {
 wilcox.split<-apply(x,FUN=wilcox.split.internal,MARGIN=2,y=y,split=split,algo="new",n=n,niter=niter)

#print(matrix(1-y[split],niter))
 n0.vector<-rep(n0,niter)-apply(matrix(1-y[split],niter),MARGIN=1,FUN=sum)
 E.vector<-n0.vector*(ntrain+1)/2
 SD.vector<-sqrt(E.vector*(ntrain-n0.vector)/6)
 result<-abs(apply(wilcox.split,FUN="-",MARGIN=2,E.vector))
 result<-apply(result,FUN="/",MARGIN=2,SD.vector)

 if (pvalue)
  {
  pvalue.split<-2*(1-pnorm(result))
  ordering.split<-t(apply(pvalue.split,FUN=order,MARGIN=1))
  return(list(ordering.split=ordering.split,pvalue.split=pvalue.split))
  }
 else
  {
  ordering.split<-t(apply(-result,FUN=order,MARGIN=1))
  return(list(ordering.split=ordering.split))
  }
 }
 

if (algo=="naive")
 {
 pvalue.split<-apply(x,FUN=wilcox.split.internal,MARGIN=2,y=y,split=split,algo="naive.pvalue",n=n,niter=niter)
 # wilcox.split<-apply(x,FUN=wilcox.split.internal,MARGIN=2,y=y,type="naive.stat",n0=n0,n1=n1,n=n)
 ordering.split<-t(apply(pvalue.split,FUN=order,MARGIN=1))

 if (pvalue)
  return(list(ordering.split=ordering.split,pvalue.split=pvalue.split))
 else
  return(list(ordering.split=ordering.split))
 }


}
################
