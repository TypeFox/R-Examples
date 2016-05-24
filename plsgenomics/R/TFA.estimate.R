### TFA.estimate.R  (2005-04-11)
###
###     Prediction of Transcription Factor Activities using PLS
###
### Copyright 2004-04 Anne-Laure Boulesteix and Korbinian Strimmer
###
### 
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



TFA.estimate<-function(CONNECdata,GEdata,ncomp=NULL,nruncv=0,alpha=2/3,unit.weights=TRUE)
{
n<-nrow(GEdata)
m<-ncol(GEdata)
p<-ncol(CONNECdata)
if (nruncv==0&length(ncomp)>1)
 stop("Since length(ncomp)>1, nruncv must be >0")

if (is.null(ncomp))
 {
 ncomp<-min(n,p)
 }
if (n!=nrow(CONNECdata))
 stop("The number of genes must be the same in the gene expression data and in the CONNEC data")
 
X<-scale(CONNECdata,center=TRUE,scale=TRUE)
X[is.na(X)]<-0
Y<-scale(GEdata,center=TRUE,scale=TRUE)
if (nruncv>0)
 {
 ncomp<-pls.regression.cv(X,Y,ncomp=ncomp,nruncv=nruncv) 
 }

Dx<-diag(1/attributes(X)$"scaled:scale")
Dy<-diag(1/attributes(Y)$"scaled:scale")
pls.out<-pls.regression(X,Y,ncomp=ncomp,Xtest=NULL,unit.weights=unit.weights)
TFA<-Dx%*%pls.out$B%*%solve(Dy)
metafactor<-pls.out$Q

return(list(TFA=TFA,metafactor=metafactor,ncomp=ncomp))
}

