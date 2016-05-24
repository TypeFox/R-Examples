###    Testing a classifier based on microarray data and clinical parameters
###
### Copyright 2007-11 Anne-Laure Boulesteix 
###
### 
###
###
### This file is part of the `MAclinical' library for R and related languages.
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


testclass<-function(x=NULL,y,z=NULL,learningsets,classifier,ncomp=0:3,varsel=NULL,nbgene=NULL,fold=10,...)
{

n<-length(y)
p<-ncol(x)
niter<-nrow(learningsets)
error<-numeric(niter)
bestncomp<-numeric(niter)
OOB<-list()
importance<-list()
y<-as.numeric(factor(y))-1


for (i in 1:niter)
 {
 print(i)
 ylearni<-y[learningsets[i,]]
 ytesti<-y[-learningsets[i,]]

 if (!is.null(x))
  {
  xtesti<-x[-learningsets[i,],]
  xlearni<-x[learningsets[i,],]
  }
 else
  {
  xlearni<-NULL
  xtesti<-NULL
  }

 if (!is.null(z))
  {
  ztesti<-z[-learningsets[i,],]
  zlearni<-z[learningsets[i,],]
  }
 else
  {
  zlearni<-NULL
  ztesti<-NULL
  }


 result<-classifier(Xlearn=xlearni,Ylearn=ylearni,Zlearn=zlearni,Xtest=xtesti,Ztest=ztesti,ordered=varsel[i,],ncomp=ncomp,nbgene=nbgene,fold=fold,...)
 prediction<-result$prediction
print(prediction)
 error[i]<-sum(ytesti!=prediction)
 if (!is.null(result$bestncomp))
  {
  bestncomp[i]<-result$bestncomp
  }
 if (!is.null(result$OOB))
  {
  OOB[[i]]<-result$OOB
  }
 }

if (is.null(result$OOB))
 {
 return(list(error=error))
 }

if (!is.null(result$OOB)&is.null(result$bestncomp))
 {
 return(list(error=error,OOB=OOB))
 }
else
 {
 return(list(error=error,bestncomp=bestncomp,OOB=OOB))
 }

}


