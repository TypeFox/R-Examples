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


testclass_simul<-function(datalist,nlearn=100,classifier,ncomp=0:3,nbgene=NULL,varsel=NULL,fold=10,...)
{
importance<-list()
niter<-length(datalist)
p<-ncol(datalist[[1]]$x)
q<-ncol(datalist[[1]]$z)
n<-nrow(datalist[[1]]$z)
error<-numeric(niter)
masig<-logical(niter)
bestncomp<-numeric(niter)
OOB<-list()

for (i in 1:niter)
 {
  print(i)

  ylearni<-as.numeric(factor(datalist[[i]]$y[1:nlearn]))-1
  ytesti<-as.numeric(factor(datalist[[i]]$y[(nlearn+1):n]))-1

  if (!is.null(datalist[[i]]$x))
   {
   xi<-datalist[[i]]$x
   xlearni<-xi[1:nlearn,]
   xtesti<-xi[(nlearn+1):n,]
   }
  else
   {
   xlearni<-NULL
   xtesti<-NULL
   }

  if (!is.null(datalist[[i]]$z))
   {
   zi<-datalist[[i]]$z
   zlearni<-zi[1:nlearn,]
   ztesti<-zi[(nlearn+1):n,]
   }
  else
   {
   zlearni<-NULL
   ztesti<-NULL
   }

 result<-classifier(Xlearn=xlearni,Ylearn=ylearni,Zlearn=zlearni,Xtest=xtesti,Ztest=ztesti,ncomp=ncomp,ordered=varsel[i,],nbgene=nbgene,fold=fold,...)
 prediction<-result$prediction
 error[i]<-sum(ytesti!=prediction)
 if (!is.null(result$bestncomp))
  {
  bestncomp[i]<-result$bestncomp
  }
 if (is.null(result$OOB))
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



