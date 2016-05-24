###    Generating learning sets 
###
### Copyright 2007-12 Anne-Laure Boulesteix 
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
 
generate.learningsets<-function(n,method=c("LOOCV","CV","MCCV","bootstrap"),fold=NULL,niter=NULL,nlearn=NULL)
{

method<-match.arg(method,c("LOOCV","CV","MCCV","bootstrap"))

if (!is.null(fold)&method!="CV")
 warning("argument fold is ignored when method is not CV ")

if (!is.null(nlearn)&is.element(method,c("LOOCV","CV")))
 warning("argument nlearn is ignored when method is LOOCV or CV")


if (method=="MCCV")
 {
 if (is.null(niter)|is.null(nlearn))
  stop("With the MCCV method, arguments niter and nlearn should be given.")
 output.matrix<-matrix(0,niter,nlearn)
 for (i in 1:niter)
  {
  output.matrix[i,]<-sample(n,nlearn,replace=FALSE)
  }
 return(output.matrix)
 }


if (method=="CV")
 {
 if (fold==1)
  {
  return(matrix(1:n,nrow=1))
  }
 if (is.null(fold))
  stop("With the CV method, argument fold should be given.")
 if (fold==n)
  {
  method<-"LOOCV"
  }
 else
  {
  if (is.null(niter))
   {
   niter<-1
   }
  size<-n/fold
  output.matrix<-matrix(0,niter*fold,n-ceiling(size))

  if (size<5)
   stop("fold is too large")

  size.int<-floor(size)
  size.vector<-rep(size.int,fold)
 
  if (size.int!=size)
   {
   size.vector[1:((size-size.int)*fold)]<-size.vector[1:((size-size.int)*fold)]+1
   }
  group.index<-c()
  for (j in 1:fold)
   {
   group.index<-c(group.index,rep(j,size.vector[j]))
   }

  for (i in 1:niter)
   {
   group.index<-group.index[sample(n,n,replace=FALSE)]
   for (j in 1:fold)
    {
    whichj<-which(group.index==j)
    output.matrix[j+(i-1)*fold,1:(n-length(whichj))]<-setdiff(1:n,whichj)
    }
   }
  }
 return(output.matrix)
 }


if (method=="LOOCV")
 {
 output.matrix<-matrix(0,n,n-1)
 for (i in 1:n)
  {
  output.matrix[i,]<-(1:n)[-i]
  }
 return(output.matrix)
 }


if (method=="bootstrap")
 {
 if (is.null(nlearn))
  {
  nlearn<-n
  }
 if (is.null(niter))
  stop("With the bootstrap method, argument niter should be given.")
 output.matrix<-matrix(0,niter,nlearn)
 for (i in 1:niter)
  {
  output.matrix[i,]<-sample(n,nlearn,replace=TRUE)
  }
 return(output.matrix)
 }





}