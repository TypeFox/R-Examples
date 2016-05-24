### preprocess.R  (2005-04-06)
###
###    a pre-processing box for standard Microarrays data set
###
### Copyright 2006-01 Sophie Lambert-Lacroix and Julie Peyre
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

preprocess <- function(Xtrain, Xtest=NULL,Threshold=c(100,16000),Filtering=c(5,500),log10.scale=TRUE,row.stand=TRUE)
{
##    INPUT VARIABLES
#########################
## Xtrain   : matrix ntrain x p
##      train data matrix
## Xtest   : matrix ntest x p
##      test data matrix
## Threshold : vector of length 2 
##  values for thresholding data in preprocess (ex = c(100,16000))
## Filtering : vector of length 2 
##  values for filtering data in preprocess (ex = c(5,500))
## log10.scale : logical
##      should data be log10-transformed in preprocess ?
## row.stand : logical
##      should data be standardized in row in preprocess ?
##
##
##  OUTPUT variables
#####################
## Structure with fields
##    pXtrain : matrix ntrain x p
##              preprocessed train data matrix
##    pXtest : matrix ntest x p
##              preprocessed test data matrix


## TEST ON THE INPUT VARIABLES
#################################
## on X train
if ((is.matrix(Xtrain)==FALSE)||(is.numeric(Xtrain)==FALSE))
   stop("Message from preprocess.R: Xtrain is not of valid type")

if (dim(Xtrain)[2]==1)
  stop("Message from preprocess.R: p=1 is not valid")

ntrain <- dim(Xtrain)[1]
p <- dim(Xtrain)[2]

## on Xtest
ntest <- 0
if (is.null(Xtest)==FALSE)
{
   if (is.vector(Xtest)==TRUE)
      Xtest <- matrix(Xtest,nrow=1)

   if ((is.matrix(Xtest)==FALSE)||(is.numeric(Xtest)==FALSE))
      stop("Message from preprocess.R: Xtest is not of valid type")

   if (dim(Xtrain)[2]!=dim(Xtest)[2])
      stop("Message from preprocess.R: columns of Xtest and columns of Xtrain must be equal")
  
   ntest <- dim(Xtest)[1] 
}   

## on Threshold     
Thresh <- FALSE
if (is.null(Threshold) == FALSE)
{
   if (is.numeric(Threshold)==FALSE | is.vector(Threshold)==FALSE)
      warning("Message from preprocess.R: no thresholding is done in preprocess.")
   else
      if ((length(Threshold)!=2) | Threshold[2] <= Threshold[1])
      {
         warning("Message from preprocess.R: threshold argument in preprocess is not correct, no thresholding is done in preprocess.")
      }
      else
         Thresh <- TRUE
}

## on Filtering
Filter <- FALSE
if (is.null(Filtering) == FALSE)
{
   if (is.numeric(Filtering)==FALSE | is.vector(Filtering)==FALSE)
      warning("Message from preprocess.R: no filtering is done in preprocess.")
   else
      if (length(Filtering)!=2)
      {
         warning("Message from preprocess.R: filtering argument in preprocess is not correct, no filtering is done in preprocess.")
      }
      else
         Filter <- TRUE
}

## on log10.scale
if ((is.logical(log10.scale) == FALSE) | (length(log10.scale) != 1))
{
   warning("Message from preprocess.R: log10.scale argument in preprocess is not correct. No log10 transform is done in preprocess",.call=FALSE)
   log10.scale <- FALSE
}
   
## on row.stand
if ((is.logical(row.stand) == FALSE) | (length(row.stand) != 1))
{
   warning("Message from preprocess.R: row.stand argument in preprocess is not correct. No standardization in row is done in preprocess")
   row.stand <- FALSE
}
   
# check colnames before creating DATA
if (ntest != 0)
{
   # if the two colnames are null
   if (is.null(colnames(Xtrain)) && is.null(colnames(Xtest)))
   {
      colnames(Xtrain) <- paste("g",1:p,sep="")
      colnames(Xtest) <- paste("g",1:p,sep="")
   }
   else
   {
      # if only one of the two is null then no problem

         # if the two colnames are not null
         if (is.null(colnames(Xtrain))==FALSE && is.null(colnames(Xtest))==FALSE)
         {
            # check the colnames are the same
            if (sum(colnames(Xtrain)==colnames(Xtest))!=p)
               stop("Message from preprocess.R: error in data. Genes must be the same (same names) in train and test sets.")
         }
      }
   }

   DATA <- rbind(Xtrain,Xtest)
   names <- colnames(DATA)
   rm(Xtrain)
   rm(Xtest)
   
## PREPROCESS
##############
## Threshold
   if (Thresh == TRUE)
   {
    # floor and ceil
    if (sum(DATA < Threshold[1]) > 0)
       DATA[DATA < Threshold[1]] <- Threshold[1]
    
    if (sum(DATA > Threshold[2]) > 0)
       DATA[DATA > Threshold[2]] <- Threshold[2]
   }
     
## Filtering
   if (Filter == TRUE)
   {
    minval <- apply(DATA[1:ntrain,],2,min)
    maxval <- apply(DATA[1:ntrain,],2,max)
    genei <- as.numeric(maxval/minval>Filtering[1])*as.numeric(maxval-minval>Filtering[2])
    genei <- which(genei==1)
    names <- names[genei]
    
    if (length(genei)==0)
    {
       stop("Message from preprocess.R: check if the data is not already preprocessed. No filtering can be executed on that data set.")
    }
    else
    {
       DATA <- DATA[,genei]
    }
   }
         
##  Log10 transform
   if (log10.scale == TRUE)
   {
    DATA <- log10(DATA)
   } 
   
## standardisation in row
   if (row.stand == TRUE)
   {
    # for DATA
    mean.att<-apply(DATA,1,mean)
    DATA<-sweep(DATA,1,mean.att,FUN="-") 
    std.att<-sqrt((length(DATA[1,])-1)*apply(DATA,1,var)/length(DATA[1,]))
    DATA<-sweep(DATA,1,std.att,FUN="/")     
   }
   
   colnames(DATA) <- names

   pXtrain <- DATA[1:ntrain,]

   if (ntest != 0)
      pXtest <- DATA[(ntrain+1):(ntrain+ntest),]
   else
      pXtest <- NULL   
     
return(list(pXtrain=pXtrain,pXtest=pXtest))
}
