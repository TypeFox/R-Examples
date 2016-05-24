
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License 
## as published by the Free Software Foundation; either version 2 
## of the License, or (at your option) any later version.
##  
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General 
## Public License for more details.
##   
## You may also obtain a copy of the GNU General Public License from 
## the Free Software Foundation by visiting their Web site or by writing to
##   
##
##  Free Software Foundation, Inc.
##  59 Temple Place - Suite 330
##  Boston, MA 02111-1307
##  USA
##
################################################################################

## Function for generating nPt by nBiom data matrix for classification analysis 
## ----------------------------------------------------------------------------

## nTrain = Number of subjects in the training set (nGr1+nGr2)
## nGr1 = Number of subjects in Group 1
## nGr2 = Number of subjects in Group 2 
## nBiom = Number of biomarkers (genes, proteins)
## nRep= Number of technical replications

## sdW= sqrt(experimental or technical variation)
## sdB= sqrt(biological variation)
## rhoMax= Maximum Pearson's correlation coefficient between biomarkers
## rhoMinx= Minimum Pearson's correlation coefficient between biomarkers
##
## foldMin= Minimum value of fold changes
## nBlock= number of blocks in the block diagonal (Hub-Toeplitz) correlation matrix
## bSizes= a vector of block sizes (should sum to nBlock)
## gamma= if NULL, assume independence. gamma>=0 specify a correlation structure. gamma=0 indicates a single block exchangeable correlation marix with constant correlation rho=0.5*(rhoMin+rhoMax). A value greater than zero (e.g., gamma=1) block diagonal (Hub-Toeplitz) correlation matrix.
## sigma= Standard deviation of the normal distribution (before truncation)
##        where fold changes are generated from  
## baseExpr = A vector of length nBiom to be used as \mu_g (log2 scale, <16) 

## Main function for simulating data
## ---------------------------------


simData<-function(nTrain=100, nGr1=floor(nTrain/2), nBiom=50,nRep=3,
                  sdW=1.0, sdB=1.0 ,rhoMax=NULL, rhoMin=NULL, nBlock=NULL,bsMin=3, bSizes=NULL, gamma=NULL,
                  sigma=0.1,diffExpr=TRUE, foldMin=2,orderBiom=TRUE,baseExpr=NULL)
{
    checkInt<-c(nTrain,nGr1,nBiom, nRep)
    if(!identical(checkInt, floor(checkInt))) stop("non-integer(s) given where argument(s) should be integer valued") 

    if (!is.null(baseExpr) && length(baseExpr)!=nBiom)
    {stop("length of 'baseExpr' should be equal to nBiom") }

    if(nTrain<2) stop("training set size can not be less than 2")

    if( any(c(nGr1, nBiom, nRep)<1)) stop("'nGr1', 'nBiom' and 'nRep' can not be less than 1") 
    nGr2<-nTrain-nGr1

    if (nGr2<1) stop("'nTrain' must be greater than 'nGr1'")
    if(!is.null(rhoMax) && (rhoMax<0 || rhoMax>0.95))stop("allowed values of 'rhoMax' are between 0 and 0.95 inclusive")
    
     if (!is.null(rhoMin) && (rhoMin<0 || rhoMin>0.95))stop("allowed values of 'rhoMin' are between 0 and 0.95 inclusive")
     if (!is.null(rhoMin) && !is.null(rhoMax) && (rhoMin>rhoMax))stop("'rhoMin' should be less than or equal to 'rhoMax'")
    if (!is.null(bSizes) && sum(bSizes)!=nBiom)stop("block sizes do not add up to the nBiom")
    if (!is.null(gamma) && gamma<0) stop("'gamma' should be non-negative.")

    ## Simulate  nTrain by nBiom by nRep array of experimental errors
    ## --------------------------------------------------------------
    if (nRep>1) e<-array(rnorm(nTrain*nBiom*nRep,mean=0,sd=sdW),dim=c(nTrain,nBiom,nRep)) else
                e<-array(rnorm(nTrain*nBiom,mean=0,sd=sdW),dim=c(nTrain,nBiom))
    
    ## Average e over replicates
    if (nRep>1) eAvg<-apply(e,c(1,2),mean) else eAvg<-e
    
    ## Simulate  nTrain by nBiom matrix  of biological errors
    ## ------------------------------------------------------
    b<-matrix(rnorm(nTrain*nBiom, mean=0,sd=sdB),nTrain)
    
    ## Add biological and experimental errors
    ## ---------------------------------------
    randErr<-b+eAvg

    ## Introduce correlation structure
    ## -------------------------------

    ## Block diagonal (Hub-Toeplitz) or Exchangeable correlation structure
   
 if (nBiom>1 && !is.null(gamma)) {
    ## Define exchangeable correlation matrix
    if (gamma==0) {
        repeat{
        if (is.null(rhoMax)) rhoMax<-runif(1,min=0.6,max=0.8)
        if (is.null(rhoMin)) rhoMin<-runif(1,min=0.2,max=0.4)
        if (rhoMin>rhoMax)stop("'rhoMin' should be less than or equal to 'rhoMax'")
        Rho<-diag(rep(1,nBiom))
        Rho[upper.tri(Rho)|lower.tri(Rho)]<-0.5*(rhoMax+rhoMin)
        if(all(eigen(Rho,only.values=TRUE)[[1]]>0)) break
    }
    }
    ## Block diagonal (Hub-Toeplitz) correlation structure
    if (gamma>0){
        ## Define block number and block sizes
        if (nBiom<5)bMax<-1 else bMax<-floor(nBiom/bsMin)
        if(is.null(nBlock)) nBlock<-sample(1:bMax,1)
        if(is.null(bSizes)){
            bs<-nBiom%/%nBlock
            mod<-nBiom%%nBlock
            bs1<-bs+mod
            bSizes<-c(bs1,rep(bs,nBlock-1))
        }
        lRho<-vector("list", nBlock)
        lRho
        for (i in 1:length(bSizes)){
            lRho[[i]]<-diag(rep(1, bSizes[i]))
        }

         fill<-function(x,rhoMax, rhoMin, gamma){
           
            k<-ncol(x)
            y<-x[1,]
            y[2:k]<-rhoMax-(((2:k)-2)/(k-2))^gamma*(rhoMax-rhoMin)
            for (i in 1:(k-1)){
                x[i,(i+1):k]<-y[2:(k+1-i)]
            }
            x[lower.tri(x)]<-t(x)[lower.tri(t(x))]
           return(x)                	
        }
        
        if (nullMax<-is.null(rhoMax)) rhoMax<-runif(1,min=0.6,max=0.8)
        if (nullMin<-is.null(rhoMin)) rhoMin<-runif(1,min=0.2,max=0.4)

        lRho<-lapply(lRho, fill, rhoMax, rhoMin,gamma)
        if (nullMax || nullMin){
            repeat{
                posd<-lapply(lRho, function(x) eigen(x, only.values=TRUE)$values>0)
                if (all(unlist(posd))) break
                rhoMax<-runif(1,min=0.6,max=0.8)
                rhoMin<-runif(1,min=0.2,max=0.4)
                lRho<-lapply(lRho, fill, rhoMax, rhoMin,gamma) 
            }
            
        }
        
        Rho<-as.matrix(Matrix::bdiag(lRho))
    }
      
      ## Existing standard deviations in the data
      
      sdMat<-diag(sqrt(diag(var(randErr))))

      ## Scale the data to have unit variances
      randErrScaled<-scale(randErr,center=FALSE,scale=TRUE)
   
      ## Covariance matrix with desired Rho


      
      covMat<-sdMat%*%Rho%*%sdMat
      cholRoot<-chol(covMat) ## property: t(cholRoot)%*%cholRoot==covMat

      ## Transformed data to have desired covariance structure

      randErrTrans<-randErrScaled%*%cholRoot

      ## Replace original data by the transformed data

      randErr<- randErrTrans
}

    ## Give column names
    ## -----------------
    colnames(randErr)<-paste("Biomarker",1:nBiom, sep="")
    
    ## Add \mu_g (baseline expression) to each row of randErr
    ## ------------------------------------------------------
    
    if (is.null(baseExpr))G<-rep(5,nBiom) else G<-baseExpr
    avgData<-randErr+matrix(G, nrow=nTrain, ncol=nBiom, byrow=T)

    
    ## Introducing differential expressions
    ##-------------------------------------
     
    if(diffExpr){
    ## Generate indicators for up/down regulation
    upDown<-sample(c(-1,1),nBiom,replace = TRUE)
        
    ## Generate log2(fold-change) from half normal distribution with mean=0, sd=sigma
    
    foldChange<-rtnorm(nBiom,mean=0,sd=sigma,lower=log2(foldMin))

    ## Generate "noChange" for the healthy group from normal distribution with mean=0, sd=sigma
    ## This is just to match the variability of the data in both groups

    noChange<-rtnorm(nBiom,mean=0,sd=sigma,lower=-Inf)
    
    if(orderBiom) foldChange<-sort(foldChange)
    
    ## Make the data in the diseased group (D) to be up/down regulated
    ## by the amount foldChange
      
    diffD<-matrix(foldChange*upDown, ncol=nBiom, nrow=nGr2, byrow=T)

    diffH<-matrix(noChange, ncol=nBiom, nrow=nGr1, byrow=T)
    
    avgData[(nGr1+1):nTrain,]<-avgData[(nGr1+1):nTrain,]+diffD
    avgData[1:nGr1,]<-avgData[1:nGr1,]+diffH
      }

    ## Prepare data for classification
    ##--------------------------------

    ## Define groups (H=Healthy, D=Diseased)

    class<-factor(c(rep("H",nGr1),rep("D",nGr2)))
    data<-data.frame(class=class,avgData)
    return(list(data=data,corrMat=if(is.null(gamma)) NULL else Rho))
}

################################################################################
## End of:                         simData.R
################################################################################
