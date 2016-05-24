#TODO: Nearest Neigbhour und RF einbauen?
###################################################################################
#' Classifier for SFA demos
#' 
#'  Train or apply a Gaussian classifier..
#'
#' @param gauss		List created by gaussCreate. Contains also the elements: \describe{
#'					\item{aligned}{
#'					     =0: do not align the Gaussian classifiers with axes, use full covariance matrix \cr
#'					     =1 (default): set the off-diagonals in covariance matrix to 0,
#' 		       		  i.e. the Gaussian classifier is forced to be aligned with the	
#'            		axes. This is more robust in the case where the data deviate
#'            		largely from a multivariate normal distribution.    }
#'    			\item{epsD}{  
#'					      [defaults to 0.04] replace diagonal elements of COV smaller than epsD with
#'            		epsD to avoid too small Gaussians }
#'          }
# %					(only for algo='gauss')
# %     algo  either 'gauss' (def.) or 'nearneig' (nearest neigbour)
#' @param y    		K x M matrix where K is the total number of patterns and M is the
#'            		number of variables used for classification. I.e. each row of
#'           		y contains the data for one pattern.
#' @param realC  	1 x K matrix with NCLASS distinct real class labels needed only
#'            		for method='train'. In case of method="apply" realC is not 
#'            		used and can have any value
#' @param method 	either "train" (default) or "apply"
#'
#' @return list \code{gauss} containing
#'    \item{gauss$predC}{ 1 x K matrix: the predicted class  }
#'    \item{gauss$prob}{ K x NCLASS matrix: prob(k,n) is the estimated probability that pattern k belongs to class m  }
#'
#' @seealso  \code{\link{gaussCreate}} 
#' @export
#' @importFrom MASS ginv
###################################################################################
gaussClassifier <- function(gauss,y,realC,method="train"){
	aligned=gauss$aligned
	epsD=gauss$epsD
	nclass = gauss$nclass; 
    if(is.vector(y)){y=t(as.matrix(y))}
	else{y=as.matrix(y)};
    if(method=="train"){
        gauss$uni = unique(realC);        #% the unique labels 
        #% (during training we must have at least one record from each
        #% class, during 'apply' there can be fewer class labels present, 
        #% therefore we store gauss$uni for later reference)
		first=1;
        for(n in 1:nclass){
            ind = which(realC==gauss$uni[n]);
            gauss$X0[n,] = colMeans(y[ind,]);
            #arg = y[ind,]-customRep(t(gauss$X0[n,]),length(ind)); 
            arg = y[ind,]-matrix(t(gauss$X0[n,]),length(ind),length(t(gauss$X0[n,])),byrow=TRUE); #MZ, 11.11.12: speedfix
            COV = t(arg) %*% arg/(length(ind)-1);
			if (epsD>0){
                dcov = diag(COV);
                if ((min(dcov)<epsD) && (first==1)){
                    first=0;
                    warning("Warning: Some diagonal elements of COV are smaller than epsD and are replaced with epsD");
                }
                dcov[dcov<epsD] = epsD;
                COV = COV - diag(diag(COV)) + diag(dcov); # replace orig. diagonal with diagonal without small value
            }
            if(aligned == 1){
                #% zero all off-diagonal elements >> ensures that the densitiy
                #% function has ellipsoidal equipotential surfaces and that
                #% the ellipsoid axes are *aligned* with the coordinate axes
                COV = diag(diag(COV));
            }
            gauss$COV[,,n] = COV; 
			iCOV<-try(solve(COV)); 
			if(class(iCOV) == "try-error"){ 
				#warning("The function solve(COV) in gaussClassifier produced an error, ginv(COV) is used instead");   
				iCOV= ginv(COV);
			}
			gauss$iCOV[,,n] = iCOV;
			gauss$f0[n] = 1/(2*pi)^(gauss$odim/2)/sqrt(det(COV));
			gauss$Pc[n] = length(ind)/length(realC);            
        }
        gauss$isTrained=1;
    }
	else if(method=="apply"){}
    else{stop(paste("method = ", method, " is not allowed"));}
	
    if (gauss$isTrained==0){
        stop(paste("Gaussian Classifier is not trained >> cannot apply it"));
    }
    #% gauss
    #% this is done for both methods, 'train' and 'apply':
    #%
    prob = matrix(0,customSize(y,1),nclass);
	predC = matrix(0,customSize(y,1),1);
    for(k in 1:customSize(y,1)){
        sumPy=0;
		Pyc<-matrix(0,1,nclass);
        for(n in 1:nclass){
            arg = y[k,]-gauss$X0[n,];
            Pyc[n] = gauss$f0[n] * exp(-arg %*% gauss$iCOV[,,n]%*%arg/2); #TODO maybe -t(arg) instead of -arg
            sumPy = sumPy+Pyc[n]*gauss$Pc[n];
        }
        if(sumPy==0){
            #% the pattern y(k,:) is so far away from all Gaussians that we
            #% can't tell the class >> random guess
            Pcy = runif(nclass);
            warning(paste("sumPy==0, pattern far outside for k=",k));
        }else{
            Pcy = Pyc * t(gauss$Pc) / sumPy; #?
        }
        predC[k] = gauss$uni[which(Pcy==max(Pcy))];
        prob[k,] = Pcy;
    }
	gauss$predC=predC;
	gauss$prob=prob; #TODO: is prob even asked for anywhere? not yet found
	return(gauss) 
}

###################################################################################
#' Create an Gaussian classifier object
#'
#' @param nclass			number of classes
#' @param dimY				dimension
#'
#' @return list of defaults for gauss classifier
#'
#' @seealso  \code{\link{gaussClassifier}}
#' @export
###################################################################################
gaussCreate <- function(nclass,dimY){  
	gauss=list();
	gauss$nclass=nclass;
	odim=dimY;
	gauss$odim=odim
	gauss$isTrained=0;
	gauss$f0=matrix(0,nclass,1);      # % normalization factor
	gauss$Pc=matrix(0,nclass,1);     # % apriori prob for class c
	gauss$X0=matrix(0,nclass,odim);
	gauss$COV=array(0,c(odim,odim,nclass));
	gauss$iCOV=array(0,c(odim,odim,nclass));
	gauss$aligned=1
	gauss$epsD=0.04
	return(gauss)
}