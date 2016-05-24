###################################################################################
#' Update a step of the SFA algorithm. 
#' 
#'  sfaStep() updates the current step of the SFA algorithm. Depending on \code{sfaList$deg}
#'  it calls either \code{\link{sfa1Step}} or \code{\link{sfa2Step}} to do the main work. 
#'  See further documentation there
#'	
#' @param sfaList 			A list that contains all information about the handled sfa-structure
#' @param arg				Input data, each column a different variable
#' @param step				Specifies the current SFA step.  Must be given in the right sequence: 
#'								for SFA1 objects:  "preprocessing", "sfa"\cr
#'								for SFA2 objects:  "preprocessing", "expansion", "sfa"		
#' 						  	Each time a new step is invoked, the previous one is closed, which
#'   						might take some time.				
#' @param method			Method to be used: For \code{sfaList$step="expansion"} the choices are "TIMESERIES" or "CLASSIF". \cr
#'							For \code{sfaList$step="sfa"} (\code{\link{sfa2Step}} only) the choices are "SVDSFA" (recommended) or "GENEIG" (unstable).
#'
#' @return list \code{sfaList} taken from the input, with new information added to this list. 
#'    See \code{\link{sfa1Step}} or \code{\link{sfa2Step}} for details.
#'
#' @examples
#'    ## Suppose you have divided your training data into two chunks,
#'    ## DATA1 and DATA2. Let the number of input dimensions be N. To apply
#'    ## SFA on them write:
#'    \dontrun{ 
#'    sfaList = sfa2Create(N,xpDim(N))
#'    sfaList = sfaStep(sfaList, DATA1, "preprocessing")
#'    sfaList = sfaStep(sfaList, DATA2)
#'    sfaList = sfaStep(sfaList, DATA1, "expansion")
#'    sfaList = sfaStep(sfaList, DATA2)
#'    sfaList = sfaStep(sfaList, NULL, "sfa")
#'    output1 = sfaExecute(sfaList, DATA1)
#'    output2 = sfaExecute(sfaList, DATA2)
#'    }
#'
#' @seealso \code{\link{sfa1Step}} \code{\link{sfa2Step}}  \code{\link{sfa1Create}} \code{\link{sfa2Create}} \code{\link{sfaExecute}}
#' @export
###################################################################################
sfaStep <- function (sfaList, arg, step=NULL, method=NULL){
	if(!is.null(arg)){
		arg<-as.matrix(arg)
	}
	if (is.null(method)){
		if (!is.null(step) && (step=="sfa")){
			method = "SVDSFA";
		}
		else{
			method = "TIMESERIES";
		}
	}
	if (sfaList$deg==1){
		sfaList<-sfa1Step(sfaList, arg, step, method);}
	else{
		sfaList<-sfa2Step(sfaList, arg, step, method);}
	return(sfaList)
}

###################################################################################
#' A step in the SFA2 algorithm. 
#' 
#'  !!! Do not use this function directly, use sfaStep instead !!!
#'
#' @param sfaList 			A list that contains all information about the handled sfa-structure
#' @param arg				Input data, each column a different variable
#' @param step				Specifies the current SFA step.  Must be given in the right sequence: 
#'								for SFA1 objects:  "preprocessing", "sfa"\cr
#'								for SFA2 objects:  "preprocessing", "expansion", "sfa"		
#' 						  	Each time a new step is invoked, the previous one is closed, which
#'   						might take some time.				
#' @param method			Method to be used: For \code{sfaList$step="expansion"} the choices are "TIMESERIES" or "CLASSIF". \cr
#'							For \code{sfaList$step="sfa"} the choices are "SVDSFA" (recommended) or "GENEIG" (unstable).
#' 						 	GENEIG is not implemented in the current version, since
#' 						 	R lacks the option to calculate generalized eigenvalues easily.
#'
#' @return list \code{sfaList} taken from the input, with new information added to this list. 
#'    Among the new items are:
#'    \item{avg0}{  mean vector in input space}
#'    \item{avg1}{  mean vector in expanded space}
#'    \item{W0}{ (ppRange x ppRange)-matrix, the whitening matrix for the input data}
#'    \item{C}{ covariance matrix of the time-diff of expanded and sphered data}
#'    \item{SF}{  (sfaRange x sfaRange)-matrix with rows which contain the directions in expanded space with slow signals. The rows are 
#'        sorted acc. to increasing eigenvalues of C}
#'
#' @seealso  \code{\link{sfaStep}} \code{\link{sfa2Create}} \code{\link{sfa1Step}}
#' @export
#' @keywords internal
###################################################################################
sfa2Step <- function (sfaList, arg=NULL, step=NULL, method=NULL){
	#if(is.null(sfaList$dbg)){dbg<-0}else{dbg<-sfaList$dbg}
	if(is.null(sfaList$opts$epsC)){epsC<-1e-7}else{epsC<-sfaList$opts$epsC}
	if(!is.null(step))
	{
		oldStep=sfaList$step
		# step: init -> preprocessing
		if (oldStep=="init" & (step=="preprocessing")){
			print("Start preprocessing");
			if (substr(sfaList$ppType, 1, 3)=="PCA"){ # check if first three chars are PCA: PCA, PCA2 or PCAVAR
				sfaList$lcov=lcovCreate(ncol(arg));
				#sfaList$diff=sfaList$lcov;
			}
			else{
				sfaList$sfa1List=sfa1Create(sfaList$ppRange);
			}

		}		
		# step: preprocessing -> expansion
		else if (oldStep=="preprocessing" & (step=="expansion")){
			print("Close preprocessing");
			if(sfaList$ppType=="SFA1"){
				sfaList$sfa1List=sfaStep(sfaList$sfa1List, NULL, "sfa")
				sfaList$W0=sfaList$sfa1List$SF;
				sfaList$D0=sfaList$sfa1List$DSF;
				sfaList$avg0=sfaList$sfa1List$avg0; #save avg and tlen from lcov
				sfaList$tlen0=sfaList$sfa1List$tlen0;
				sfaList$sfa1List=NULL; # clear sfa1List
			}
			else{#use PCA if not SFA1  
				sfaList$lcov=lcovFix(sfaList$lcov)	
				if(sfaList$ppType=="PCA"){ 
					print("Whitening and dimensionality reduction (PCA)");
					pcaResult=lcovPca(sfaList$lcov,sfaList$ppRange)
					sfaList$W0=pcaResult$W;
					sfaList$DW0=pcaResult$DW;
					sfaList$D0=pcaResult$D;
					sfaList$avg0=sfaList$lcov$avg; #save avg and tlen from lcov
					sfaList$tlen0=sfaList$lcov$tlen;
					#additional check: is covariance matrix illconditioned?
					sfaCheckCondition(sfaList$lcov$COVMTX, "input")
				}
				else if(sfaList$ppType=="PCA2"){
					# the improved preprocessing sphering by Konen, using SVD.
					# Redundant dimensions with eigenvalue close to zero are detected
					# and the corresponding rows in W0 removed.
					print("Whitening and dimensionality reduction (PCA2)");
					pcaResult=lcovPca2(sfaList$lcov,sfaList$ppRange)
					sfaList$W0=pcaResult$W;
					sfaList$DW0=pcaResult$DW;
					sfaList$D0=pcaResult$D;
					# lcovPca2 will null the rows of SFA_STRUCTS{hdl}.W0 with too
					# small eigenvalues. Here we reduce the rows of W0 and the numbers
					# pp_range  and xp_range accordingly:
					ppRange=length(which(colSums(t(sfaList$W0))!=0)); 
					sfaList$ppRange=ppRange
					sfaList$xpRange=sfaList$xpDimFun(ppRange) 
					sfaList$sfaRange=min(cbind(sfaList$xpRange,sfaList$sfaRange)); # ??
					sfaList$W0=sfaList$W0[1:ppRange,];
					
					sfaList$avg0=sfaList$lcov$avg;
					sfaList$tlen0=sfaList$lcov$tlen;					
				}
				else if(sfaList$ppType=="PCAVAR"){
					# another preprocessing as done in [Wiskott&Sejnowski'02] which
					# does not use PCA, but simply shifts and scales the input data to
					# have zero mean and unit variance
					#
					print("unit variance w/o dimensionality reduction (PCAVAR)");
					varmat = diag(diag(sfaList$lcov$COVMTX)); 
					sfaList$W0 = varmat^(-0.5);
					sfaList$avg0=sfaList$lcov$avg;
					sfaList$tlen0=sfaList$lcov$tlen;	
				}			
				sfaList$lcov=NULL; # clear lcov
			}
			print("Init expansion step");
			#inSize=sfaList$ppRange; #used nowhere.. why?
			#if (length(inSize)==2){
			#	inSize=inSize[2]-inSize[1]+1; 
			#}
			xpSize=sfaList$xpRange;
			sfaList$xp=lcovCreate(xpSize);
			sfaList$diff=lcovCreate(xpSize);
		}
		# step: expansion -> sfa
		else if (oldStep=="expansion" & (step=="sfa")){
			print("Close expansion step");
			sfaList$xp=lcovFix(sfaList$xp);
			sfaList$avg1=sfaList$xp$avg;
			sfaList$tlen1=sfaList$xp$tlen;
			xpsize=sfaList$xpRange
			sfaList$diff=lcovFix(sfaList$diff);
			print("Perform Slow Feature Analysis")
			sfaInt=sfaGetIntRange(sfaList$sfaRange);
			################################################################################
			#First check method
			if(method=="GENEIG" ){#|| dbg>0){
				stop("GENEIG method is not implemented in rSFA package. 
					Please choose method SVDSFA instead.") #see note below
				#sfaCheckCondition(sfaList$xp$COVMTX, "expanded")
				#
				# Original Code
				#
				#Bm<-1*sfaList$xp$COVMTX
				#Am<-1*sfaList$diff$COVMTX
				#res<-sfaDggev(Am,Bm) #Please note: sfaDggev is not running properly, thus deprecated, not part of package. Code not working.
				#D=res$val;
				#sfaList$SF<-res$vec;
				#
				# End Originial Code
				#
			}
			#TODO WHY IF INSTEAD OF ELSE IF ???
			if(method=="SVDSFA"){			
				# extension /WK/08/2009: first sphere expanded data with
				# LCOV_PCA2, taking care of zero or very small eigenvalues in B
				# by using the SVD approach
				# 
				print("Using alternate [WisSej02] approach for SFA-calculation ..."); #TODO verbosity einführen
				pcaResult<-lcovPca2(sfaList$xp);
				S<-pcaResult$W    # S: sphering matrix for expanded data (xphdl)
				#not used anywhere?#DS<-pcaResult$DW   # DS: de-sphering matrix, BD: eigenvalues  
				BD<-pcaResult$D  # of B (covariance matrix of expanded data)
				C = S %*% sfaList$diff$COVMTX %*% t(S); 
				#res= eigen(C);
				#W1=res$vectors
				#D1=res$values
				resvd=svd(C,nu=0,nv=ncol(C))
				W1=resvd$v;
				D1=resvd$d;
				SF1 = t(S)%*%W1;
				sfaList$SF = SF1;
				sfaList$BD = BD;
				sfaList$myS=S;				
				D=D1;
			}		
			#always calculate rank(B) (for diagnostics only)
			B = sfaList$xp$COVMTX;
			rankB = qr(B)$rank;    #TODO: this might not be completely the same like matlabs rank(B);
			print(paste("rank of B = ",rankB));
			sfaList$rankB = rankB;
			sfaList$myB = B; # needed for nl_regress only
			
			idx=t(order(D));	# % idx(1): index to smallest eigenvalue, idx(2): to 2nd-smallest, ...
			lammax=max(D);
			print(paste("epsC*lammax= ",epsC*lammax)); #TODO maybe remove this print? or only do with high verbosity (implement later)
			
			if(method=="SVDSFA"){ 
				#rankC = qr(C)$rank  #only used in sfatk for a print, skipped.
				#print(paste("rank of C = ",rankC)); #see above
				#idx = idx[which(D[idx]!=0)]; # 'SVDSFA': exclude eigenvalues which
										  # are *exactly* zero from further
										  # analysis, since they correspond to
										  # zeros in the sphering matrix
										  # (degenerate dimensions)
				idx = idx[which(abs(D[idx])>rep(epsC*lammax,length(D[idx])))];    #TODO: ugly solution ?
				sfaInt = 1:length(idx);
				sfaList$sfaRange = length(idx);				
				# REMARK: These statement were also beneficial for 'GENEIG', because
				# there it may happen in the degenerate case that some eigenvalues of
				# D become negative (??, and the corresponding eigenvectors contain
				# only noisy signals). However, we do not apply it here, because it
				# lets 'GENEIG' deviate from the original [Berkes03] code. And the
				# slow signals would still have the wrong variance.
			}
			sfaList$DSF<-t(D[idx[sfaInt]]);
			sfaList$SF<-t(sfaList$SF[,idx[sfaInt]]); 
			################################################################################
			#clear unneeded parts
			sfaList$cp=NULL;
			sfaList$diff=NULL;
			print("SFA2 closed");				
		}
		else if (!(oldStep==step)){ #oldStep and step should only be different for well defined sequences like above
			warning("Unknown Step Sequence in sfa2Step")
			return(sfaList)
		}
		sfaList$step=step;
	}	
	#
	# things to do always when sfaList$step is either 'preprocessing' or 'expansion' 
	# (no matter whether it is invoked for the first time or once again)
	#	
	if(sfaList$step=="preprocessing"){
		if(substr(sfaList$ppType, 1, 3)=="PCA"){
			sfaList$lcov=lcovUpdate(sfaList$lcov,arg);
		}
		else{ #else SFA1
			sfaList$sfa1List=sfaStep(sfaList$sfa1List, arg, "preprocessing")
		}
	}
	if(sfaList$step=="expansion"){
		#arg=arg-customRep(sfaList$avg0,customSize(arg,1)); 
		arg=arg-matrix(sfaList$avg0,customSize(arg,1),length(sfaList$avg0),byrow=T) #MZ, 11.11.12: speedfix
		arg=sfaList$sfaExpandFun(sfaList, arg %*%  t(sfaList$W0));
		sfaList$xp=lcovUpdate(sfaList$xp,arg);		
        if(method=="TIMESERIES"){
            sfaList$diff=lcovUpdate(sfaList$diff, sfaTimediff(arg,sfaList$axType));
        }
		else if (method=="CLASSIF"){
        # extension /WK/08/2009: generate the difference of all pattern
        # pairs in 'pdiff'
        K = customSize(arg,1);
		lt = customSize(arg,2);
			  if(K<2){
				  stop("This class has less than two training records. Expansion can not run, pattern difference can not be calculated")
			  }
            pdiff = NULL;
            for (k in 1:(K-1)){ #TODO: check and maybe improve
                #pdiff = rbind(pdiff, customRep(t(arg[k,]),K-k) - arg[(k+1):K,]);
                pdiff = rbind(pdiff, matrix(t(arg[k,]),K-k,lt,byrow=TRUE) - arg[(k+1):K,]);#MZ, 11.11.12: speedfix
                if (k%%100==0) {       # Time and Mem optimization: do not let pdiff grow too large /WK/01/2012
                  sfaList$diff=lcovUpdate(sfaList$diff, pdiff);
                  pdiff=NULL;
                  #cat("zeroing pdiff\n");
                }
                #cat(k,"\n");flush.console();
            }
            sfaList$diff=lcovUpdate(sfaList$diff, pdiff);
		    }
        else{
		      warning(paste(method," is not an allowed method in expansion step"));
        }		
	}	
	return(sfaList)
}

###################################################################################
#' A step in the SFA1 algorithm. 
#' 
#'  !!! Do not use this function directly, use sfaStep instead !!!
#'
#' @param sfaList 			A list that contains all information about the handled sfa-structure
#' @param arg				Input data, each column a different variable
#' @param step				Specifies the current SFA step.  Must be given in the right sequence: 
#'								for SFA1 objects:  "preprocessing", "sfa"\cr
#'								for SFA2 objects:  "preprocessing", "expansion", "sfa"		
#' 						  	Each time a new step is invoked, the previous one is closed, which
#'   						might take some time.				
#' @param method			Method to be used: For \code{sfaList$step="expansion"} the choices are "TIMESERIES" or "CLASSIF". \cr
#'							For \code{sfaList$step="sfa"} currently no choices.
#'
#' @return list \code{sfaList} taken from the input, with new information added to this list. 
#'    Among the new items are:
#'    \item{avg0}{  mean vector in input space}
#'    \item{SF}{  (sfaRange x sfaRange)-matrix with rows which contain the directions in expanded space with slow signals. The rows are 
#'        sorted acc. to increasing eigenvalues of time-diff covariance matrix}
#'
#' @seealso  \code{\link{sfaStep}} \code{\link{sfa1Create}} \code{\link{sfa2Step}}
#' @export
#' @keywords internal
###################################################################################
sfa1Step <- function (sfaList, arg=NULL, step=NULL, method=NULL){
	
	if(!is.null(step))
	{
		oldStep=sfaList$step
		if (oldStep=="init" & (step=="preprocessing")){
			print("Start preprocessing");
			sfaList$lcov=lcovCreate(ncol(arg));
			sfaList$diff=sfaList$lcov;
		}
		else if (oldStep=="preprocessing" & (step=="sfa")){
			print("Close preprocessing");
			sfaList$lcov=lcovFix(sfaList$lcov);
			sfaList$avg0=sfaList$lcov$avg;
			sfaList$tlen0=sfaList$lcov$tlen;
			print("Perform slow feature analysis");
			
			if(length(sfaList$sfaRange)==1){
				sfaInt=1:sfaList$sfaRange;
			}
			else{ 
				sfaInt=sfaList$sfaRange[1]:sfaList$sfaRange[2];
			}	
			################################################################################
			#
			# original code: with generalized eigenvalues. unstable and bad implementation
			#
			#Bm<-1*sfaList$lcov$COVMTX
			#Am<-1*sfaList$diff$COVMTX
			#res<-sfaDggev(Am,Bm) #will not work for complex inputs	
			#D=res$val; #CAREFULL: This only works for non complex outputs, complex outputs are difficult
			#idx=t(order(D));			
			#sfaList$DSF<-t(D[idx[sfaInt]]); 
			#sfaList$SF<-t(res$vec[,idx[sfaInt]]); 
			#
			# end of original code
			#
				# svd approach instead:
				pcaResult<-lcovPca2(sfaList$lcov);
				S<-pcaResult$W    # S: sphering matrix for expanded data
				C = S %*% sfaList$diff$COVMTX %*% t(S); 
				resvd=svd(C,nu=0,nv=ncol(C))
				W1=resvd$v;
				D=resvd$d;
				sfaList$SF = t(S)%*%W1;
				idx=t(order(D));	# % idx(1): index to smallest eigenvalue, idx(2): to 2nd-smallest, ...
				lammax=max(D);
				if(is.null(sfaList$opts$epsC)){epsC<-0}else{epsC<-sfaList$opts$epsC}
				print(paste("epsC*lammax= ",epsC*lammax)); #TODO maybe remove this print? or only do with high verbosity (implement later)
				idx = idx[which( abs(D[idx])>rep(epsC*lammax,length(D[idx])))]; 
				sfaInt = 1:length(idx);
				sfaList$DSF<-t(D[idx[sfaInt]]); 
				sfaList$SF<-t(sfaList$SF[,idx[sfaInt]]); 
			################################################################################
			#clean up
			sfaList$lcov=NULL;
			sfaList$diff=NULL;
			print("SFA1 closed");	
		}
		else if (!(oldStep==step)){
			warning("Unknown Step Sequence in sfa1Step")
			return(sfaList)
		}		
		sfaList$step=step;
	}	
	if(sfaList$step=="preprocessing"){
		sfaList$lcov=lcovUpdate(sfaList$lcov,arg);
		if(method=="TIMESERIES"){
           sfaList$diff=lcovUpdate(sfaList$diff, sfaTimediff(arg,sfaList$axType));
    }
		else if (method=="CLASSIF"){
            #% extension /WK/12/2009: generate the difference of all pattern
            #% pairs in 'pdiff'
            K = customSize(arg,1);
			lt = customSize(arg,2);
            pdiff = NULL;
            for (k in 1:(K-1)){ #TODO: check and maybe improve
                #pdiff = rbind(pdiff, customRep(arg[k,],K-k) - arg[(k+1):K,]); 
				pdiff = rbind(pdiff, matrix(t(arg[k,]),K-k,lt,byrow=TRUE) - arg[(k+1):K,]);#MZ, 11.11.12: speedfix
            }
            sfaList$diff=lcovUpdate(sfaList$diff, pdiff);
		}
        else{
		    stop(paste(method," is not an allowed method in expansion step"));
        }		
	}
	return(sfaList)
}