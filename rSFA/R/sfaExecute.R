###################################################################################
#' Execute learned function for input data 
#' 
#' After completion of the learning phase (step="sfa") this function can be used
#' to apply the learned function to the input data. \cr
#'    The execution is completed in 4 steps:\cr
#'     1. projection on the input principal components (dimensionality
#'     reduction)\cr
#'     2. expansion (if necessary)\cr
#'     3. projection on the whitened (expanded) space\cr
#'     4. projection on the slow functions
#'
#' @param sfaList 			A list that contains all information about the handled sfa-structure
#' @param DATA				Input data, each column a different variable
#' @param prj				If not NULL, the preprocessing step 1 is skipped for SFA2				
#' @param ncomp				number of learned functions to be used
#'
#' @return matrix \code{DATA} containing the calculated output \cr
#'
#' @seealso  \code{\link{sfa2}} \code{\link{sfa1}} \code{\link{sfaStep}}
#' @export
###################################################################################
sfaExecute <- function (sfaList, DATA, prj=NULL, ncomp=NULL){
  if (!is.null(ncomp))
    if (ncomp<1 | ncomp>dim(sfaList$SF)[1])
      stop(sprintf("argument ncomp=%d not in allowed range [1,%d]",ncomp,dim(sfaList$SF)[1]));
      
	if(is.vector(DATA)){DATA=t(as.matrix(DATA))}
	else{DATA=as.matrix(DATA)};
	if (sfaList$deg>=2){
		if(is.null(prj)){#TODO: why is prj not used anywhere else here in this function? whats the use of it?			
			DATA=(DATA-matrix(sfaList$avg0,customSize(DATA,1),length(sfaList$avg0),byrow=TRUE))%*%t(sfaList$W0); #MZ, 11.11.12: speedfix
		}
		DATA=sfaList$sfaExpandFun(sfaList, DATA);
		DATA=DATA-matrix(sfaList$avg1,customSize(DATA,1),length(sfaList$avg1),byrow=TRUE);#MZ, 11.11.12: speedfix
		if(is.null(ncomp)){
			DATA=DATA%*%t(sfaList$SF);
		}
		else{			
			DATA=DATA%*%t(sfaList$SF[1:ncomp,]);
		}
	}
	else{   #deg==1
		DATA=DATA-matrix(sfaList$avg0,customSize(DATA,1),length(sfaList$avg0),byrow=TRUE);#MZ, 11.11.12: speedfix
		if (!is.null(sfaList$SFWt)){
			DATA=DATA%*%sfaList$SFWt;
		}
		else{
			if(is.null(ncomp)){
				DATA=DATA%*%t(sfaList$SF);
			}
			else{			
				DATA=DATA%*%t(sfaList$SF[1:ncomp,]);
			}
		}
	}
	return(DATA)
}



