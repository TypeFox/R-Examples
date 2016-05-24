###################################################################################
#' Return a SFA function as a quadratic form.
#' 
#'    sfaGetHf returns function number NR in the
#'    sfa object referenced by HDL in the form of a quadratic form\cr
#'             q(x) = 1/2*x'*H*x + f'*x + c\cr
#'    Of course, this only works if a quadratic expansion was used during
#'    training.
#'    The quadratic form can lie in different spaces, i.e. it can receive
#'    as input preprocessed or non-preprocessed vectors. This is specified
#'    by setting the argument WHERE. The quadratic form lies
#'     - in the preprocessed space for WHERE==0 (e.g. the whitened space if
#'       the preprocessing type is PCA)
#'     - in the PCA space (i.e. projected on the principal components but
#'       not whitened, works only if PCA was used for preprocessing) for
#'       WHERE==1
#'     - in the input, mean-free space for WHERE==2
#'     - in the input space for WHERE==3
#'    In general you will need to set WHERE to 2 or 3, but working in the
#'    preprocessed spaces can often drastically improve the speed of
#'    analysis.
#'
#' @param sfaList 			A list that contains all information about the handled sfa-structure
#' @param nr				function number
#' @param where				WHERE parameter
#'
#' @return list \code{res} \cr
#' - \code{res} contains:
#'	\code{res$H}
#'	\code{res$f}
#'	\code{res$c}
#'
#' @seealso  \code{\link{sfa2Create}}
#' @keywords internal
###################################################################################
#TODO needs testing, not exported, internal
sfaGetHf <- function(sfaList, nr, where){
	#%%% check arguments
	if(sfaList$deg==1){
		stop("sfaGetHf: sfaList is SFA1-object")
	}
	if(where==1 && sfaList$pp_type=="SFA1"){
		stop("sfaGetHf: sfaList preprocessing type is SFA1")
	}
	if(where>3 || where<0){
		stop("sfaGetHf: wrong where argument.")
	}
	sf=sfaList$SF[nr,];
	c=-sfaList$avg1%*%t(sf);
	
	pcaDim=sfaList$pp_range;
	if(length(pcaDim)>1){pcaDim=pcaDim[2]-pcaDim[1];}
	
	#%--- split linear and quadratic part
	#% f is linear part
	#% H is matrix of the quadratic part
	f=t(sf[1:pcaDim]);  
	H=matrix(0,pcaDim,pcaDim);
	k=pcaDim;
	for(i in 1:pcaDim){
		for(j in 1:pcaDim){
			if(j>i){ 
				k=k+1; 
				H[i,j]=sf[k];
			}else if(j==i){ 
				k=k+1; 
				H[i,j]=2*sf[k];
			}else{H[i,j]=H[j,i];}
		}
	}
	#% transform H and f according to 'where'
	if(where==1){
		D=diag(sfaList$D0^(-0.5));
		H=t(D)%*%H%*%D; 
		f=t(D)%*%f;
	}
	else if(where>=2){
		W0=sfaList$W0;
		H=t(W0)%*%H%*%W0; 
		f=t(W0)%*%f;	
		if(where==3){
			avg = t(sfaList$avg0);
			c = 0.5%*%t(avg)%*%H%*%avg - t(f)%*%avg + c;
			f = -H%*%avg + f; 
		}
	}

	res$H=H
	res$f=f
	res$c=c
	return(res)
} #end of classify
