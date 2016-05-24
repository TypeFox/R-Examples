###################################################################################
#' Add noisy copies for parametric bootstrap 
#' 
#' Given training data X with true labels REALCLASS, 
#' add new records to X and REALCLASS, which are noisy copies of the training data. 
#'
#' @param x 			a matrix containing the training data
#' @param realclass		true class of training data (can be vector, numerics, integers, factors)
#' @param pars			list of parameters:\cr
#'						\code{pars$ncopies}: Number of new records to add\cr
#'						\code{pars$ncsort}: Defines if training data should be sorted by class. Default is FALSE\cr
#' 						\code{pars$ncsigma}: The noise in each column of x has the std.dev. pars$ncsigma*(standard deviation of column). Default Value: 0.8\cr
#'						\code{pars$ncmethod}:  =1: each 'old' record from X in turn is the centroid for a new pattern;\cr
#'                    					=2: the centroid is the average of all records from the same class, the std.dev. is the same for all classes;\cr
#'                    					=3: centroid as in '2', the std.dev. is the std.dev. of all records from the same class  (*recommended*)\cr
#'
#' @return list \code{res} \cr
#' - \code{res} contains two list entries: realclass and x (including added copies)
#'
#' @references  \code{\link{sfaPBootstrap}}
#' @export
###################################################################################
addNoisyCopies <- function(realclass, x, pars) {
    ncopies = pars$ncopies;
    if(is.null(pars$ncsigma)) pars$ncsigma=1.0;                           
    if(is.null(pars$ncmethod)) pars$ncmethod=3; 
    if(is.null(pars$ncsort)) pars$ncsort=FALSE;                           

    if (ncopies>0){
        #% add noisy copies
        x0 = x;
        r0 = realclass;
        c0 = getCentroid(x0,r0,pars$ncmethod);
        s0 = getSigma(x0,r0,pars$ncmethod);

        sz1 = customSize(x0,1);           #% number of training records in x0
        sz2 = customSize(x0,2);
        nn1 = trunc(ncopies/sz1)+1;  # % example: ncopies=220, sz1=100 --> nn1=3
        #% We divide the generation of the ncopies bootstrap patterns
        #% into nn1 loop steps. This division is only needed pars.ncmethod=1
        #% (each training record in turn serves as centroid). If
        #% pars.ncmethod=2, we could as well generate the bootstrap patterns
        #% in one step.
        for (k in 1:nn1){
            x2 = c0 + pars$ncsigma * s0 * matrix(rnorm(sz1*sz2),sz1,sz2);  # /WK/ runif -> rnorm
            #% i.e., calculate for each column of x0 the standard deviation,
            #% repeat it sz1 times, multiply each entry by a random number
            #% with standard normal distribution, multiply by pars.ncsigma
            #% --> as a result we have
            #%       sd(x2-c0) = pars.ncsigma * sd(x0)
            #% (approximately, up to statistical fluctuations).
            #% The tacit assumption is here that the standard deviations in
            #% each column are the same for all classes. 
            
            if (k==nn1){
                #% in the last pass, add only that many randomly selected
                #% records from x2 that the total number of added records 
                #% amounts exactly to pars.ncopies
                idx = sample(sz1);
                idx = idx[1:(ncopies%%sz1)];
                x2 = x2[idx,];
                r0 = r0[idx];
            }
            x = rbind(x, x2);
            realclass = c(realclass, r0);
        }
        if (pars$ncsort) {                                  # /WK/
          #% sort the training data class-by-class
          idx = order(realclass);
          realclass = realclass[idx];
          x = x[idx,];
        }
    }
	res<-list(realclass=realclass, x=x)
	return(res)

}
  
###################################################################################
# Get centroid for adding noisy copies
#
# @param x0 			training data
# @param r0				true class of training data
# @param ncmethod		possible settings:\cr
#						=1: c0 = x0 (each pattern is its own centroid). Faster, but not exactly 'parametric bootstrap'\cr
#                    	=2: c0 = centroid(m(x0)) where m(x0) is the class to which x0 belongs. We recommend this method since
#						it adds new patterns which are statistically independent from x0 (true parametric bootstrap)
#
# @return matrix \code{c0} \cr
# - \code{c0} are the centroids for adding noisy copies
#
# @references  \code{\link{sfaPBootstrap}} \code{\link{addNoisyCopies}}
###################################################################################
#NOT EXPORTED, NOT DOCUMENTED
getCentroid <- function(x0,r0,ncmethod) {
    if (ncmethod==1){c0=x0;} 
    else{
        classes=unique(r0);
        c0 = x0*0; 
        for(i in 1:length(classes)){
            ind = which(r0==classes[i]);
            c0[ind,] = customRepmat (colMeans(x0[ind,]), length(ind), 1); # /WK/ col mean # /MZ/ colMeans (faster)
        }
    }
	return(c0)
}

###################################################################################
# Get sigma for adding noisy copies
#
# @param x0 			training data
# @param r0				true class of training data
# @param ncmethod		possible settings:\cr
#						<3: s0 = std(x0) (each column has a common sigma, std over all classes) faster, but perhaps not precise modeling\cr
#                    	=3: s0 = std(m(x0)) where m(x0) is the class to which x0 belongs
#
# @return matrix \code{s0} \cr
# - \code{s0} sigma values
#
# @references  \code{\link{sfaPBootstrap}} \code{\link{addNoisyCopies}}
###################################################################################
#NOT EXPORTED, NOT DOCUMENTED
getSigma <- function(x0,r0,ncmethod) {
    sz1 = customSize(x0,1);           #% number of training records in x0
    if (ncmethod<3){ 
        s0=customRepmat (sd(x0),sz1,1); 
    }else{
        classes=unique(r0);
        s0 = x0*0;
        for(i in 1:length(classes)){
            ind = which(r0==classes[i]);
            s0[ind,] = customRepmat (apply(x0[ind,],2,sd), length(ind), 1);  # /WK/ col std  # /MZ/ sd uses cols by default, faster than apply  #MZ: changed back to apply, slower but sd(matrix) is deprecated in R14.X
        }
    }
	return(s0)
}

###################################################################################
#' Parametric Bootstrap
#'
#' If training set too small, augment it with parametric bootstrap
#'
#' @param realclass		true class of training data (can be vector, numerics, integers, factors)
#' @param x 			matrix containing the training data
#' @param sfaList		list with several parameter settings, e.g. as created by \code{\link{sfa2Create}}\cr
#'						\code{sfaList$xpDimFun} (=xpDim by default) calculated dimension of expaned SFA space\cr
#'						\code{sfaList$deg} degree of expansion (should not be 1, not implemented)\cr
#'						\code{sfaList$ppRange} ppRange for SFA algorithm\cr
#'						\code{sfaList$nclass} number of unique classes\cr
#'						\code{sfaList$doPB} do (1) or do no (0) param. bootstrap.
#'
#' @return a list \code{list} containing:
#'  \item{x}{ training set extended to minimu number of recors1.5*(xpdim+nclass), if necessary }
#'  \item{realclass}{ training class labels, extended analogously  }
#'
#' @seealso  \code{\link{addNoisyCopies}}
#' @export
###################################################################################
sfaPBootstrap <- function(realclass,x,sfaList){
	if(is.null(sfaList$xpDimFun)){ #set a default, but should allready be set in sfaClassify
		sfaList$xpDimFun=xpDim
	}
	
    if (sfaList$deg==1){ warning('sfaPBootstrap not yet applicable for sfaList$deg==1 --> nothing done');}          
	else{
		NTrig = ceiling(1.2*(sfaList$xpDimFun(sfaList$ppRange)+sfaList$nclass));
		NDesi = ceiling(1.5*(sfaList$xpDimFun(sfaList$ppRange)+sfaList$nclass));
		N = customSize(x,1);
		if (N<=NTrig){
			#warning("# of training patterns = %d is less than the necessary 1.2*(D_x+K) = %d',N,NTrig));        
			if (sfaList$doPB==0){
				warning("Warning: NOTHING done, because sfaList$doPB==0");
			}else {
				#disp(sprintf('--> augmenting the training set with Parametric Bootstrap to 1.5*(D_x+K) = %d',NDesi));
				pars=list()
				pars$ncopies=NDesi-N;
				pars$ncsigma=1.0;
				pars$ncmethod=3;
				res = addNoisyCopies(realclass,x,pars);
				realclass=res$realclass
				x=res$x
			}
		}
	}
	return(list(realclass=realclass,x=x))
}



