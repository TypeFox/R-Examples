### Documentation ####
#' Linear predictive models estimation with Online Adaptation of Coordinate 
#' Frequencies based on the LIBLINEAR C/C++ Library.
#' 
#' @description
#' \code{LiblineaR.ACF} is a modification of the LiblineaR package that
#' uses the idea of adaptive coordinate frequencies (ACF) method.
#' Solving the linear SVM problem with coordinate descent
#' is very efficient and is implemented in one of the most often used packages,
#' LIBLINEAR (available at http://www.csie.ntu.edu.tw/~cjlin/liblinear). 
#' It has been shown that the uniform selection of coordinates can be 
#' accelerated by using an online adaptation of coordinate frequencies (ACF).
#' This package implements ACF and is based on LIBLINEAR as well as
#' the LiblineaR package (https://cran.r-project.org/package=LiblineaR).
#' It currently supports L2-regularized L1-loss as well as L2-loss linear SVM. 
#' Similar to LIBLINEAR multi-class classification (one-vs-the rest, and 
#' Crammer & Singer method) and cross validation for model selection is 
#' supported. The training of the models based on ACF is much faster than 
#' standard LIBLINEAR on many problems.
#' 
#' @details
#' For details for the implementation of LIBLINEAR, see the README file of the
#' original c/c++ LIBLINEAR library at
#' \url{http://www.csie.ntu.edu.tw/~cjlin/liblinear}.
#' The ACF code can be found at
#' \url{http://www.ini.rub.de/PEOPLE/glasmtbl/code/acf-cd}.
#' 
#' @param data 	a nxp data matrix. Each row stands for an example (sample,
#'   point) and each column stands for a dimension (feature, variable). A sparse
#'   matrix (from SparseM package) will also work.
#' @param target a response vector for prediction tasks with one value for 
#'   each of the n rows of \code{data}. For classification, the values 
#'   correspond to class labels and can be a 1xn matrix, a simple vector or a 
#'   factor. 
#' @param type \code{LiblineaR} can produce several types of (generalized) linear 
#'   models, by combining several types of loss functions and regularization 
#'   schemes. The regularization is L2, and the losses can be the regular L2-loss 
#'   or L1-loss. The default value for \code{type} is 1. Valid options are:
#'   \describe{
#'     \item{for multi-class classification}{
#'       \itemize{
#'         \item 1 -- L2-regularized L2-loss support vector classification (dual)
#'         \item 3 -- L2-regularized L1-loss support vector classification (dual)
#'         \item 4 -- support vector classification by Crammer and Singer
#'       }
#'     }
#'   }
#' 
#' @param cost cost of constraints violation (default: 1). Rules the trade-off
#'   between regularization and correct classification on \code{data}. It can be
#'   seen as the inverse of a regularization constant. See information on the
#'   'C' constant in details below. A usually good baseline heuristics to tune
#'   this constant is provided by the \code{heuristicC} function in the LiblineaR
#'   package.
#' @param epsilon set tolerance of termination criterion for optimization.
#'   If \code{NULL}, the LIBLINEAR defaults are used, which are:
#'   \describe{
#'     \item{if \code{type} is 1, 3 or 4}{\code{epsilon}=0.1}
#'   }
#'   
#'   The meaning of \code{epsilon} is as follows:
#'            Dual maximal violation \eqn{\le \code{epsilon}}{\le epsilon} 
#'            (default 0.1)
#' @param bias if \code{bias} is \code{TRUE} (default), instances of \code{data} becomes [\code{data}; 1].
#' @param wi a named vector of weights for the different classes, used for
#'   asymmetric class sizes. Not all factor levels have to be supplied (default
#'   weight: 1). All components have to be named according to the corresponding
#'   class label. 
#' @param cross if an integer value k>0 is specified, a k-fold cross validation
#'   on \code{data} is performed to assess the quality of the model via a
#'   measure of the accuracy. Note that this metric might not be appropriate if
#'   classes are largely unbalanced. Default is 0.
#' @param change_rate learning rate of the preference adaptation, default is 0.2
#' @param pref_min lower bound on the preference adaptation, default is 1/20
#' @param pref_max upper bound on the preference adaptation, default is 20
#' @param max_iter the maximum number of iterations, default (from original LIBLINEAR code) is 1000.
#' @param verbose if \code{TRUE}, information are printed. Default is
#'   \code{FALSE}.
#' @param ... for backwards compatibility, parameter \code{labels} may be
#'   provided instead of \code{target}. A warning will then be issued, or an
#'   error if both are present. Other extra parameters are ignored.
#' 
#' @return
#' 	If \code{cross}>0, the average accuracy (classification) computed over \code{cross} runs of cross-validation is returned.\cr\cr
#' Otherwise, an object of class \code{"LiblineaR"} containing the fitted model is returned, including:
#' \item{TypeDetail}{A string decsribing the type of model fitted, as determined by \code{type}.}
#' \item{Type}{An integer corresponding to \code{type}.}
#' \item{W}{A matrix with the model weights. If \code{bias} is TRUE, \code{W} contains p+1 columns, the last being the bias term. The columns are named according to the names of \code{data}, if provided, or \code{"Wx"} where \code{"x"} ranges from 1 to the number of dimensions. The bias term is named \code{"Bias"}.If the number of classes is 2, the matrix only has one row. If the number of classes is k>2 (classification), it has k rows. Each row i corresponds then to a linear model discriminating between class i and all the other classes. If there are more than 2 classes, rows are named according to the class i which is opposed to the other classes.}
#' \item{Bias}{TRUE or FALSE, according to the value of \code{bias}}
#' \item{ClassNames}{A vector containing the class names.}
#' 
#' @references
#' 	\itemize{
#' \item 
#' For more information on LIBLINEAR itself, refer to:\cr
#' R.-E. Fan, K.-W. Chang, C.-J. Hsieh, X.-R. Wang, and C.-J. Lin.\cr
#' \emph{LIBLINEAR: A Library for Large Linear Classification,}\cr
#' Journal of Machine Learning Research 9(2008), 1871-1874.\cr
#' \url{http://www.csie.ntu.edu.tw/~cjlin/liblinear}
#' }
#' 
#' @author Aydin Demircioglu \email{aydin.demircioglu@@ini.rub.de}
#'   Based on LiblineaR package by Thibault Helleputte \email{thibault.helleputte@@dnalytics.com} and\cr
#'   Pierre Gramme \email{pierre.gramme@@dnalytics.com}.\cr 
#'   Based on C/C++-code by Chih-Chung Chang and Chih-Jen Lin
#'   Based on C/C++-code by Tobias Glasmachers and Urun Dogan
#' 
#' @note Classification models usually perform better if each dimension of the data is first centered and scaled.
#' 
#' @seealso \code{\link{predict.LiblineaR.ACF}}, \code{\link[LiblineaR]{heuristicC}}
#' 
#' @examples
#' data(iris)
#' attach(iris)
#' 
#' x=iris[,1:4]
#' y=factor(iris[,5])
#' train=sample(1:dim(iris)[1],100)
#' 
#' xTrain=x[train,]
#' xTest=x[-train,]
#' yTrain=y[train]
#' yTest=y[-train]
#' 
#' # Center and scale data
#' s=scale(xTrain,center=TRUE,scale=TRUE)
#' 
#' # Find the best model with the best cost parameter via 3-fold cross-validations
#' tryTypes=c(1,3,4)
#' tryCosts=c(1000,1,0.001)
#' bestCost=NA
#' bestAcc=0
#' bestType=NA
#' 
#' for(ty in tryTypes){
#' 	for(co in tryCosts){
#' 		acc=LiblineaR.ACF(data=s,target=yTrain,type=ty,cost=co,
#'				bias=TRUE,cross=3,verbose=FALSE)
#' 		cat("Results for C=",co," : ",acc," accuracy.\n",sep="")
#' 		if(acc>bestAcc){
#' 			bestCost=co
#' 			bestAcc=acc
#' 			bestType=ty
#' 		}
#' 	}
#' }
#' 
#' cat("Best model type is:",bestType,"\n")
#' cat("Best cost is:",bestCost,"\n")
#' cat("Best accuracy is:",bestAcc,"\n")
#' 
#' # Re-train best model with best cost value.
#' m=LiblineaR.ACF(data=s,target=yTrain,type=bestType,cost=bestCost,bias=TRUE,verbose=FALSE)
#' 
#' # Scale the test data
#' s2=scale(xTest,attr(s,"scaled:center"),attr(s,"scaled:scale"))
#' 
#' # Make prediction
#' pr=FALSE
#' if(bestType==0 || bestType==7) pr=TRUE
#' 
#' p=predict(m,s2,proba=pr,decisionValues=TRUE)
#' 
#' # Display confusion matrix
#' res=table(p$predictions,yTest)
#' print(res)
#' 
#' # Compute Balanced Classification Rate
#' BCR=mean(c(res[1,1]/sum(res[,1]),res[2,2]/sum(res[,2]),res[3,3]/sum(res[,3])))
#' print(BCR)
#' 
#' #' #############################################
#' 
#' # Example of the use of a sparse matrix:
#' 
#' if(require(SparseM)){
#' 
#'  # Sparsifying the iris dataset:
#'  iS=apply(iris[,1:4],2,function(a){a[a<quantile(a,probs=c(0.25))]=0;return(a)})
#'  irisSparse<-as.matrix.csr(iS)
#' 
#'  # Applying a similar methodology as above:
#'  xTrain=irisSparse[train,]
#'  xTest=irisSparse[-train,]
#' 
#'  # Re-train best model with best cost value.
#'  m=LiblineaR.ACF(data=xTrain,target=yTrain,type=bestType,cost=bestCost,bias=TRUE,verbose=FALSE)
#' 
#'  # Make prediction
#'  p=predict(m,xTest,proba=pr,decisionValues=TRUE)
#' 
#'  # Display confusion matrix
#'  res=table(p$predictions,yTest)
#'  print(res)
#' }
#' 
#'
#' 
#' 
#' @keywords classification multivariate models optimize classes
#'
#' @export

### Implementation ####
LiblineaR.ACF<-function(data, target, type=0, cost=1, epsilon=0.01, bias=TRUE, wi=NULL, cross=0, 
    change_rate = 0.2, pref_min = 0.05, pref_max = 20.0, max_iter = 1000, verbose=FALSE, ...) {
	# <Arg preparation>
  
  if(sparse <- inherits(data, "matrix.csr")){
    if(requireNamespace("SparseM",quietly=TRUE)){
  		# trying to handle the sparse martix case
	  	data = SparseM::t(SparseM::t(data)) # make sure column index are sorted
		  n = data@dimension[1]
		  p = data@dimension[2]
    } else {
      stop("newx inherits from 'matrix.csr', but 'SparseM' package is not available. Cannot proceed further. Use non-sparse matrix or install SparseM.")
    }
	} else {
		# Nb samples
		n=dim(data)[1]
		# Nb features
		p=dim(data)[2]
	}

	# Backwards compatibility after renaming 'labels' to 'target'
	cc = match.call()
	if ("labels" %in% names(cc)) {
		if("target" %in% names(cc)) {
			stop("Argument 'labels' is deprecated and was renamed to 'target'. Stopping due to ambiguous call.")
		} else {
			warning("Argument 'labels' is deprecated and was renamed to 'target'.")
			target=eval(cc$labels, parent.frame())
		}
	}
	rm(cc)
	
	# Bias
	b = if(bias) 1 else -1
	
	# Type 
	typesLabels = c("L2-regularized logistic regression (primal)", "L2-regularized L2-loss support vector classification (dual)", "L2-regularized L2-loss support vector classification (primal)", "L2-regularized L1-loss support vector classification (dual)",
									"support vector classification by Crammer and Singer", "L1-regularized L2-loss support vector classification", "L1-regularized logistic regression", "L2-regularized logistic regression (dual)",
									"", "", "",
									"L2-regularized L2-loss support vector regression (primal)", "L2-regularized L2-loss support vector regression (dual)", "L2-regularized L1-loss support vector regression (dual)")
	typesLabels = gsub("[()]","",typesLabels)
	typesCodes = c("", "L2R_L2LOSS_SVC_DUAL", "", "L2R_L1LOSS_SVC_DUAL",
								 "MCSVM_CS", "", "", "",
								 "", "", "",
								 "", "", "")
	types=ifelse(typesCodes=="", "", paste0(typesLabels, " (", typesCodes, ")"))
	if(!type %in% (which(types!="")-1))
		stop("Unknown model type ",type,". Expecting one of: ", paste(which(types!="")-1, collapse=", "))
	isRegression = type>=11
	if (isRegression) 
        stop("Internal error. Regression was not allowed. Please contact the maintainer.")
	
	# Epsilon
	if(is.null(epsilon) || epsilon<0){
		# Will use LIBLINEAR default value for epsilon
		epsilon = -1
	}
	
	# Target
	if(!is.null(dim(target))) {
		if(identical(dim(target), c(n,1)))
			target = target[,1]
		else
			stop("Wrong dimension for target")
	}
	if(length(target)!=n){
		stop("Number of target elements disagrees with number of data instances.")
	}
	
    if(is.character(target))
        target = factor(target)
    
    # targetLabels are sorted by first occurrence in target ; if target is a factor, targetLabels too, with the same levels.
    targetLabels = unique(target)
    nbClass = length(targetLabels)
    
    yC = as.integer(target)
    
    if (nbClass<2)
        stop("Wrong number of classes ( < 2 ).")
    isMulticlass = (nbClass>2) || (nbClass==2 && type==4)
    
    # Different class penalties? Default to 1
    nrWi=nbClass
    Wi=rep(1,times=nrWi)
    names(Wi)=as.character(targetLabels)
    WiLabels=as.integer(targetLabels)
    
    if(!is.null(wi)) {
        if(is.null(names(wi)))
            stop("wi has to be a named vector!")
        
        if( !all(names(Wi)%in% names(wi)) )
            stop("Mismatch between provided names for 'wi' and class labels.")
        
        for(i in 1:length(wi)){
            Wi[names(wi)[i]]=wi[i]
        }
    }
	
	# Cross-validation?
	if(cross<0){
		stop("Cross-validation argument 'cross' cannot be negative!")
	}
	else if(cross>n){
		stop("Cross-validation argument 'cross' cannot be larger than the number of samples (",n,").",sep="")
	}
	
	# Return storage preparation
	labels_ret = rep(0L, nbClass)
	nr_W = if(isMulticlass) nbClass else 1
	nc_W = if(bias) p+1 else p
	W_ret=matrix(data=0, nrow=nr_W, ncol=nc_W)
	
	# ACF check
	if (pref_min <= 0) {
        stop ("pref_min must be positive!")
    }
	if (pref_max <= 0) {
        stop ("pref_max must be positive!")
    }
    if (pref_min >= pref_max) {
        stop ("pref_max must be greater than pref_min!")
    }
	if (change_rate <= 0) {
        stop ("Change rate must be positive!")
    }
	if (max_iter <= 0) {
        stop ("Maximum iteration must be positive!")
    }
    
	
	#
	# </Arg preparation>
	# as.double(t(X)) corresponds to rewrite X as a nxp-long vector instead of a n-rows and p-cols matrix. Rows of X are appended one at a time.
	ret <- .C("trainLinear",
			as.double(W_ret),
			as.integer(labels_ret),
			as.double(if(sparse) data@ra else t(data)),
			as.double(yC),
			as.integer(n),
			as.integer(p),
			# sparse index info
			as.integer(sparse),
			as.integer(if(sparse) data@ia else 0),
			as.integer(if(sparse) data@ja else 0),

			as.double(b),
			as.integer(type),
			as.double(cost),
			as.double(epsilon),

			as.integer(nrWi),
			as.double(Wi),
			as.integer(WiLabels),
			as.integer(cross),
	
            as.integer(max_iter),
			
			# acf stuff
			as.double(change_rate),
			as.double(pref_min),
			as.double(pref_max),
			
			as.integer(verbose),
			PACKAGE="LiblineaR.ACF"
			)
	
	if(cross==0){
		labels_ret = ret[[2]]
		
		# classNames is a lookup table for conversion between outputs of C code and initial levels of target
		classNames = {
            if(is.logical(targetLabels)) as.logical(labels_ret)
			else if(is.null(levels(targetLabels))) labels_ret 
			else factor(levels(targetLabels)[labels_ret], levels=levels(targetLabels))
		}
		
		W_ret = matrix(data=ret[[1]], nrow=nr_W, ncol=nc_W)
		colnames(W_ret) = c(
			if(!is.null(colnames(data))) colnames(data) else paste0("W",1:p),
			if(bias) "Bias" else c()
		)
		
		if(isMulticlass)
			rownames(W_ret) = classNames
		
		m=list()
		class(m)="LiblineaR.ACF"
		m$TypeDetail=types[type+1]
		m$Type=type
		m$W=W_ret
		m$Bias=bias
		m$ClassNames=classNames
		m$NbClass=nbClass
		return(m)
	}
	else{
		return(ret[[1]][1])
	}
}
