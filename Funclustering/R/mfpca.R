#'
#' This function will run a weighted functional pca in the two cases of uni, and multivariate cases. 
#' If the observations (the curves) are given with weights, set up the parameter tik. 
#' @title Multivariate functional pca
#' 
#' 
#' @param fd in the univariate case fd is an object from a class fd. 
#' Otherwise in the multivariate case fd is a list of fd object (fd=list(fd1,fd2,..)). 
#' 
#' @param nharm number of harmonics or principal component to be retain.
#' 
#' @param tik the weights of the functional pca which corresponds to the weights of the curves. 
#' If don't given, then we will run a classic functional pca (without weighting the curves). 
#' 
#' @return When univarite functional data, the function are returning an object of calss "pca.fd", 
#' When multivariate a list of "pca.fd" object by dimension. The "pca.fd" class contains the folowing parameter:
#' 			harmonics: functional data object storing the eigen function
#' 			values: the eigenvalues
#' 			varprop: the normalized eigenvalues (eigenvalues divide by their sum)
#' 			scores: the scores matrix  
#' 			meanfd: the mean of the functional data object
#' 			
#' @export 
#' @examples 
#' 
#' data(growth)
#' data=cbind(matrix(growth$hgtm,31,39),matrix(growth$hgtf,31,54));
#' t=growth$age;
#' splines <- create.bspline.basis(rangeval=c(1, max(t)), nbasis = 20,norder=4);
#' fd <- Data2fd(data, argvals=t, basisobj=splines);
#' pca=mfpca(fd,nharm=2)
#' summary(pca)
#' 
#' @useDynLib Funclustering

mfpca <- function(fd,nharm,tik=numeric(0)) {
	#cheking for functional data object
	if(missing(fd)){
		stop("fd is missing.")
	}
	
	#cheking if fd is an univariate or multivariate fd obj
	if(class(fd)=="list"){
		for(i in 1:length(fd)){
			if(class(fd[[i]])!="fd"){
				stop("the dimension ", i , " of the multivarite object is not a functional data object")
			}
		}
		# getting the coefs and basisProd
		fdData=cppMultiData(fd)
		coefs=fdData$coefs
		basisProd=fdData$basisProd
	}
	else{
		if(class(fd)=="fd"){
			# getting the coefs and basisProd
			fdData=cppUniData(fd)
			coefs=fdData$coefs
			basisProd=fdData$basisProd
		}
		else{
			stop("fd is not a functional data object")
		}
	}
	# getting the number of curves
	n=nrow(coefs)
	
	#cheking for tik 
	if(length(tik)==0){
		tik=rep(1/n, n)
	}
	
	# here the input object, to the C++ code. NB: nbClust and thd=0.05 has no effect, it's just 
	  # to further use to a C++ constructor  of IModel 
	inpobj=new("Input",coefs=coefs, basisProd=basisProd, K=2,thd=0.05)
	
	### call the c++ function which compute functional pca in c++: 
	nharm=as.integer(nharm)
	
	pca=.Call("mfpcaCpp",inpobj,nharm,tik,PACKAGE="Funclustering")
	
	# if fd is a univariate functional data, then pca is exactly the result
	if(class(fd)=="fd"){
		harmonics=fd(coef=pca$harmonics,basisobj=fd$basis,fdnames=fd$fdnames)
		mfpcaReturn=list(harmonics=harmonics,values=pca$values,scores=pca$scores,
				varprop=pca$varprop,meanfd=mean.fd(fd))
		class(mfpcaReturn)="pca.fd"
	}
	else{
		## at the moment, we have the harmonics of all dimensions in the same matrix, we have to cut 
		#  this matrix and make the harmonics matrix for each dimension in a specific matrix
		
		# dim of the functional data
		dimFd=length(fd)
		# nbasis for each dimension
		nbasis=c()
		for(i in 1:dimFd){
			nbasis[i]=fd[[i]]$basis$nbasis
		}
		
		# cut the harmonicis coefficients 
		harmcoefs=harmsCut(pca$harmonics,nbasis)
		
		harmonics=list()#list of functional harmonics: the eigen functions
		mfpcaReturn=list()
		# construction de l'objet de class pca pour chaque dim
		for(i in 1:dimFd){
			# Build of the eigen function for the dimension i
			harmonics[[i]]=fd(coef=harmcoefs[[i]],basisobj=fd[[i]]$basis,fdnames=fd[[i]]$fdnames)
			
			# build of the functional pca object for the dimension i
		    mfpcaReturn[[i]]=list(harmonics=harmonics[[i]],values=pca$values,scores=pca$scores,
					varprop=pca$varprop,meanfd=mean.fd(fd[[i]]))
			class(mfpcaReturn[[i]])="pca.fd"
		}
	}
	return(mfpcaReturn)
}
