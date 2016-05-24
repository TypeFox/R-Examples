# A note on dependencies: 
	#parallel is used for mclappy, in As2Vs.
	#ff is used to store large matrices
		# bootstrap PCs are stored as a B-length list
		# of pxK matrices




#      _       _        
#     | |     | |       
#   __| | __ _| |_ __ _ 
#  / _` |/ _` | __/ _` |
# | (_| | (_| | || (_| |
#  \__,_|\__,_|\__\__,_|
                      
#' Functional mean from EEG dataset
#'
#' This package is based on (Fisher et al., 2014), which uses as an example a subset of the electroencephalogram (EEG) measurements from the Sleep Heart Health Study (SHHS) (Quan et al. 1997). Since we cannot publish the EEG recordings from SHHS participants in this package, we instead include the summary statistics of the PCs from our subsample of the processed SHHS EEG data. These summary statistics were generated from measurements of smoothed Normalized Delta Power. This data is used by the \code{\link{simEEG}} to simulate data examples to demonstrate our functions.
#'
#' @details Specifically, \code{EEG_mu} is a vector containing the mean normalized delta power function across all subjects, for the first 7.5 hours of sleep.
#' @seealso \code{\link{EEG_leadingV}}, \code{\link{EEG_score_var}}
#'
#'
#' @references
#' Aaron Fisher, Brian Caffo, and Vadim Zipunnikov. \emph{Fast, Exact Bootstrap Principal Component Analysis for p>1 million}. 2014. http://arxiv.org/abs/1405.0922
#'
#' Stuart F Quan, Barbara V Howard, Conrad Iber, James P Kiley, F Javier Nieto, George T O'Connor, David M Rapoport, Susan Redline, John Robbins, JM Samet, et al.\emph{ The sleep heart health study: design, rationale, and methods}. Sleep, 20(12):1077-1085, 1997. 1.1
#'
#'
#' @name EEG_mu
#' @docType data
#' @keywords data
NULL

#' Leading 5 Principal Components (PCs) from EEG dataset
#'
#' This package is based on (Fisher et al., 2014), which uses as an example a subset of the electroencephalogram (EEG) measurements from the Sleep Heart Health Study (SHHS) (Quan et al. 1997). Since we cannot publish the EEG recordings from SHHS participants in this package, we instead include the summary statistics of the PCs from our subsample of the processed SHHS EEG data. These summary statistics were generated from measurements of smoothed Normalized Delta Power. This data is used by the \code{\link{simEEG}} to simulate data examples to demonstrate our functions.
#' @details Specifically, \code{EEG_leadingV} is a matrix whose columns contain the leading 5 principal components of the EEG dataset.
#' @seealso \code{\link{EEG_mu}}, \code{\link{EEG_score_var}}
#'
#' @references
#' Aaron Fisher, Brian Caffo, and Vadim Zipunnikov. \emph{Fast, Exact Bootstrap Principal Component Analysis for p>1 million}. 2014. http://arxiv.org/abs/1405.0922
#'
#' Stuart F Quan, Barbara V Howard, Conrad Iber, James P Kiley, F Javier Nieto, George T O'Connor, David M Rapoport, Susan Redline, John Robbins, JM Samet, et al.\emph{ The sleep heart health study: design, rationale, and methods}. Sleep, 20(12):1077-1085, 1997. 1.1
#'
#' @name EEG_leadingV
#' @docType data
#' @keywords data
NULL

#' Empirical variance of the first 5 score variables from EEG dataset
#'
#' This package is based on (Fisher et al., 2014), which uses as an example a subset of the electroencephalogram (EEG) measurements from the Sleep Heart Health Study (SHHS) (Quan et al. 1997). Since we cannot publish the EEG recordings from SHHS participants in this package, we instead include the summary statistics of the PCs from our subsample of the processed SHHS EEG data. These summary statistics were generated from measurements of smoothed Normalized Delta Power. This data is used by the \code{\link{simEEG}} to simulate data examples to demonstrate our functions.
#'
#' @details Specifically, \code{EEG_score_var} is a vector containing the variances of the first 5 empirical score variables. Here, we refer to the score variables refer to the \eqn{n}-dimensional, uncorrelated variables, whose coordinate vectors are the principal components \code{\link{EEG_leadingV}}.
#' @seealso \code{\link{EEG_mu}}, \code{\link{EEG_leadingV}}
#'
#'
#' @references
#' Aaron Fisher, Brian Caffo, and Vadim Zipunnikov. \emph{Fast, Exact Bootstrap Principal Component Analysis for p>1 million}. 2014. http://arxiv.org/abs/1405.0922
#'
#' Stuart F Quan, Barbara V Howard, Conrad Iber, James P Kiley, F Javier Nieto, George T O'Connor, David M Rapoport, Susan Redline, John Robbins, JM Samet, et al.\emph{ The sleep heart health study: design, rationale, and methods}. Sleep, 20(12):1077-1085, 1997. 1.1
#'
#' @name EEG_score_var
#' @docType data
#' @keywords data
NULL




#   __                  _   _                 
#  / _|                | | (_)                
# | |_ _   _ _ __   ___| |_ _  ___  _ __  ___ 
# |  _| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# | | | |_| | | | | (__| |_| | (_) | | | \__ \
# |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
                                            
                                            

#' Simulation functional EEG data
#'
#' Our data from (Fisher et al. 2014) consists of EEG measurements from the Sleep Heart Health Study (SHHS) (Quan et al. 1997). Since we cannot publish the EEG recordings from the individuals in the SHHS, we instead include the summary statistics of the PCs from our subsample of the processed SHHS EEG data. This data is used by the \code{simEEG} to simulate functional data that is approximately similar to the data used in our work. The resulting simulated vectors are always of length 900, and are generated from 5 basis vectors (see \code{\link{EEG_leadingV}}).
#'
#' @param n the desired sample size
#' @param centered if TRUE, the sample will be centered to have mean zero for each dimension. If FALSE, measurements will be simulated from a population where the mean is equal to that observed in the sample used in (Fisher et al. 2014) (see \code{\link{EEG_mu}}).
#' @param propVarNoise the approximate proportion of total sample variance attributable to random noise.
#' @param wide if TRUE, the resulting data is outputted as a \code{n} by 900 matrix, with each row corresponding to a different subject. If FALSE, the resulting data is outputted as a 900 by \code{n} matrix, with each column corresponding to a different subject.
#'
#' @return A matrix containing \code{n} simulated measurement vectors of Normalized Delta Power, for the first 7.5 hours of sleep. These vectors are generated according to the equation:
#'
#' \eqn{y = \sum_{j=1}^{5} B_j * s_j + e}
#'
#' Where \eqn{y} is the simulated measurement for a subject, \eqn{B_j} is the \eqn{j^{th}} basis vector, \eqn{s_j} is a random normal variable with mean zero, and e is a vector of random normal noise. The specific values for \eqn{B_j} and \eqn{var(s_j)} are determined from the EEG data sample studied in (Fisher et al., 2014), and are respectively equal to the \eqn{j^{th}} empirical principal component vector (see \code{\link{EEG_leadingV}}), and the empirical variance of the \eqn{j^{th}} score variable (see \code{\link{EEG_score_var}}).
#'
#'
#' @export
#'
#'
#' @references
#' Aaron Fisher, Brian Caffo, and Vadim Zipunnikov. \emph{Fast, Exact Bootstrap Principal Component Analysis for p>1 million}. 2014. http://arxiv.org/abs/1405.0922
#'
#' Stuart F Quan, Barbara V Howard, Conrad Iber, James P Kiley, F Javier Nieto, George T O'Connor, David M Rapoport, Susan Redline, John Robbins, JM Samet, et al.\emph{ The sleep heart health study: design, rationale, and methods}. Sleep, 20(12):1077-1085, 1997. 1.1
#'
#' @examples
#' set.seed(0)
#'
#' #Low noise example, for an illustration of smoother functions
#' Y<-simEEG(n=20,centered=FALSE,propVarNoise=.02,wide=FALSE)
#' matplot(Y,type='l',lty=1)
#'
#' #Higher noise example, for PCA
#' Y<-simEEG(n=100,centered=TRUE,propVarNoise=.5,wide=TRUE)
#' svdY<-fastSVD(Y)
#' V<-svdY$v #since Y is wide, the PCs are the right singular vectors (svd(Y)$v). 
#' d<-svdY$d
#' head(cumsum(d^2)/sum(d^2),5) #first 5 PCs explain about half the variation
#'
#' # Compare fitted PCs to true, generating basis vectors
#' # Since PCs have arbitrary sign, we match the sign of 
#' # the fitted sample PCs to the population PCs first
#' V_sign_adj<- array(NA,dim=dim(V))
#' for(i in 1:5){
#' 	V_sign_adj[,i]<-V[,i] * sign(crossprod(V[,i],EEG_leadingV[,i]))
#' }
#' par(mfrow=c(1,2))
#' matplot(V_sign_adj[,1:5],type='l',lty=1,
#' 		main='PCs from simulated data,\n sign adjusted')
#' matplot(EEG_leadingV,type='l',lty=1,main='Population PCs')
#Up to sign changes, the fitted sample PCs resemble the population PCs
simEEG<-function(n=100, centered=TRUE, propVarNoise=.45,wide=TRUE){
	#start by building Y as tall, then flip if wide==TRUE


	signalVar<-sum(bootSVD::EEG_score_var)
	noiseVar<- signalVar * propVarNoise/(1-propVarNoise)

	noise<- matrix(rnorm(n*900,sd=sqrt(noiseVar/900)),nrow=900,ncol=n) #900 is the dimension of the sample
	scores<-  diag(sqrt(bootSVD::EEG_score_var)) %*% matrix(rnorm(n*5),nrow=5,ncol=n)

	Y<- tcrossprod(bootSVD::EEG_mu,rep(1,n)) + 
		bootSVD::EEG_leadingV%*%scores + noise
	
	if(centered) Y <- t(scale(t(Y), center=TRUE, scale=FALSE)) 

	if(wide) Y<-t(Y)
	return(Y)
}


#' Quickly print an R object's size
#'
#' @param x an object of interest
#' @param units measure to print size in
#'
#' @return print(object.size(x),units=units)
#' @export
#' @examples
#' Y<-simEEG(n=50)
#' os(Y)
os<-function(x,units='Mb') print(object.size(x),units=units)



#' Fast SVD of a wide or tall matrix
#'
#' \code{fastSVD} uses the inherent low dimensionality of a wide, or tall, matrix to quickly calculate its SVD. For a matrix \eqn{A}, this function solves \eqn{svd(A)=UDV'}. 
#' This function can be applied to either standard matrices, or, when the data is too large to be stored in memeory, to matrices with class \code{\link{ff}}. \code{\link{ff}} objects have a representation in memory, but store their contents on disk. In these cases, \code{fastSVD} will implement block matrix algebra to compute the SVD.
#'
#' @param A matrix of dimension (\eqn{n} by \eqn{m}). This can be either of class \code{matrix} or \code{ff}.
#' @param nv number of high dimensional singular vectors to obtain. If \eqn{n>m}, this is the number of \eqn{n}-dimensional left singular vectors to be computed. If \eqn{n<m}, this is the number of \eqn{m}-dimensional right singular vectors to be computed.
#' @param warning_type passed to \code{\link{qrSVD}}, which calculates either \code{svd(tcrossprod(A))} or \code{svd(crossprod(A))}, whichever is of lower dimension.
#' @param center_A  Whether the matrix \code{A} should be centered before taking it's SVD. Centering is done along whichever dimension of \code{A} is larger. For example, if \code{A} is tall, then setting \code{center_A=TRUE} will return the SVD of \code{A} after centering the rows of \code{A}. This centering is implemented as a low dimensional matrix operation that does not require creating a copy of the original matrix \code{A}.
#' @param pattern passed to \code{\link{ff}}. When \code{A} has class \code{ff}, the returned high dimensional singular vectors will also have class \code{ff}. The argument \code{pattern} is passed to \code{\link{ff}} when creating the files on disk for the high dimensional singular vectors.
#'
#' @details Users might also consider changing the global option \code{ffbatchbytes}, from the \code{ff} package. When a \code{ff} object is entered, the \code{ffbatchbytes} option determines the maximum block size in the block matrix algebra used to calculate the SVD.
#'
#' @return Let \eqn{r} be the rank of the matrix \code{A}. \code{fastSVD} solves \eqn{svd(A)=UDV'}, where \eqn{U} is an (\eqn{n} by \eqn{r}) orthonormal matrix, \eqn{D} is an (\eqn{r} by \eqn{r}) diagonal matrix; and \eqn{V} is a (\eqn{m} by \eqn{r}) orthonormal matrix. When \code{A} is entered as an \code{ff} object, the high dimensional singular vectors of \code{A} will be returned as an \code{ff} object as well. For matrices where one dimension is substantially large than the other, calculation times are considerably faster than the standard \code{svd} function.
#'
#' @export
#'
#' @examples
#'
#' 
#' Y<-simEEG(n=100,centered=TRUE,wide=TRUE)
#' svdY<-fastSVD(Y)
#' V<-svdY$v #sample PCs for a wide matrix are the right singular vectors
#' matplot(V[,1:5],type='l',lty=1) #PCs from simulated data
#'
#' #Note: For a tall, demeaned matrix Y, with columns corresponding 
#' #to subjects and rows to measurements, 
#' #the PCs are the high dimensional left singular vectors.
#' 
#' #Example with 'ff'
#' dev.off()
#' library(ff)
#' Yff<-as.ff(Y)
#' Vff<-fastSVD(Yff)$v
#' matplot(Vff[,1:5],type='l',lty=1) 
fastSVD<-function(A,nv=min(dim(A)),warning_type='silent', center_A=FALSE, pattern=NULL){ 
	N<-min(dim(A))
	p<-max(dim(A))
	tall<- dim(A)[1] == p #tall, or big square.
	if(p==N) warning('for square matrix A, this function offers no speed improvement')
	if(center_A & all(dim(A)==p)) warning('Not clear whether to center rows or columns of A. Currently centering rows of A.')
	
	#get AA', or A'A, whichever is (N by N)
	#Let Vn denote N-dim singular vectors, Vp denote p-dim singular vectors, and d denote the singular values.
	if(!tall) AA_NxN<- ffmatrixmult(A,xt=FALSE,yt=TRUE)
	if( tall) AA_NxN<- ffmatrixmult(A,xt=TRUE,yt=FALSE)
	#For example, suppose AA_NxN = AA', then:
	#svd for AA' = UDV'VDU'=U D^2 U'= Vn D^2 Vn'
	#note, we have D=diag(d), where d is a vector and D is a matrix
	if( center_A) M <- diag(N)-matrix(1/N,N,N)
	if(!center_A) M <- diag(N)

	AA_NxN <- M %*% AA_NxN %*% M 
	#above step is equivalent to de_meaning A before taking SVD
	#unless center_A=FALSE, in which case there's no effect.

	svdAA_NxN<-qrSVD(as.matrix(AA_NxN),warning_type=warning_type) 

	#Solve for Vp:
	#wide A: A=Vn D Vp' => (D^-1)(Vn')A=V' => A'(Vn)(D^-1) = Vp
	#tall A: A=Vp D Vn'			 => 		  A (Vn)(D^-1) = Vp 
	#(same results whether A is wide or tall)
	#Or, if we de_mean A with matrix M=de_meaner:
	#wide A:  A'M Vn (D^-1)= Vp
	#tall A:  A M Vn (D^-1)= Vp 
	# arbitrary sign changes in columns of U won't affect anything, they'll just translate into arbitrary sign changes of V columns.
	d<-sqrt(svdAA_NxN$d) 
	Vn<-as.matrix(svdAA_NxN$u[,1:nv])
	if(nv >1) Dinv <-      diag(1/d[1:nv])
	if(nv==1) Dinv <- as.matrix(1/d[1:nv])
	Vp <- ffmatrixmult(A, (M %*% Vn %*% Dinv), xt=(!tall),yt=FALSE, pattern=pattern)
	
	# In svd() convention: A=UDV'
	# (regardless of whether A is wide or tall)
	if(!tall){
		return(list(v=Vp,u=svdAA_NxN$u,d=d)) #U is short, V is long
	}
	if(tall){ 
		return(list(v=svdAA_NxN$u, u=Vp, d=d)) #U is long, V is short
	}
	
}


#' Matrix multiplication with "ff_matrix" or "matrix" inputs
#'
#' A function for \code{crossprod(x,y)}, for \code{tcrossprod(x,y)}, or for regular matrix multiplication, that is compatible with \code{ff} matrices. Multiplication is done without creating new matrices for the transposes of \code{x} or \code{y}. Note, the crossproduct function can't be applied directly to objects with class \code{ff}. 
#'
#' @param x a matrix or ff_matrix
#' @param y a matrix or ff_matrix. If NULL, this is set equal to x, although a second copy of the matrix x is not actually stored.
#' @param xt should the x matrix be transposed before multiplying
#' @param yt should the y matrix be transposed before multiplying (e.g. \code{xt=TRUE}, \code{yt=FALSE} leads to \code{crossprod(x,y)}).
#' @param ram.output force output to be a normal matrix, as opposed to an object with class \code{ff}.
#' @param override.big.error If the dimension of the final output matrix is especially large, \code{ffmatrixmult} will abort, giving an error. This is meant to avoid the accidental creation of very large matrices. Set override.big.error=TRUE to bypass this error.
#' @param ... passed to \code{\link{ff}}.
#' @export
#' @import ff
#' @return A standard matrix, or a matrix with class \code{ff} if one of the input matrices has class \code{ff}.
#' @examples \dontrun{
#'  library(ff)
#' 	#Tall data
#' 	y_tall<-matrix(rnorm(5000),500,10) #y tall
#' 	x_tall<-matrix(rnorm(5000),500,10)
#' 	y_wide<-t(y_tall)
#' 	x_wide<-t(x_tall)
#' 	y_tall_ff<-as.ff(y_tall) #y tall and ff
#' 	x_tall_ff<-as.ff(x_tall) 
#' 	y_wide_ff<-as.ff(y_wide) #y tall and ff
#' 	x_wide_ff<-as.ff(x_wide) 
#'
#'  #Set options to ensure that block matrix algebra is actually done,
#'  #and the entire algebra isn't just one in one step.
#'  #Compare ffmatrixmult against output from standard methods
#'  options('ffbytesize'=100)
#'
#'  #small final matrices
#' 	#x'x
#' 	range(  crossprod(x_tall) - ffmatrixmult(x_tall_ff, xt=TRUE)  )
#' 	range(  tcrossprod(x_wide) - ffmatrixmult(x_wide_ff, yt=TRUE)  )
#' 	range(  crossprod(x_tall,y_tall) - ffmatrixmult(x_tall_ff,y_tall_ff, xt=TRUE)  )
#' 	range(  tcrossprod(x_wide,y_wide) - ffmatrixmult(x_wide_ff,y_wide_ff, yt=TRUE)  )
#' 	range(  (x_wide%*%y_tall) - ffmatrixmult(x_wide_ff,y_tall_ff)  )
#'
#' 	#ff + small data
#' 	s_tall <- matrix(rnorm(80),10,8) 
#' 	s_wide <- matrix(rnorm(80),8,10) 
#'
#' 	#tall output
#' 	range(  crossprod(x_wide, s_tall) - ffmatrixmult(x_wide_ff, s_tall,xt=TRUE)[]  )
#' 	range(  tcrossprod(x_tall, s_wide) - ffmatrixmult(x_tall_ff, s_wide,yt=TRUE)[]  )
#' 	range( x_tall%*%s_tall - ffmatrixmult(x_tall_ff, s_tall)[])
#'
#' 	#Wide output
#' 	range(  crossprod(s_tall, y_wide) - ffmatrixmult( s_tall, y_wide_ff,xt=TRUE)[]  )
#' 	range(  tcrossprod(s_wide, y_tall) - ffmatrixmult( s_wide,y_tall_ff,yt=TRUE)[]  )
#' 	range( s_wide%*%y_wide - ffmatrixmult(s_wide,y_wide_ff)[])
#'
#'  #Reset options for more practical use
#'  options('ffbytesize'=16777216)
#'
#'}
ffmatrixmult <- function(x,y=NULL,xt=FALSE,yt=FALSE,ram.output=FALSE, override.big.error=FALSE,...) {	
	{i1<-NULL; i2<- NULL} #To avoid errors in R CMD check

	dimx<-dim(x)
	if(!is.null(y)) dimy<-dim(y)
	if(is.null(y))  dimy<-dimx

	p <- max(c(dimx,dimy))
	n <- max(min(dimx),min(dimy)) #We assume that both in put matrices have at least one dimension that is managably small. Will allow only 1 dimension of the final output to be greater than n. For two wide or tall matrices as inputs, this ensures that the output is also wide or tall (not large in both dimenions).

	#Find the dimension of the final output (outDim)
	#Check to make sure the "inner dimension" (inDim) that we multiply on matches in these two matrices.
	outDim <-
	inDim  <- rep(NA,2)
	outDim[1] <- dimx[xt+1]
	outDim[2] <- dimy[2-yt]
	inDim[1] <- dimx[2-xt]
	inDim[2] <- dimy[yt+1]
	if(inDim[1]!=inDim[2]) stop('non-conformable arguments')

	if(all(outDim>n) & (!override.big.error)) stop('Returned value is at risk of being extremely large. Both dimensions of output will be fairly large.')

	if(xt & yt) stop('For ff matrix algebra, set only one of xt or yt to TRUE')


	if(all(outDim==n) | (!'ff'%in% c(class(x),class(y)))|ram.output){ 
		out <- matrix(0,outDim[1],outDim[2])
	}else{
		out <- ff(0,dim=outDim,...)
	}

	if(all(outDim==n)){ #Output is square, inDim has the HD part
		if( (xt) &(!yt)) ffapply({
			out<-out+crossprod(x[i1:i2,], y[i1:i2,])
			},X=x,MARGIN=1)
		if((!xt) & (yt)) ffapply({
			out<-out+tcrossprod(x[,i1:i2], y[,i1:i2])
			},X=x,MARGIN=2)
		if((!xt) &(!yt)) ffapply({
			out<-out+x[,i1:i2]%*% y[i1:i2,]
			},X=x,MARGIN=2)
		#Note, if y=NULL, then y[i1:i2,]=NULL, and crossprod will ignore it
	}
	if(outDim[1]>outDim[2] | (outDim[1]==p & outDim[2]==p)){ #output is tall, or big & square
		if( (xt) & (!yt)) ffapply({
			out[i1:i2,]<-crossprod(x[,i1:i2], y)
			},X=x,MARGIN=2)
		if((!xt) &  (yt)) ffapply({
			out[i1:i2,]<-tcrossprod(x[i1:i2,], y)
			},X=x,MARGIN=1)
		if((!xt) & (!yt)) ffapply({
			out[i1:i2,]<-x[i1:i2,]%*% y
			},X=x,MARGIN=1)
		#Here, if y=NULL, we would've already gotten an error
	}
	if(outDim[1]< outDim[2]){ 
		if( (xt) & (!yt))  ffapply({
			out[,i1:i2]<-crossprod(x, y[,i1:i2])
			},X=y,MARGIN=2)
		if((!xt) &  (yt))  ffapply({
			out[,i1:i2]<-tcrossprod(x, y[i1:i2,])
			},X=y,MARGIN=1)
		if((!xt) & (!yt))  ffapply({
			out[,i1:i2]<- x %*% y[,i1:i2]
			},X=y,MARGIN=2)
		#Here, if y=NULL, we would've already gotten an error
	}
	return(out)
}








#' Generate random orthonormal matrix
#'
#' \code{genQ} generates a square matrix of random normal noise, and then takes the QR decomposition to return Q, a random orthogonal square matrix.
#'
#' @param n the dimension of the desired random orthonormal matrix
#' @param lim_attempts the random matrix of normal noise must be full rank to generate the appropriate QR decomposition. \code{lim_attempts} gives the maximum number of attempts for generating a full rank matrix of normal noise.
#'
#' @return a random orthonormal (\eqn{n} by \eqn{n}) matrix
#' @export
#'
#' @examples
#' A<-genQ(3)
#' round(crossprod(A),digits=10)
genQ<-function(n,lim_attempts=200){
	#get full rank noise matrix
	full_rank<-FALSE
	attempt_no<-0
	while((!full_rank) & attempt_no<lim_attempts){
		attempt_no <- attempt_no+1
		normal_mat<-matrix(rnorm(n^2),n,n)
		qr_normal_mat<-qr(normal_mat)
		full_rank<- qr_normal_mat$rank == dim(normal_mat)[1]
	}
	#Output Q matrix
	return(qr.Q(qr_normal_mat))	
}


#' Wrapper for \code{\link{svd}}, which uses random preconditioning to restart when svd fails to converge
#'
#' In order to generate the SVD of the matrix \code{x}, \code{\link{qrSVD}}  calls \code{\link{genQ}} to generate a random orthonormal matrix, and uses this random matrix to precondition \code{x}. The svd of the preconditioned matrix is calculated, and adjusted to account for the preconditioning process in order to find \code{svd(x)}.
#'
#' @param x a matrix to calculate the svd for
#' @param lim_attempts the number of tries to randomly precondition x. We generally find that one preconditioning attempt is sufficient.
#' @param warning_type controls whether the user should be told if an orthogonal preconditioning matrix is required, or if \code{\link{svd}} gives warnings. 'silent' ignores these warnings, 'print' prints the warning to the console, and 'file' saves the warnings in a text file.
#' @param warning_file gives the location of a file to print warnings to, if \code{warning_type} is set to 'file'.
#' @param ... parameters passed to \code{\link{svd}}, such as \code{nv} and \code{nu}.
#'
#' @return Solves \eqn{svd(x)=UDV'}, where \eqn{U} is an matrix containing the left singular vectors of \eqn{x}, \eqn{D} is a diagonal matrix containing the singular values of \eqn{x}; and \eqn{V} is a matrix containing the right singular vectors of \eqn{x} (output follows the same notation convention as the \code{\link{svd}} function).
#' 
#' \code{qrSVD} will attempt the standard \code{\link{svd}} function before preconditioning the matrix \eqn{x}.
#'
#' @seealso \code{\link{fastSVD}}
#'
#' @export
#'
#' @examples
#' x <-matrix(rnorm(3*5),nrow=3,ncol=5)
#' svdx <- qrSVD(x)
#' svdx
qrSVD<-function(x,lim_attempts=50, warning_type='silent',warning_file='qrSVD_warnings.txt', ...){ 

	gotSvd<-FALSE
	dx <- dim(x)
    n <- dx[1]
    p <- dx[2]
	silenceTry<-warning_type %in% c('silent','file')

	#try basic attempt -- if successful, skip rest of function
	try({
		out<-svd(x, ...)
		gotSvd<-TRUE
	},silent=silenceTry)

	#if this fails, try with Q preconditioning matrices
	#Get UDV'= svd(t(Q_n) %*% x %*% Q_p)
	#(Q_nU)D(Q_pV)'= svd(x)
	# nu and nv will translate appropriately (first nv columns of V will yeild first nv columns of Q_p %*% V=Q_p %*% V[,1:nv]).
	attempt_svd <-0
	while( (!gotSvd) & (attempt_svd<lim_attempts) ){
		attempt_svd <- attempt_svd+1
		
		#Generate Q matrices
		Q_n<-genQ(n,lim_attempts= 200)
		if(n==p) Q_p<-Q_n
		if(n!=p) Q_p<-genQ(p,lim_attempts= 200)

		#Try to find svd on the transformed space, then map back
		try({
			UDVt<-svd( crossprod(Q_n,x) %*% Q_p, ...)
			out<-list()
			out$d<- UDVt$d
			out$u<- Q_n %*% UDVt$u
			out$v<- Q_p %*% UDVt$v

			gotSvd<-TRUE
		},silent=silenceTry )

	}

	#give warnings if we had to do extra attempts
	if(attempt_svd>0){
		err_message<- paste('SVD Attained. Number of preconditioning attempts needed =',attempt_svd)

		#print warnings
		if(warning_type=='print') print(err_message)
		
		#Save warnings
		#if file size is already to big, don't print
		print2file<-TRUE
		if(file.exists (warning_file)) # if it exists but it's too big
			if( file.info(warning_file)$size > 50000) print2file<-FALSE
		if(warning_type=='file' & print2file) 
			cat(paste(date(),'-',err_message,'\n'),file=warning_file,append=TRUE)	
	}

	return(out)
}

#' Generate a random set of bootstrap resampling indeces
#'
#' Let \eqn{n} be the original sample size, \eqn{p} be the number of measurements per subject, and \eqn{B} be the number of bootstrap samples. \code{genBootIndeces} generates a (\eqn{B} by \eqn{n}) matrix containing \eqn{B} indexing vectors that can be used to create \eqn{B} bootstrap samples, each of size \eqn{n}.
#' @param B number of desired bootstrap samples
#' @param n size of original sample from which we'll be resampling.
#' @return A (\eqn{B} by \eqn{n}) matrix of bootstrap indeces. Let \code{bInds} denote the output of \code{getBootIndeces}, and \code{Y} denote the original (\eqn{p} by \eqn{n}) sample. Then \code{Y[,bInds[b,]]} is the \eqn{b^{th}} bootstrap sample.
#' @export
#' @examples
#' bInds<-genBootIndeces(B=50,n=200)
genBootIndeces<-function(B,n){
	bInds<-array(NA,dim=c(B,n)) #bootstrap indeces
	for(b in 1:B) bInds[b,]<-sample(n,replace=TRUE)
	return(bInds)
}


#' Calculate bootstrap distribution of \eqn{n}-dimensional PCs
#'
#' \code{bootSVD_LD} Calculates the bootstrap distribution of the principal components (PCs) of a low dimensional matrix. If the score matrix is inputted, the output of \code{bootSVD_LD} can be used to to calculate bootstrap standard errors, confidence regions, or the full bootstrap distribution of the high dimensional components. Most users may want to instead consider using \code{\link{bootSVD}}, which also calculates descriptions of the high dimensional components. Note that \code{\link{bootSVD}} calls \code{bootSVD_LD}.
#'
#' @param UD (optional) a (\eqn{n} by \eqn{n}) matrix of scores, were rows denote individuals, and columns denote measurements in the PC space.
#' @param DUt the transpose of \code{UD}. If both \code{UD} and \code{UDt} are entered and \code{t(UD)!=DUt}, the \code{DUt} argument will override the \code{UD} argument.
#' @param bInds a (\eqn{B} by \eqn{n}) matrix of bootstrap indeces, where \code{B} is the number of bootstrap samples, and \code{n} is the sample size. Each row should be an indexing vector that can be used to generate a new bootstrap sample (i.e. \code{sample(n,replace=TRUE)}). The matrix of bootstrap indeces is taken as input, rather than being calculated within \code{bootSVD_LD}, so that this method can be more easily compared against traditional bootstrap SVD methods on the exact same bootstrap samples. The \code{bInds} matrix can be calculated using the helper function \code{\link{genBootIndeces}}).
#' @param K the number of PCs to be estimated.
#' @param warning_type passed to \code{\link{qrSVD}}, when taking the SVD of the low dimensional bootstrap score matrices.
#' @param verbose if \code{TRUE}, a progress bar will appear.
#' @param centerSamples whether each bootstrap sample should be centered before calculating the SVD.
#'
#' @return For each bootstrap matrix \eqn{(DU')^b}, let \eqn{svd(DU')=:A^b D^b U^b}, where \eqn{A^b} and \eqn{U^b} are (\eqn{n} by \eqn{n}) orthonormal matrices, and \eqn{D^b} is a (\eqn{n} by \eqn{n}) diagonal matrix \eqn{K}. Here we calculate only the first \code{K} columns of \eqn{A^b}, but all \code{n} columns of \eqn{U^b}. The results are stored as a list containing
#'	\item{As}{a \code{B}-length list of the (\code{n} by \code{K}) matrices containing the first \code{K} PCs from each bootstrap sample. This list is indexed by \code{b}, with the \eqn{b^{th}} element containing the results from the \eqn{b^{th}} bootstrap sample.}
#'	\item{ds}{a \code{B}-length list of vectors, indexed by the bootstrap index \code{b}, with each vector containing the singular values of the corresponding bootstrap sample.}
#'	\item{Us}{a \code{B}-length list, indexed by the bootstrap index \code{b}, of the (\eqn{n} by \eqn{n}) matrices \eqn{U^b}.}
#'	\item{time}{The computation time required for the procedure, taken using \code{\link{system.time}}.}
#' If the score matrix is inputted to \code{bootSVD_LD}, the results can be transformed to get the PCs on the original space by multiplying each matrix \eqn{A^b} by the PCs of the original sample, \eqn{V} (see \code{\link{As2Vs}}). The bootstrap scores of the original sample are equal to \eqn{U^b D^b}.
#' @export
#'
#' @examples
#' #use small n, small B, for a quick illustration
#' set.seed(0)
#' Y<-simEEG(n=100, centered=TRUE, wide=TRUE) 
#' svdY<-fastSVD(Y)
#' DUt<- tcrossprod(diag(svdY$d),svdY$u)
#' bInds<-genBootIndeces(B=50,n=dim(DUt)[2])
#' bootSVD_LD_output<-bootSVD_LD(DUt=DUt,bInds=bInds,K=3,verbose=interactive())
bootSVD_LD<-function(UD,DUt=t(UD),bInds=genBootIndeces(B=1000,n=dim(DUt)[2]),K,warning_type='silent',verbose=getOption('verbose'),centerSamples=TRUE){
	B<-dim(bInds)[1]
	n<-dim(DUt)[2]

	#Objects to store results
	dbs<- #diagonals of the matrix D_b
	Ubs<-
	Abs<-list()

	#Do SVDs
	if(verbose) pb<-txtProgressBar(min = 1, max = B,  char = "=", style = 3)
	timeSVD<-system.time({
	for(b in 1:B){		
		#generally we won't include the subscript b here, in our variable names
		DUtP<-DUt[,bInds[b,]] 
		if(centerSamples) DUtP <- t(scale(t(DUtP),center=TRUE,scale=FALSE))
		svdDUtP<-qrSVD(DUtP, warning_type=warning_type)

		#store d_b, switch sign and then store A_b and U_b
		dbs[[b]]<-svdDUtP$d
		signSwitcher <- sign(diag(svdDUtP$u))
		signSwitcher[signSwitcher==0]<-1 #sign() can possibly give us a zero, so we use as.numeric(logic)
		Abs[[b]] <- matrix(NA,nrow=min(dim(DUt)),ncol=K) 
		Ubs[[b]] <- matrix(NA,nrow=n,ncol=K) 
		for(i in 1:K){
			Abs[[b]][,i]<-signSwitcher[i]*svdDUtP$u[,i]
			Ubs[[b]][,i]<-signSwitcher[i]*svdDUtP$v[,i]
		}
		if(verbose) setTxtProgressBar(pb,b)
	}
	})

	return(list(As=Abs,ds=dbs,Us=Ubs,time=timeSVD))
}


#' Convert low dimensional bootstrap components to high dimensional bootstrap components
#'
#' Let \eqn{B} be the number of bootstrap samples, indexed by \eqn{b=1,2,...B}.
#' \code{As2Vs} is a simple function converts the list of principal component (PC) matrices for the bootstrap scores to a list of principal component matrices on the original high dimensional space. Both of these lists, the input and the output of \code{As2Vs}, are indexed by \eqn{b}.
#' @param AsByB a list of the PCs matrices for each bootstrap sample, indexed by \eqn{b}. Each element of this list should be a (\eqn{n} by \eqn{K}) matrix, where \eqn{K} is the number of PCs of interest, and \eqn{n} is the sample size.
#' @param V a tall (\eqn{p} by \eqn{n}) matrix containing the PCs of the original sample, where \eqn{n} is sample size, and \eqn{p} is sample dimension.
#' @param pattern if \code{V} is a class \code{ff} object, the returned value will also be a class \code{ff} object. \code{pattern} is passed to \code{\link{ff}} in creation of the output.
#' @param ... passed to \code{\link{mclapply}}.
#'
#' @return a \code{B}-length list of (\code{p} by \code{K}) PC matrices on the original sample coordinate space (denoted here as \eqn{V^b}). This is achieved by the matrix multiplication \eqn{V^b=VA^b}. Note that here, \eqn{V^b} denotes the \eqn{b^{th}} bootstrap PC matrix, not \eqn{V} raised to the power \eqn{b}. This notation is the same as the notation used in (Fisher et al., 2014).
#' @export
#'
#' @references
#' Aaron Fisher, Brian Caffo, and Vadim Zipunnikov. \emph{Fast, Exact Bootstrap Principal Component Analysis for p>1 million}. 2014. http://arxiv.org/abs/1405.0922
#'
#' @import parallel
#' @examples
#' #use small n, small B, for a quick illustration
#' set.seed(0)
#' Y<-simEEG(n=100, centered=TRUE, wide=TRUE) 
#' svdY<-fastSVD(Y)
#' DUt<- tcrossprod(diag(svdY$d),svdY$u)
#' bInds<-genBootIndeces(B=50,n=dim(DUt)[2])
#' bootSVD_LD_output<-bootSVD_LD(DUt=DUt,bInds=bInds,K=3,verbose=interactive())
#'
#' Vs<-As2Vs(As=bootSVD_LD_output$As,V=svdY$v)
#' # Yields the high dimensional bootstrap PCs (left singular 
#' # vectors of the bootstrap sample Y), 
#' # indexed by b = 1,2...B, where B is the number of bootstrap samples
As2Vs<-function(AsByB, V, pattern=NULL, ...){
	B<-length(AsByB)

	bumpUp<-function(b){
		 	#default ffmatrixmult uses a random tempfile name to store each matrix (for big matrices).
			out <- ffmatrixmult(V,AsByB[[b]],xt=FALSE,yt=FALSE, pattern=paste0(pattern,b,'_'))
			if('ff' %in% class(out)) close(out)
			out
		}

	VsByB<-mclapply(as.list(1:B), FUN=bumpUp, mc.cores=getOption("mc.cores", 1), ...)

	if(length(VsByB)==0) stop('error using mclappy. Parallelization failed.')
	return(VsByB)
}

#' Used for calculation of low dimensional standard errors & percentiles, by re-indexing the \eqn{A^b} by PC index (\eqn{k}) rather than bootstrap index (\eqn{b}).
#'
#' This function is used as a precursor step for calculate bootstrap standard errors, or percentiles. For very high dimensional data, we recommend that the this function be applied to the low dimensional components \eqn{A^b}, but the function can also be used to reorder a list of high dimensional bootstrap PCs. It can equivalently be used to reorder a list of scores. In general, we recommend that as many operations as possible be applied to the low dimensional components, as opposed to their high dimensional counterparts.  This function is called by \code{\link{getMomentsAndMomentCI}}.
#'
#' @param matricesByB a \code{B}-length list of (\code{r} by \code{K}) matrices from each bootstrap sample. If the list elements have class \code{ff}, the returned matrices will also have class \code{ff}.
#' @param pattern (optional) passed to \code{\link{ff}}.
#' @return a \code{K}-length list of (\eqn{B} by \eqn{r}) matrices. If elements of \code{matricesByB} have class \code{ff}, then the returned, reordered matrices will also have class \code{ff}.
#'
#' @export
#' @import ff
#'
#' @examples
#' #use small n, small B, for a quick illustration
#' set.seed(0)
#' Y<-simEEG(n=100, centered=TRUE, wide=TRUE) 
#' svdY<-fastSVD(Y)
#' V<- svdY$v #original sample PCs
#' DUt<- tcrossprod(diag(svdY$d),svdY$u)
#' bInds<-genBootIndeces(B=50,n=dim(DUt)[2])
#' bootSVD_LD_output<-bootSVD_LD(DUt=DUt,bInds=bInds,K=3,verbose=interactive())
#'
#' ########
#' # to get 'low dimensional PC' moments and lower percentiles
#' AsByB<-bootSVD_LD_output$As
#' AsByK<-reindexMatricesByK(AsByB)
#'
#' meanA1<-	apply(AsByK[[1]],2,mean)
#' seA1<-	apply(AsByK[[1]],2,sd)
#' pA1<-	apply(AsByK[[1]],2,function(x) quantile(x,.05))
#' #can also use lapply to get a list (indexed by k=1,...K) of 
#' #the means, standard errors, or percentiles for each PC. 
#' #See example below, for high dimensional bootstrap PCs.
#'
#' #Alternatively, moments can be calculated with
#' seA1_v2<- getMomentsAndMomentCI(As=AsByK,
#' 		V=diag(dim(AsByK[[1]])[2]))$sdPCs[[1]]
#' all(seA1_v2==seA1)
#'
#' #Additional examples of exploring the low dimensional bootstrap 
#' #PC distribution are given in the documentation for 
#' #the 'bootSVD' function.
#' #########
#'
#' #########
#' #High dimensional percentiles for each PC
#' VsByB<-As2Vs(As=AsByB,V=V)
#' VsByK<-reindexMatricesByK(VsByB)
#' percentileCI_Vs<-lapply(VsByK,function(mat_k){
#' 	apply(mat_k,2,function(x) quantile(x,c(.025,.975)))
#' })
#' k=2 # the 2nd PC is a little more interesting here.
#' matplot(t(percentileCI_Vs[[k]]),type='l',lty=1,col='blue')
#' lines(V[,k])
#' ########
#'
#' # Note: This function can also be used to reorganize the
#' #   high dimensional PCs. For 'ff' matrices, this will
#' #   create a new set of files on disk. 
reindexMatricesByK<-function(matricesByB, pattern){
	ff_mat<- ('ff' %in% class(matricesByB[[1]]))

	K<-dim(matricesByB[[1]])[2]
	r<-dim(matricesByB[[1]])[1]
	B<-length(matricesByB)
	matricesByK<-list() 

	#each item in this list has B rows and K cols
	for(k in 1:K){
		if(!ff_mat) matricesByK[[k]] <- matrix(NA,B,r)
		if( ff_mat){
			matricesByK[[k]] <- ff(0,dim=c(B,r), pattern=paste0(pattern,k,'_'))
			close(matricesByK[[k]])
		}
		for(b in 1:B){
			if(ff_mat) {
				open(matricesByB[[b]])
				open(matricesByK[[k]])
			}

			matricesByK[[k]][b,]<-matricesByB[[b]][,k]

			if(ff_mat) {
				close(matricesByB[[b]])
				close(matricesByK[[k]])
			}
		}
	}
	return(matricesByK)
}

#' Used to study of the bootstrap distribution of the k^th singular values, by re-indexing the list of \eqn{d^b} vectors to be organized by PC index (\eqn{k}) rather than bootstrap index (\eqn{b}).
#' @param vectorsByB a \code{B}-length list, containing vectors with the \code{n} values from each bootstrap sample.
#' @return a \code{K}-length list of (\eqn{B} by \eqn{n}) matrices, where each matrices' rows refers to the values from a different bootstrap sample.
#' @export
#'
#' @examples
#' #use small n, small B, for a quick illustration
#' set.seed(0)
#' Y<-simEEG(n=100, centered=TRUE, wide=TRUE) 
#' svdY<-fastSVD(Y)
#' DUt<- tcrossprod(diag(svdY$d),svdY$u)
#' bInds<-genBootIndeces(B=50,n=dim(DUt)[2])
#' bootSVD_LD_output<-bootSVD_LD(DUt=DUt,bInds=bInds,K=3,verbose=interactive())
#' 
#' dsByK<-reindexVectorsByK(bootSVD_LD_output$ds)
#' 
#' boxplot(dsByK[[1]],main='Bootstrap distribution of 1st singular value')
reindexVectorsByK<-function(vectorsByB){
	n<-length(vectorsByB[[1]])
	B<-length(vectorsByB)
	vectorsByK<-list() 
	#each item in this list has b on the rows, k on the cols
	for(k in 1:n){
		dk<- rep(NA,B) #b on the rows
		for(b in 1:B) dk[b]<-vectorsByB[[b]][k]
		vectorsByK[[k]]<-dk
	}
	return(vectorsByK)
}


#' Calculate bootstrap moments and moment-based confidence intervals for the PCs.
#'
#' Let \eqn{K} be the number of PCs of interest, let \eqn{B} be the number of bootstrap samples, and let \eqn{p} be the number of measurements per subject, also known as the dimension of the sample. In general, we use \eqn{k} to refer to the principal component (PC) index, where \eqn{k=1,2,...K}, and use \eqn{b} to refer to the bootstrap index, where \eqn{b=1,2,...B}.
#' @param AsByK a list of the bootstrap PC matrices. This list should be indexed by \eqn{k}, with the \eqn{k^{th}} element of the lsit containing a \eqn{b} by \eqn{p} matrix of results for the \eqn{k^{th}} PC, across bootstrap samples.
#' @param V a (\eqn{p} by \eqn{n}) matrix containing the coordinate vectors for the matrices within the \code{AsByK} list, where \eqn{n} is sample size and \eqn{p} is sample dimension. Generally for bootstrap PCA, \code{AsByK} should contain the PCs for the bootstrap scores, and \code{V} should be the matrix of PCs from the original sample. The argument \code{V} may also be a \code{\link{ff}} object.
#' @param K the number of leading PCs for which moments and confidence intervals should be obtained.
#' @param verbose setting to \code{TRUE} will cause the function to print its progress in calculating the bootstrap variance for each PC.
#' @return a list containing
#'	\item{EVs}{a list containing element-wise bootstrap means for each of the \code{K} fitted PCs, indexed by \code{k}.}
#'	\item{varVs}{a list containing element-wise bootstrap variances for each of the \code{K} fitted PCs, indexed by \code{k}.}
#'	\item{sdVs}{a list containing element-wise bootstrap standard errors for each of the \code{K} fitted PCs, indexed by \code{k}.}
#'	\item{momentCI}{a list of (\eqn{p} by \eqn{2}) matrices, indexed by \code{k}, where \code{momentCI[[k]][j,]} is the pointwise moment-based CI for the \eqn{j^{th}} element of the \eqn{k^{th}} PC.}
#' @export
#' @import ff
#' @examples
#'
#' #use small n, small B, for a quick illustration
#' set.seed(0)
#' Y<-simEEG(n=100, centered=TRUE, wide=TRUE) 
#' svdY<-fastSVD(Y)
#' V<-svdY$v #right singular vectors of the wide matrix Y
#' DUt<- tcrossprod(diag(svdY$d),svdY$u)
#' bInds<-genBootIndeces(B=50,n=dim(DUt)[2])
#' bootSVD_LD_output<-bootSVD_LD(DUt=DUt,bInds=bInds,K=3,verbose=interactive())
#' 
#' AsByB<-bootSVD_LD_output$As
#' AsByK<-reindexMatricesByK(AsByB)
#' moments<-getMomentsAndMomentCI(AsByK,V,verbose=interactive())
#' plot(V[,1],type='l',ylim=c(-.1,.1),main='Original PC1, with CI in blue')
#' matlines(moments$momentCI[[1]],col='blue',lty=1)
#'
#' #Can also use this function to get moments for low dimensional
#' #vectors A^b[,k], by setting V to the identity matrix.
#' moments_A<- getMomentsAndMomentCI(As=AsByK,V=diag(ncol(AsByK[[1]])))
getMomentsAndMomentCI<-function(AsByK,V,K=length(AsByK),verbose=FALSE){
	{i1<-NULL; i2<- NULL} #To avoid errors in R CMD check

	EAs<-lapply(AsByK, colMeans) #EAs is indexed by k
	if(verbose) cat('...calculating expected value for PCs...\n')
	EVs<-lapply(EAs,function(EA) as.matrix(ffmatrixmult(V, matrix(EA,ncol=1),xt=FALSE,yt=FALSE,ram.output=TRUE)) )#V is pxn, EA_k is nx1, so this is overall complexity O(pnK)
	varAs<-lapply(AsByK,var) #indexed by k
	varVs<-list()
	#varTime<-system.time({
	for(k in 1:length(AsByK)){
		# Calculate the diagonals without doing loops using:
		# diag(VA')=rowSums(V*A); or diag(V'A)=colSums(V*A)
		if(verbose) cat(paste0('...calculating variance for PC #',k,'...\n'))

		varVs[[k]]<-ffapply({
			rowSums((V[i1:i2,]%*%varAs[[k]])*V[i1:i2,])
		},X=V, MARGIN=1, RETURN=TRUE, FF_RETURN=FALSE)
	}
	#})#end system.time
	sdVs<-lapply(varVs,sqrt)

	momentCI<-lapply(1:K,function(k){
		return(cbind(EVs[[k]],EVs[[k]])+cbind(-1.96*sdVs[[k]],+1.96*sdVs[[k]]))
	})

	return(list(EPCs=EVs,varPCs=varVs,sdPCs=sdVs,momentCI=momentCI))#varTime=varTime
}


#Note, if VsByB=NULL, it still has 
getHDpercentiles<-function(AsByK,V,K=length(AsByK),percentiles=c(.025,.975),VsByB=NULL,verbose=getOption('verbose')){
	{i1<-NULL; i2<- NULL} #To avoid errors in R CMD check
	
	HDpercentiles<-list()
	B <- dim(AsByK[[1]])[1]
	n <- dim(AsByK[[1]])[2]
	p <- dim(V)[1]

	for(k in 1:K){
		if(verbose) cat(paste0('...calculating percentiles for PC #',k,'...\n'))

		HDpercentiles[[k]]<-matrix(NA,p,length(percentiles))
		if(is.null(VsByB)){ #if we don't have the saved, full HD bootstrap PC distribution (VsByB)
			if(k==1 & verbose) cat('...generating temporary block matrices of HD bootstrap components...\n')
			ffapply({
				Vk_i1_i2 <- tcrossprod(V[i1:i2,], AsByK[[k]])
				HDpercentiles[[k]][i1:i2,] <- t(apply(Vk_i1_i2,1, function(x) quantile(x,percentiles)))
			},DIM=c(p,max(B,n)),VBYTES=.rambytes['double'],MARGIN=1,VMODE='double',FF_RETURN=FALSE) 
		}

		if(!is.null(VsByB)){ #if we have the HD bootstrap PC distribution (VsByB)
		#Note the difference in the DIM argument of ffapply, depending on
		#whether VsByB is null.
			if(k==1 & verbose) cat('...aggregating over calculated HD bootstrap components...\n')
			ff_data <- 'ff' %in% class(VsByB[[1]])
	
			ffapply({
				Vk_i1_i2 <- matrix(0,B,i2-i1+1) #Blank matrix to be filled.
				for(b in 1:B){
					if(ff_data) open(VsByB[[b]]) #open manually, or else we get a warning message
					Vk_i1_i2[b,]<-VsByB[[b]][i1:i2,k]
					if(ff_data) close(VsByB[[b]])
				}
				HDpercentiles[[k]][i1:i2,]<-t(apply(Vk_i1_i2,2,function(x) quantile(x,percentiles)))
			},DIM=c(p,max(B,n)),VBYTES=.rambytes['double'],MARGIN=1,VMODE='double',FF_RETURN=FALSE)
		}
	}

	return(HDpercentiles)
}



#' Calculates bootstrap distribution of PCA (i.e. SVD) results
#'
#' Applies fast bootstrap PCA, using the method from (Fisher et al., 2014). Dimension of the sample is denoted by \eqn{p}, and sample size is denoted by \eqn{n}, with \eqn{p>n}.
#'
#' @param Y initial data sample, which can be either a matrix or a \code{\link{ff}} matrix. \code{Y} can be either tall (\eqn{p} by \eqn{n}) or wide (\eqn{n} by \eqn{p}). If \code{Y} is entered and \code{V}, \code{d} and \code{U} (see definitions below) are not entered, then \code{bootSVD} will also compute the SVD of \code{Y}. In this case where the SVD is computed, \code{bootSVD} will assume that the larger dimension of \code{Y} is \eqn{p}, and the smaller dimension of code {Y} is \eqn{n} (i.e. \code{bootSVD} assumes that (\eqn{p>n}). This assumption can be overriden by manually entering \code{V}, \code{U} and \code{d}.\cr
#' For cases where the entire data matrix can be easily stored in memory (e.g. \eqn{p<50000}), it is generally appropriate to enter \code{Y} as a standard matrix. When \code{Y} is large enough that matrix algebra on \code{Y} is too demanding for memory though, \code{Y} should be entered as a \code{\link{ff}} object, where the actual data is stored on disk. If \code{Y} has class \code{ff}, and \code{V}, \code{d} or \code{U} is not entered, then block matrix algebra will be used to calculate the PCs and bootstrap PCs. The results of these calculations will be returned as \code{\link{ff}} objects as well.
#' @param K number of PCs to calculate the bootstrap distribution for.
#' @param V (optional) the (\eqn{p} by \eqn{n}) full matrix of \eqn{p}-dimensional PCs for the sample data matrix. If \code{Y} is wide, these are the right singular vectors of \code{Y} (i.e. \eqn{Y=UDV'}). If \code{Y} is tall, these are the left singular vectors of \code{Y} (i.e. \eqn{Y=VDU'}). In general it is assumed that \eqn{p>n}, however, this can be overridden by setting \code{V} and \code{U} appropriately.\cr
#' Like \code{Y}, the argument \code{V} can be either a standard matrix or a \code{\link{ff}} matrix. If \code{V} is a \code{ff} object, the bootstrap PCs, if requested, will be returned as \code{\link{ff}} objects as well.
#' @param U (optional) the (\eqn{n} by \eqn{n}) full set of \eqn{n}-dimensional singular vectors of \code{Y}. If \code{Y} is wide, these are the left singular vectors of \code{Y} (i.e. \eqn{Y=UDV'}). If \code{Y} is tall, these are the right singular vectors of \code{Y} (i.e. \eqn{Y=VDU'}).
#' @param d (optional) \eqn{n}-length vector of the singular values of \code{Y}. For example, if \code{Y} is tall, then we have \eqn{Y=VDU'} with \code{D=diag(d)}.
#' @param B number of bootstrap samples to compute.
#' @param output a vector telling which descriptions of the bootstrap distribution should be calculated. Can include any of the following: 'initial_SVD', 'HD_moments', 'full_HD_PC_dist', and 'HD_percentiles'. See below for explanations of these outputs.\cr
#' For especially high dimensional cases, caution should be used if requesting 'full_HD_PC_dist' due to potential storage limitations.
#' @param verbose if TRUE, the function will print progress during calculation procedure.
#' @param bInds a (\eqn{B} by \eqn{n}) matrix of bootstrap indeces, where \code{B} is the number of bootstrap samples, and \code{n} is the sample size. The purpose of setting a specific bootstrap sampling index is to allow the results to be more precisely compared against standard bootstrap PCA calculations. If entered, the \code{bInds} argument will override the \code{B} argument.
#' @param percentiles a vector containing percentiles to be used to calculate element-wise percentiles across the bootstrap distribution (both across the distribution of  \eqn{p}-dimensional components and the distribution of \eqn{n}-dimensional components). For example, \code{percentiles=c(.025,.975)} will return the 2.5 and 97.5 percentiles, which can be used as \eqn{95} percent bootstrap percentile CIs. Alternatively, a longer vector of percentiles can be entered.
#' @param centerSamples whether each bootstrap sample should be centered before calculating the SVD.
#' @param pattern_V if \code{Y} is a class \code{ff} object, then the returned PCs of \code{Y} will also be a class \code{ff} object. \code{pattern_V} is passed to \code{\link{ff}} in creation of the \code{initial_SVD} output. Specifically, \code{pattern_V} is a filename prefix used for storing the high dimensional PCs of the original sample.
#' @param pattern_Vb if \code{Y} or \code{V} is a class \code{ff} object, then the returned bootstrap PCs will also be class \code{ff} objects. \code{pattern_Vb} is passed to \code{\link{ff}} in creation of the \code{full_HD_PC_dist} output. Specifically, \code{pattern_Vb} is a filename prefix used for storing the high dimensional bootstrap PCs.
#'
#' @details Users might also consider changing the global options \code{ffbatchbytes}, from the \code{ff} package, and \code{mc.cores}, from the \code{parallel} package. When \code{ff} objects are entered as arguments for \code{bootSVD}, the required matrix algebra is done using block matrix alebra. The \code{ffbatchbytes} option determines the size of the largest block matrix that will be held in memory at any one time. The \code{mc.cores} option (set to 1 by default) determines the level of parallelization to use when calculating the high dimensional distribution of the bootstrap PCs (see \code{\link{mclapply}}).
#'
#' @return \code{bootSVD} returns a list that can include any of the following elements, depending on what is specified in the \code{output} argument:
#' \describe{
#' 	\item{initial_SVD}{The singular value decomposition of the centered, original data matrix. \code{initial_SVD} is a list containing \code{V}, the matrix of \eqn{p}-dimensional principal components, \code{d}, the vector of singular values of \code{Y}, and \code{U}, the matrix of \eqn{n}-dimensional singular vectors of \code{Y}.}
#'	\item{HD_moments}{A list containing the bootstrap expected value (\code{EPCs}), element-wise bootstrap variance (\code{varPCs}), and element-wise bootstrap standard deviation (\code{sdPCs}) for each of the \eqn{p}-dimensional PCs. Each of these three elements of \code{HD_moments} is also a list, which contains \eqn{K} vectors, one for each PC. \code{HD_moments} also contains \code{momentCI}, a \eqn{K}-length list of (\eqn{p} by 2) matrices containing element-wise moment based confidence intervals for the PCs.}
#'	\item{full_HD_PC_dist}{A \eqn{B}-length list of matrices (or \code{ff} matrices), with the \eqn{b^{th}} list element equal to the (\eqn{p} by \eqn{K}) matrix of high dimensional PCs for the \eqn{b^{th}} bootstrap sample. \cr
#' For especially high dimensional cases when the output is returned as \code{\link{ff}} matrices, caution should be used if requesting 'full_HD_PC_dist' due to potential storage limitations. \cr
#' To reindex these PCs by \eqn{k} (the PC index) as opposed to \eqn{b} (the bootstrap index), see \code{\link{reindexMatricesByK}}. Again though, caution shoulded be used when reindexing PCs stored as \code{ff} objects, as this will double the number of files stored.}
#'	\item{HD_percentiles}{A list of \eqn{K} matrices, each of dimension (\eqn{p} by \eqn{q}), where \eqn{q} is the number of percentiles requested (i.e. \eqn{q} = \code{length(percentiles)}). The \eqn{k^{th}} matrix in \code{HD_percentiles} contains element-wise percentiles for the \eqn{k^{th}}, \eqn{p}-dimensional PC.}
#' }
#'
#'
#'
#' In addition, the following results are always included in the output, regardless of what is specified in the \code{output} argument:
# \describe{
#' \item{full_LD_PC_dist}{A \eqn{B}-length list of matrices, with the \eqn{b^{th}} list element equal to the (\eqn{p} by \eqn{K}) matrix of PCs of the scores in the \eqn{b^{th}} bootstrap sample. To reindex these vectors by \eqn{k} (the PC index), see \code{\link{reindexMatricesByK}}.}
#' \item{d_dist}{A \code{B}-length list of vectors, with the \eqn{b^{th}} element of \code{d_dist} containing the \eqn{n}-length vector of singular values from the \eqn{b^{th}} bootstrap sample. To reindex these values by \eqn{k} (the PC index), see \code{\link{reindexVectorsByK}}.}
#' \item{U_dist}{A \code{B}-length list of (\eqn{n} by \eqn{K}) matrices, with the columns of the \eqn{b^{th}} matrix containing the \eqn{n}-length singular vectors from the \eqn{b^{th}} bootstrap sample. To reindex these vectors by \eqn{k} (the PC index), see \code{\link{reindexMatricesByK}}.}
#' \item{LD_moments}{A list that is comparable to \code{HD_moments}, but that instead describes the variability of the \eqn{n}-dimensional principal components of the resampled score matrices. \code{LD_moments} contains the bootstrap expected value (\code{EPCs}), element-wise bootstrap variances (\code{varPCs}), and element-wise bootstrap standard deviations (\code{sdPCs}) for each of the \eqn{n}-dimensional PCs. Each of these three elements of \code{LD_moments} is also a list, which contains \eqn{K} vectors, one for each PC. \code{LD_moments} also contains \code{momentCI}, a list of \eqn{K} (\eqn{n} by 2) matrices containing element-wise, moment-based confidence intervals for the PCs.}
#' \item{LD_percentiles}{A list of \eqn{K} matrices, each of dimension (\eqn{p} by \eqn{q}), where \eqn{q} is the number of percentiles requested (i.e. \eqn{q} = \code{length(percentiles)}). The \eqn{k^{th}} matrix in \code{LD_percentiles} contains element-wise percentiles for the \eqn{k^{th}} \eqn{n}-dimensional PC.
#' }
#'
#'
#' @references
#' Aaron Fisher, Brian Caffo, and Vadim Zipunnikov. \emph{Fast, Exact Bootstrap Principal Component Analysis for p>1 million}. 2014. http://arxiv.org/abs/1405.0922
#'
#' @export
#' @import ff
#' @examples
#' #use small n, small B, for a quick illustration
#' set.seed(0)
#' Y<-simEEG(n=100, centered=TRUE, wide=TRUE) 
#' b<-bootSVD(Y, B=50, K=2, output= 
#'  	c('initial_SVD', 'HD_moments', 'full_HD_PC_dist',
#'  	'HD_percentiles'), verbose=interactive())
#' 
#' #explore results
#' matplot(b$initial_SVD$V[,1:4],type='l',main='Fitted PCs',lty=1)
#' legend('bottomright',paste0('PC',1:4),col=1:4,lty=1,lwd=2)
#'
#' ######################
#' # look specifically at 2nd PC
#' k<-2
#'
#' ######
#' #looking at HD variability
#'
#' #plot several draws from bootstrap distribution
#' VsByK<-reindexMatricesByK(b$full_HD_PC_dist)
#' matplot(t(VsByK[[k]][1:20,]),type='l',lty=1,
#' 		main=paste0('20 Draws from bootstrap\ndistribution of HD PC ',k))
#'
#' #plot pointwise CIs
#' matplot(b$HD_moments$momentCI[[k]],type='l',col='blue',lty=1,
#' 		main=paste0('CIs For HD PC ',k))
#' matlines(b$HD_percentiles[[k]],type='l',col='darkgreen',lty=1)
#' lines(b$initial_SVD$V[,k])
#' legend('topright',c('Fitted PC','Moment CIs','Percentile CIs'),
#' 		lty=1,col=c('black','blue','darkgreen'))
#' abline(h=0,lty=2,col='darkgrey')
#'
#' ######
#' # looking at LD variability
#'
#' # plot several draws from bootstrap distribution
#' AsByK<-reindexMatricesByK(b$full_LD_PC_dist)
#' matplot(t(AsByK[[k]][1:50,]),type='l',lty=1,
#' 		main=paste0('50 Draws from bootstrap\ndistribution of LD PC ',k),
#'		xlim=c(1,10),xlab='PC index (truncated)')
#'
#' # plot pointwise CIs
#' matplot(b$LD_moments$momentCI[[k]],type='o',col='blue',
#' 		lty=1,main=paste0('CIs For LD PC ',k),xlim=c(1,10),
#' 		xlab='PC index (truncated)',pch=1)
#' matlines(b$LD_percentiles[[k]],type='o',pch=1,col='darkgreen',lty=1)
#' abline(h=0,lty=2,col='darkgrey')
#' legend('topright',c('Moment CIs','Percentile CIs'),lty=1,
#' 		pch=1,col=c('blue','darkgreen'))
#' #Note: variability is mostly due to rotations with the third and fourth PC.
#' 
#' # Bootstrap eigenvalue distribution
#' dsByK<-reindexVectorsByK(b$d_dist)
#' boxplot(dsByK[[k]]^2,main=paste0('Covariance Matrix Eigenvalue ',k),
#' 		ylab='Bootstrap Distribution',
#' 		ylim=range(c(dsByK[[k]]^2,b$initial_SVD$d[k]^2)))
#' points(b$initial_SVD$d[k]^2,pch=18,col='red')
#' legend('bottomright','Sample Value',pch=18,col='red')
#'
#'
#' ##################
#' #Example with ff input
#' library(ff)
#' Yff<-as.ff(Y, pattern='Y_')
#' # If desired, change options in 'ff' package to
#' # adjust the size of matrix blocks held in RAM.
#' # For example:
#' # options('ffbatchbytes'=100000)
#' ff_dir<-tempdir()
#' pattern_V <- paste0(ff_dir,'/V_')
#' pattern_Vb <- paste0(ff_dir,'/Vb_')
#' bff <- bootSVD(Yff, B=50, K=2, output=c('initial_SVD', 'HD_moments',
#'  	'full_HD_PC_dist', 'HD_percentiles'), pattern_V= pattern_V,
#'  	pattern_Vb=pattern_Vb, verbose=interactive())
#' 
#' 
#' # Note that elements of full_HD_PC_dist and initial_SVD
#' # have class 'ff'
#' lapply(bff,function(x) class(x[[1]]))
#' #Show some results of bootstrap draws
#' plot(bff$full_HD_PC_dist[[1]][,k],type='l')
#' #Reindexing by K will create a new set of ff files.
#' VsByKff<-reindexMatricesByK(bff$full_HD_PC_dist,
#'  	pattern=paste0(ff_dir,'/Vk_'))
#' physical(bff$full_HD_PC_dist[[1]])$filename
#' physical(VsByKff[[1]])$filename
#' matplot(t(VsByKff[[k]][1:10,]),type='l',lty=1,
#' main=paste0('Bootstrap Distribution of PC',k))
#' 
#' 
#' # Saving and moving results:
#' saveRDS(bff,file=paste0(ff_dir,'/bff.rds'))
#' close(bff$initial_SVD$V)
#' physical(bff$initial_SVD$V)$filename
#' # If the 'ff' files on disk are moved or renamed,
#' # this filename attribute can be changed:
#' old_ff_path <- physical(bff$initial_SVD$V)$filename
#' new_ff_path <- paste0(tempdir(),'/new_V_file.ff')
#' file.rename(from= old_ff_path, to= new_ff_path)
#' physical(bff$initial_SVD$V)$filename <- new_ff_path
#' matplot(bff$initial_SVD$V[,1:4],type='l',lty=1)
#' 
bootSVD<-function(Y=NULL,K,V=NULL,d=NULL,U=NULL,B=50,output='HD_moments',verbose=getOption('verbose'),bInds=NULL,percentiles=c(.025,.975),centerSamples=TRUE, pattern_V='V_', pattern_Vb='Vb_'){

	#Input checks and warnings
	if( (object.size(Y) > getOption('ffbatchbytes',16777216)) & (!'ff' %in% c(class(Y),class(V))) )
		stop(paste0("Very large sample data matrices, here ",object.size(Y)," bytes, can lead to very large memory requirements. For especially high dimensional data, either input Y and/or V as 'ff' objects, or increase the 'ffbatchbytes' option. This option is currently set to ",getOption("ffbatchbytes", 16777216)," bytes, which is less than the size of Y (or V). Inputting an object with class 'ff' will implement a block matrix algebra procedure."))
	svd_args<- FALSE
	if(!is.null(V) & !is.null(U) & !is.null(d)) svd_args<- TRUE
	if(svd_args){
		if(dim(V)[2] != dim(U)[2] | dim(U)[2]!=length(d)) stop("'V', 'U', and 'd' have incompatible dimensions. Please include all PC basis vectors in 'V'.")
	}

	# Set up timer
	timer<-list()
	

	#get initial SVD, if needed
	timer$svd<-NA
	if(any(is.null(V),is.null(d),is.null(U))){
		timer$svd<-system.time({
			if(verbose) cat('Getting initial svd(Y)...\n')
			tall <- dim(Y)[1]>dim(Y)[2]

			#Get SVD
			svdY<-fastSVD(Y, center_A=TRUE, pattern=pattern_V)		

			if(!tall){
				V<-svdY$v
				U<-svdY$u
			}
			if( tall){
				U<-svdY$v
				V<-svdY$u
			}
			d<-svdY$d
			rm('svdY')
		})
	}
	n<-dim(U)[1] 
	p<-dim(V)[1]
	DUt<- tcrossprod(diag(d),U)

	timer$LD_svds<-system.time({
		#if desired, use inputted bInds matrix
		if(!is.null(bInds)) B<-dim(bInds)[1]
		if(is.null(bInds)) bInds<-genBootIndeces(B=B,n=dim(DUt)[2])

		if(verbose) cat('...calculating n-dimensional bootstrap SVDs...\n')

		bootSVD_LD_output<-bootSVD_LD(DUt=DUt,bInds=bInds,K=K,verbose=verbose,centerSamples=centerSamples)
		AsByB<-bootSVD_LD_output$As
		AsByK<-reindexMatricesByK(matricesByB=AsByB)
	})

	# Build output
	out_contents<-list()

	############
	#always include
	timer$LD_moments<-system.time({
		out_contents[['LD_moments']]<-getMomentsAndMomentCI(AsByK,diag(min(n,p)),verbose=FALSE) #it's too false to bother reporting in `verbose`.
	})

	timer$LD_percentiles<-system.time({
		out_contents[['LD_percentiles']]<-lapply(AsByK,function(mat_k){
			t(apply(mat_k,2,function(x) quantile(x,percentiles)))
		})
	})

	out_contents[['full_LD_PC_dist']]<-AsByB
	
	out_contents[['d_dist']]<-bootSVD_LD_output$ds
	
	out_contents[['U_dist']]<-bootSVD_LD_output$Us
	

	#Test whether we're using 'ff' objects
	ff_data <- 'ff' %in% class(V)

	##############

	##############
	# Sometimes include
	if('initial_SVD' %in% output){
		out_contents[['initial_SVD']]<-list(V=V,d=d,U=U)
	}
	
	timer$HD_moments <- NA
	if('HD_moments' %in% output){
		if(.rambytes['double']*p > getOption('ffbatchbytes',16777216))
			warning(paste0('HD moments will contain vectors with length equal to ',p,' (sample dimension) and size equal to ',.rambytes['double']*p," bytes. These vectors will be stored as RAM objects despite their size exceeding getOption('ffbatchbytes',16777216) = ",getOption('ffbatchbytes',16777216),' bytes.'))

		timer$HD_moments <- system.time({
			out_contents[['HD_moments']]<-getMomentsAndMomentCI(AsByK,V,verbose=verbose)
		})
	}
	

	#### Full HD Distribution 
	timer$full_HD_PC_dist <- NA
	if('full_HD_PC_dist' %in% output){

		if(verbose) cat('...calculating HD Bootstrap PC distribution...\n')
		timer$full_HD_PC_dist<-system.time({
			out_contents[['full_HD_PC_dist']]<- As2Vs(AsByB,V=V, pattern=pattern_Vb)#list of B, p by K matrices
		})
	}

	#Make sure you've gotten the full HD PC dist (if required) before the percentiles.
	timer$HD_percentiles <- NA
	if('HD_percentiles' %in% output){

		#Warning
		if(.rambytes['double']*p*length(percentiles) > getOption('ffbatchbytes',16777216)){
			warning(paste0('HD percentiles will contain matrices with dimension equal to (',p,' by ',length(percentiles),') and size equal to ',.rambytes['double']*p*length(percentiles)," bytes. These matrices will be stored as RAM objects despite their size exceeding getOption('ffbatchbytes',16777216) = ",getOption('ffbatchbytes',16777216),' bytes.'))
		}

		timer$HD_percentiles<-system.time({


			#If the HD PC distribution hasn't been saved
			if(is.null(out_contents[['full_HD_PC_dist']]))
				out_contents[['HD_percentiles']] <- getHDpercentiles(AsByK=AsByK, V=V, K=K, VsByB=NULL, percentiles=percentiles, verbose=verbose)

			#If the HD PC distribution has been saved, 
			#scanning over that stored data to get 
			#percentiles was shown in some limited tests 
			#to be faster.
			if(!is.null(out_contents[['full_HD_PC_dist']]))
				out_contents[['HD_percentiles']] <- getHDpercentiles(AsByK=AsByK, V=V, K=K, VsByB=out_contents[['full_HD_PC_dist']], percentiles=percentiles, verbose=verbose)
		})
	}


	#Clean up:
	if(!'initial_SVD' %in% output){
		if(ff_data) attr(attr(V, 'physical'),'finalizer') <- "delete"
		rm(V)
	}
	if(ff_data) gc()

	out_contents[['timer']] <- timer

	return(out_contents)
}

#' Quickly calculates bootstrap PCA results (wrapper for bootSVD)
#'
#' All arguments are passed to \code{\link{bootSVD}}. This function should be used in exactly the same way as \code{\link{bootSVD}}. The only difference is that PCA typically involves re-centering each bootstrap sample, whereas calculations involving the SVD might not.
#'
#' @param centerSamples whether each bootstrap sample should be centered before computing the bootstrap principal components.
#' @param ... passed to \code{\link{bootSVD}}
#' @return \code{bootSVD(...)}
#'
#' @export
#' 
bootPCA<-function( centerSamples=TRUE, ... ) bootSVD(...)






