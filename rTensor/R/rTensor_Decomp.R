###Tensor Decompositions

#'(Truncated-)Higher-order SVD
#'
#'Higher-order SVD of a K-Tensor. Write the K-Tensor as a (m-mode) product of a core Tensor (possibly smaller modes) and K orthogonal factor matrices. Truncations can be specified via \code{ranks} (making them smaller than the original modes of the K-Tensor will result in a truncation). For the mathematical details on HOSVD, consult Lathauwer et. al. (2000).
#'@export
#'@details A progress bar is included to help monitor operations on large tensors.
#'@name hosvd
#'@rdname hosvd
#'@aliases hosvd
#'@param tnsr Tensor with K modes
#'@param ranks a vector of desired modes in the output core tensor, default is \code{tnsr@@modes}
#'@return a list containing the following:\describe{
#'\item{\code{Z}}{core tensor with modes speficied by \code{ranks}}
#'\item{\code{U}}{a list of orthogonal matrices, one for each mode}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)} - if there was no truncation, then this is on the order of mach_eps * fnorm. }
#'}
#'@seealso \code{\link{tucker}}
#'@references L. Lathauwer, B.Moor, J. Vanderwalle "A multilinear singular value decomposition". Journal of Matrix Analysis and Applications 2000.
#'@note The length of \code{ranks} must match \code{tnsr@@num_modes}.
#'@examples
#'tnsr <- rand_tensor(c(6,7,8))
#'hosvdD <-hosvd(tnsr)
#'hosvdD$fnorm_resid
#'hosvdD2 <-hosvd(tnsr,ranks=c(3,3,4))
#'hosvdD2$fnorm_resid
hosvd <- function(tnsr,ranks=NULL){
	stopifnot(is(tnsr,"Tensor"))
	if (sum(ranks<=0)!=0) stop("ranks must be positive")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")
	
	num_modes <- tnsr@num_modes
	#no truncation if ranks not provided
	if(is.null(ranks)){
		ranks <- tnsr@modes
	}else{
		if (sum(ranks>tnsr@modes)!=0) stop("ranks must be smaller than the corresponding mode")
	}
	#progress bar
	pb <- txtProgressBar(min=0,max=num_modes,style=3)
	#loops through and performs SVD on mode-m matricization of tnsr
	U_list <- vector("list",num_modes)
	for(m in 1:num_modes){
		temp_mat <- rs_unfold(tnsr,m=m)@data
		U_list[[m]] <- svd(temp_mat,nu=ranks[m])$u
		setTxtProgressBar(pb,m)
	}
	close(pb)
	#computes the core tensor
	Z <- ttl(tnsr,lapply(U_list,t),ms=1:num_modes)
	est <- ttl(Z,U_list,ms=1:num_modes)
	resid <- fnorm(est-tnsr)
	#put together the return list, and returns
	list(Z=Z,U=U_list,est=est,fnorm_resid=resid)	
}

#'Canonical Polyadic Decomposition
#'
#'Canonical Polyadic (CP) decomposition of a tensor, aka CANDECOMP/PARAFRAC. Approximate a K-Tensor using a sum of \code{num_components} rank-1 K-Tensors. A rank-1 K-Tensor can be written as an outer product of K vectors. There are a total of \code{num_compoents *tnsr@@num_modes} vectors in the output, stored in \code{tnsr@@num_modes} matrices, each with \code{num_components} columns. This is an iterative algorithm, with two possible stopping conditions: either relative error in Frobenius norm has gotten below \code{tol}, or the \code{max_iter} number of iterations has been reached. For more details on CP decomposition, consult Kolda and Bader (2009).
#'@export
#'@details Uses the Alternating Least Squares (ALS) estimation procedure. A progress bar is included to help monitor operations on large tensors.
#'@name cp
#'@rdname cp
#'@aliases cp
#'@param tnsr Tensor with K modes
#'@param num_components the number of rank-1 K-Tensors to use in approximation
#'@param max_iter maximum number of iterations if error stays above \code{tol} 
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following \describe{
#'\item{\code{lambdas}}{a vector of normalizing constants, one for each component}
#'\item{\code{U}}{a list of matrices - one for each mode - each matrix with \code{num_components} columns}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{tucker}}
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009.
#'@examples
#'subject <- faces_tnsr[,,14,]
#'cpD <- cp(subject,num_components=10) 
#'cpD$conv 
#'cpD$norm_percent 
#'plot(cpD$all_resids) 
cp <- function(tnsr, num_components=NULL,max_iter=25, tol=1e-5){
	if(is.null(num_components)) stop("num_components must be specified")
	stopifnot(is(tnsr,"Tensor"))
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	#initialization via truncated hosvd
	num_modes <- tnsr@num_modes
	modes <- tnsr@modes
	U_list <- vector("list",num_modes)
	unfolded_mat <- vector("list",num_modes)
	tnsr_norm <- fnorm(tnsr)
	for(m in 1:num_modes){
		unfolded_mat[[m]] <- rs_unfold(tnsr,m=m)@data
		U_list[[m]] <- matrix(rnorm(modes[m]*num_components), nrow=modes[m], ncol=num_components)
	}
	est <- tnsr
	curr_iter <- 1
	converged <- FALSE
	#set up convergence check
	fnorm_resid <- rep(0, max_iter)
	CHECK_CONV <- function(est){
		curr_resid <- fnorm(est - tnsr)
		fnorm_resid[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
		else{ return(FALSE)}
	}	
	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop (until convergence or max_iter)
	norm_vec <- function(vec){
	norm(as.matrix(vec))
	}
	while((curr_iter < max_iter) && (!converged)){
	setTxtProgressBar(pb,curr_iter)
		for(m in 1:num_modes){
			V <- hadamard_list(lapply(U_list[-m],function(x) {t(x)%*%x}))
			V_inv <- solve(V)			
			tmp <- unfolded_mat[[m]]%*%khatri_rao_list(U_list[-m],reverse=TRUE)%*%V_inv
			lambdas <- apply(tmp,2,norm_vec)
			U_list[[m]] <- sweep(tmp,2,lambdas,"/")	
			Z <- .superdiagonal_tensor(num_modes=num_modes,len=num_components,elements=lambdas)
			est <- ttl(Z,U_list,ms=1:num_modes)
		}
		#checks convergence
		if(CHECK_CONV(est)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)
		}else{
			curr_iter <- curr_iter + 1
			 }
	}
	if(!converged){setTxtProgressBar(pb,max_iter)}
	close(pb)
	#end of main loop
	#put together return list, and returns
	fnorm_resid <- fnorm_resid[fnorm_resid!=0]
	norm_percent<- (1-(tail(fnorm_resid,1)/tnsr_norm))*100
	invisible(list(lambdas=lambdas, U=U_list, conv=converged, est=est, norm_percent=norm_percent, fnorm_resid = tail(fnorm_resid,1),all_resids=fnorm_resid))
}

#'Tucker Decomposition
#'
#'The Tucker decomposition of a tensor. Approximates a K-Tensor using a n-mode product of a core tensor (with modes specified by \code{ranks}) with orthogonal factor matrices. If there is no truncation in one of the modes, then this is the same as the MPCA, \code{\link{mpca}}. If there is no truncation in all the modes (i.e. \code{ranks = tnsr@@modes}), then this is the same as the HOSVD, \code{\link{hosvd}}. This is an iterative algorithm, with two possible stopping conditions: either relative error in Frobenius norm has gotten below \code{tol}, or the \code{max_iter} number of iterations has been reached. For more details on the Tucker decomposition, consult Kolda and Bader (2009).
#'@export
#'@details Uses the Alternating Least Squares (ALS) estimation procedure also known as Higher-Order Orthogonal Iteration (HOOI). Intialized using a (Truncated-)HOSVD. A progress bar is included to help monitor operations on large tensors.
#'@name tucker
#'@rdname tucker
#'@aliases tucker
#'@param tnsr Tensor with K modes
#'@param ranks a vector of the modes of the output core Tensor
#'@param max_iter maximum number of iterations if error stays above \code{tol} 
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following:\describe{
#'\item{\code{Z}}{the core tensor, with modes specified by \code{ranks}}
#'\item{\code{U}}{a list of orthgonal factor matrices - one for each mode, with the number of columns of the matrices given by \code{ranks}}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{hosvd}}, \code{\link{mpca}}
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009.
#'@note The length of \code{ranks} must match \code{tnsr@@num_modes}.
#'@examples
#'tnsr <- rand_tensor(c(4,4,4,4))
#'tuckerD <- tucker(tnsr,ranks=c(2,2,2,2))
#'tuckerD$conv 
#'tuckerD$norm_percent
#'plot(tuckerD$all_resids)
tucker <- function(tnsr,ranks=NULL,max_iter=25,tol=1e-5){
	stopifnot(is(tnsr,"Tensor"))
	if(is.null(ranks)) stop("ranks must be specified")
	if (sum(ranks>tnsr@modes)!=0) stop("ranks must be smaller than the corresponding mode")
	if (sum(ranks<=0)!=0) stop("ranks must be positive")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	#initialization via truncated hosvd
	num_modes <- tnsr@num_modes
	U_list <- vector("list",num_modes)
	for(m in 1:num_modes){
		temp_mat <- rs_unfold(tnsr,m=m)@data
		U_list[[m]] <- svd(temp_mat,nu=ranks[m])$u
	}
	tnsr_norm <- fnorm(tnsr)
	curr_iter <- 1
	converged <- FALSE
	#set up convergence check
	fnorm_resid <- rep(0, max_iter)
	CHECK_CONV <- function(Z,U_list){
		est <- ttl(Z,U_list,ms=1:num_modes)
		curr_resid <- fnorm(tnsr - est)
		fnorm_resid[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
		else{return(FALSE)}
	}
	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop (until convergence or max_iter)
	while((curr_iter < max_iter) && (!converged)){
	setTxtProgressBar(pb,curr_iter)	
	modes <- tnsr@modes
	modes_seq <- 1:num_modes
		for(m in modes_seq){
			#core Z minus mode m
			X <- ttl(tnsr,lapply(U_list[-m],t),ms=modes_seq[-m])
			#truncated SVD of X
			#U_list[[m]] <- (svd(rs_unfold(X,m=m)@data,nu=ranks[m],nv=prod(modes[-m]))$u)[,1:ranks[m]]
			U_list[[m]] <- svd(rs_unfold(X,m=m)@data,nu=ranks[m])$u
		}
		#compute core tensor Z
		Z <- ttm(X,mat=t(U_list[[num_modes]]),m=num_modes)

		#checks convergence
		if(CHECK_CONV(Z, U_list)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)	
		}else{
			curr_iter <- curr_iter + 1
			}
	}
	close(pb)
	#end of main loop
	#put together return list, and returns
	fnorm_resid <- fnorm_resid[fnorm_resid!=0]
	norm_percent<-(1-(tail(fnorm_resid,1)/tnsr_norm))*100
	est <- ttl(Z,U_list,ms=1:num_modes)
	invisible(list(Z=Z, U=U_list, conv=converged, est=est, norm_percent = norm_percent, fnorm_resid=tail(fnorm_resid,1), all_resids=fnorm_resid))
}

#'Multilinear Principal Components Analysis
#'
#'This is basically the Tucker decomposition of a K-Tensor, \code{\link{tucker}}, with one of the modes uncompressed. If K = 3, then this is also known as the Generalized Low Rank Approximation of Matrices (GLRAM). This implementation assumes that the last mode is the measurement mode and hence uncompressed. This is an iterative algorithm, with two possible stopping conditions: either relative error in Frobenius norm has gotten below \code{tol}, or the \code{max_iter} number of iterations has been reached. For more details on the MPCA of tensors, consult Lu et al. (2008).
#'@export
#'@details Uses the Alternating Least Squares (ALS) estimation procedure. A progress bar is included to help monitor operations on large tensors.
#'@name mpca
#'@rdname mpca
#'@aliases mpca
#'@param tnsr Tensor with K modes
#'@param ranks a vector of the compressed modes of the output core Tensor, this has length K-1
#'@param max_iter maximum number of iterations if error stays above \code{tol} 
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following:\describe{
#'\item{\code{Z_ext}}{the extended core tensor, with the first K-1 modes given by \code{ranks}}
#'\item{\code{U}}{a list of K-1 orthgonal factor matrices - one for each compressed mode, with the number of columns of the matrices given by \code{ranks}}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{tucker}}, \code{\link{hosvd}}
#'@references H. Lu, K. Plataniotis, A. Venetsanopoulos, "Mpca: Multilinear principal component analysis of tensor objects". IEEE Trans. Neural networks, 2008.
#'@note The length of \code{ranks} must match \code{tnsr@@num_modes-1}.
#'@examples
#'subject <- faces_tnsr[,,21,]
#'mpcaD <- mpca(subject,ranks=c(10,10))
#'mpcaD$conv
#'mpcaD$norm_percent
#'plot(mpcaD$all_resids)
mpca <- function(tnsr, ranks = NULL, max_iter = 25, tol=1e-5){
	if(is.null(ranks)) stop("ranks must be specified")
	stopifnot(is(tnsr,"Tensor"))
	if (sum(ranks>tnsr@modes)!=0) stop("ranks must be smaller than the corresponding mode")
	if (sum(ranks<=0)!=0) stop("ranks must be positive")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	#initialization via hosvd of M-1 modes
	num_modes <- tnsr@num_modes
	stopifnot(length(ranks)==(num_modes-1))
	ranks <- c(ranks,1)
	modes <- tnsr@modes
	U_list <- vector("list",num_modes)
	unfolded_mat <- vector("list",num_modes)
	for(m in 1:(num_modes-1)){
		unfolded_mat <- rs_unfold(tnsr,m=m)@data
		mode_m_cov <- unfolded_mat%*%t(unfolded_mat)
		U_list[[m]] <- svd(mode_m_cov, nu=ranks[m])$u
	}
	Z_ext <- ttl(tnsr,lapply(U_list[-num_modes],t),ms=1:(num_modes-1))
	tnsr_norm <- fnorm(tnsr)
	curr_iter <- 1
	converged <- FALSE
	#set up convergence check
	fnorm_resid <- rep(0, max_iter)
	CHECK_CONV <- function(Z_ext,U_list){
		est <- ttl(Z_ext,U_list[-num_modes],ms=1:(num_modes-1))
		curr_resid <- fnorm(tnsr - est)
		fnorm_resid[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
		else{return(FALSE)}
	}
	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop (until convergence or max_iter)
	while((curr_iter < max_iter) && (!converged)){
	setTxtProgressBar(pb,curr_iter)
	modes <-tnsr@modes
	modes_seq <- 1:(num_modes-1)
		for(m in modes_seq){
			#extended core Z minus mode m
			X <- ttl(tnsr,lapply(U_list[-c(m,num_modes)],t),ms=modes_seq[-m])
			#truncated SVD of X
			U_list[[m]] <- svd(rs_unfold(X,m=m)@data,nu=ranks[m])$u
		}
		#compute core tensor Z_ext
		Z_ext <- ttm(X,mat=t(U_list[[num_modes-1]]),m=num_modes-1)
		#checks convergence
		if(CHECK_CONV(Z_ext, U_list)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)
		}else{
			curr_iter <- curr_iter + 1
			}
	}
	close(pb)
	#end of main loop
	#put together return list, and returns
	est <- ttl(Z_ext,U_list[-num_modes],ms=1:(num_modes-1))
	fnorm_resid <- fnorm_resid[fnorm_resid!=0]
	norm_percent<-(1-(tail(fnorm_resid,1)/tnsr_norm))*100
	invisible(list(Z_ext=Z_ext, U=U_list, conv=converged, est=est, norm_percent = norm_percent, fnorm_resid=tail(fnorm_resid,1), all_resids=fnorm_resid))
}

#'Population Value Decomposition
#'
#'The default Population Value Decomposition (PVD) of a series of 2D images. Constructs population-level matrices P, V, and D to account for variances within as well as across the images. Structurally similar to Tucker (\code{\link{tucker}}) and GLRAM (\code{\link{mpca}}), but retains crucial differences. Requires \code{2*n3 + 2} parameters to specified the final ranks of P, V, and D, where n3 is the third mode (how many images are in the set). Consult Crainiceanu et al. (2013) for the construction and rationale behind the PVD model.
#'@export
#'@details The PVD is not an iterative method, but instead relies on \code{n3 + 2}separate PCA decompositions. The third mode is for how many images are in the set.
#'@name pvd
#'@rdname pvd
#'@aliases pvd
#'@param tnsr 3-Tensor with the third mode being the measurement mode
#'@param uranks ranks of the U matrices
#'@param wranks ranks of the W matrices
#'@param a rank of \code{P = U\%*\%t(U)}
#'@param b rank of \code{D = W\%*\%t(W)}
#'@return a list containing the following:\describe{
#'\item{\code{P}}{population-level matrix \code{P = U\%*\%t(U)}, where U is constructed by stacking the truncated left eigenvectors of slicewise PCA along the third mode}
#'\item{\code{V}}{a list of image-level core matrices}
#'\item{\code{D}}{population-leve matrix \code{D = W\%*\%t(W)}, where W is constructed by stacking the truncated right eigenvectors of slicewise PCA along the third mode}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'}
#'@references C. Crainiceanu, B. Caffo, S. Luo, V. Zipunnikov, N. Punjabi, "Population value decomposition: a framework for the analysis of image populations". Journal of the American Statistical Association, 2013.
#'@examples
#'subject <- faces_tnsr[,,8,]
#'pvdD<-pvd(subject,uranks=rep(46,10),wranks=rep(56,10),a=46,b=56)
pvd <- function(tnsr,uranks=NULL,wranks=NULL,a=NULL,b=NULL){
	if(tnsr@num_modes!=3) stop("PVD only for 3D")
	if (sum(uranks<=0)!=0) stop("uranks must be positive")
	if (sum(wranks<=0)!=0) stop("wranks must be positive")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	if(is.null(uranks)||is.null(wranks)) stop("U and V ranks must be specified")
	if(is.null(a)||is.null(b)) stop("a and b must be specified")
	modes <- tnsr@modes
	n <- modes[3]
	if(length(uranks)!=n||length(wranks)!=n) stop("ranks must be of length n3")
	pb <- txtProgressBar(min=0,max=(n+3),style=3)
	x <- tnsr@data
	Us <- vector('list',n)
	Vs <- vector('list',n)
	S <- vector('list',n)
	for(i in 1:n){
		svdz <- svd(x[,,i],nu=uranks[i],nv=wranks[i])
		Us[[i]] <- svdz$u
		Vs[[i]] <- svdz$v
		S[[i]] <- svdz$d[1:min(uranks[i],wranks[i])]
		setTxtProgressBar(pb,i)
	}
	U <- matrix(unlist(Us),nrow=modes[1],ncol=sum(uranks)*n)
	#eigenU <- eigen(U%*%t(U))
	P <- eigen(U%*%t(U))$vectors[,1:a] #E-vecs of UU^T
	setTxtProgressBar(pb,n+1)
	V <- matrix(unlist(Vs),nrow=modes[2],ncol=sum(wranks)*n)
	#eigenV <- eigen(V%*%t(V))
	Dt <- eigen(V%*%t(V))$vectors[,1:b] #E-vecs of VV^T
	D <- t(Dt)
	setTxtProgressBar(pb,n+2)
	V2 <- vector('list',n)
	est <- array(0,dim=modes)
	for(i in 1:n){
		V2[[i]] <- (t(P)%*%Us[[i]])%*%diag(S[[i]],nrow=uranks[i],ncol=wranks[i])%*%(t(Vs[[i]])%*%Dt)
		est[,,i] <- P%*%V2[[i]]%*%D
	}
	est <- as.tensor(est)
	fnorm_resid <- fnorm(est-tnsr)	
	setTxtProgressBar(pb,n+3)
	norm_percent<-(1-(fnorm_resid/fnorm(tnsr)))*100
	invisible(list(P=P,D=D,V=V2,est=est,norm_percent=norm_percent,fnorm_resid=fnorm_resid))
}

#'Tensor Singular Value Decomposition
#'
#'TSVD for a 3-Tensor. Constructs 3-Tensors \code{U, S, V} such that \code{tnsr = t_mult(t_mult(U,S),t(V))}. \code{U} and \code{V} are orthgonal 3-Tensors with orthogonality defined in Kilmer et al. (2013), and \code{S} is a 3-Tensor consists of facewise diagonal matrices. For more details on the TSVD, consult Kilmer et al. (2013).
#'@export
#'@name t_svd
#'@rdname t_svd
#'@aliases t_svd
#'@param tnsr 3-Tensor to decompose via TSVD
#'@return a list containing the following:\describe{
#'\item{\code{U}}{the left orthgonal 3-Tensor}
#'\item{\code{V}}{the right orthgonal 3-Tensor}
#'\item{\code{S}}{the middle 3-Tensor consisting of face-wise diagonal matrices}
#'}
#'@seealso \code{\link{t_mult}}, \code{\link{t_svd_reconstruct}}
#'@references M. Kilmer, K. Braman, N. Hao, and R. Hoover, "Third-order tensors as operators on matrices: a theoretical and computational framework with applications in imaging". SIAM Journal on Matrix Analysis and Applications 2013.
#'@note Computation involves complex values, but if the inputs are real, then the outputs are also real. Some loss of precision occurs in the truncation of the imaginary components during the FFT and inverse FFT.
#'@examples
#'tnsr <- rand_tensor()
#'tsvdD <- t_svd(tnsr)
t_svd<-function(tnsr){
	if(tnsr@num_modes!=3) stop("T-SVD only implemented for 3d so far")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	modes <- tnsr@modes
	n1 <- modes[1]
	n2 <- modes[2]
	n3 <- modes[3]
	#progress bar
	pb <- txtProgressBar(min=0,max=n3,style=3)
	#define ifft
	#.ifft <- function(x){suppressWarnings(as.numeric(fft(x,inverse=TRUE))/length(x))}
	#fft for each of the n1n2 vectors (of length n3) along mode 3
	fftz <- aperm(apply(tnsr@data,MARGIN=1:2,fft),c(2,3,1))
	#svd for each face (svdz is a list of the results)
	U_arr <- array(0,dim=c(n1,n1,n3))
	V_arr <- array(0,dim=c(n2,n2,n3))
	m <- min(n1,n2)		
	S_arr <- array(0,dim=c(n1,n2,n3))
	#Think of a way to avoid a loop in the beginning
	#Problem is that svd returns a list but ideally we want 3 arrays
	#Even with unlist this doesn't seem possible
	for (j in 1:n3){
		setTxtProgressBar(pb,j)
		decomp <- svd(fftz[,,j],nu=n1,nv=n2)
		U_arr[,,j] <- decomp$u
		V_arr[,,j] <- decomp$v
		S_arr[,,j] <- diag(decomp$d,nrow=n1,ncol=n2) #length is min(n1,n2)
	}	
	close(pb)
	#for each svd result, we want to apply ifft
	U <- as.tensor(aperm(apply(U_arr,MARGIN=1:2, .ifft),c(2,3,1)))
	V <- as.tensor(aperm(apply(V_arr,MARGIN=1:2, .ifft),c(2,3,1)))
	S <- as.tensor(aperm(apply(S_arr,MARGIN=1:2, .ifft),c(2,3,1)))
	invisible(list(U=U,V=V,S=S))
}

#'Reconstruct Tensor From TSVD
#'
#'Reconstruct the original 3-Tensor after it has been decomposed into \code{U, S, V} via \code{\link{t_svd}}.
#'@export
#'@name t_svd_reconstruct
#'@rdname t_svd_reconstruct
#'@aliases t_svd_reconstruct
#'@param L list that is an output from \code{\link{t_svd}}
#'@return a 3-Tensor 
#'@seealso \code{\link{t_svd}}
#'@examples
#'tnsr <- rand_tensor(c(10,10,10))
#'tsvdD <- t_svd(tnsr)
#'1 - fnorm(t_svd_reconstruct(tsvdD)-tnsr)/fnorm(tnsr)
t_svd_reconstruct <- function(L){
	t_mult(t_mult(L$U,L$S),t(L$V))
}

#####
.is_zero_tensor <- function(tnsr){
	if (sum(tnsr@data==0)==prod(tnsr@modes)) return(TRUE)
	return(FALSE)
}



###t-compress (Not Supported)
.t_compress <- function(tnsr,k){
	modes <- tnsr@modes
	n1 <- modes[1]
	n2 <- modes[2]
	n3 <- modes[3]
	#progress bar
	pb <- txtProgressBar(min=0,max=n3,style=3)
	#define ifft
	#.ifft <- function(x){suppressWarnings(as.numeric(fft(x,inverse=TRUE))/length(x))}
	#fft for each of the n1n2 vectors (of length n3) along mode 3
	fftz <- aperm(apply(tnsr@data,MARGIN=1:2,fft),c(2,3,1))
	#svd for each face (svdz is a list of the results)
	U_arr <- array(0,dim=c(n1,n1,n3))
	V_arr <- array(0,dim=c(n2,n2,n3))
	m <- min(n1,n2)		
	S_arr <- array(0,dim=c(n1,n2,n3))
	#Think of a way to avoid a loop in the beginning
	#Problem is that svd returns a list but ideally we want 3 arrays
	#Even with unlist this doesn't seem possible
	for (j in 1:n3){
		setTxtProgressBar(pb,j)
		decomp <- svd(fftz[,,j],nu=n1,nv=n2)
		U_arr[,,j] <- decomp$u
		V_arr[,,j] <- decomp$v
		S_arr[,,j] <- diag(decomp$d,nrow=n1,ncol=n2) #length is min(n1,n2)
	}	
	close(pb)
	#for each svd result, we want to apply ifft
	U <- as.tensor(aperm(apply(U_arr,MARGIN=1:2, .ifft),c(2,3,1)))
	V <- as.tensor(aperm(apply(V_arr,MARGIN=1:2, .ifft),c(2,3,1)))
	S <- as.tensor(aperm(apply(S_arr,MARGIN=1:2, .ifft),c(2,3,1)))
	
	est <- as.tensor(array(0,dim=modes))
	for (i in 1:k){
		est <- est + t_mult(t_mult(U[,i,,drop=FALSE],S[i,i,,drop=FALSE]),t(V[,i,,drop=FALSE]))
	}
	resid <- fnorm(est-tnsr)
	invisible(list(est=est, fnorm_resid = resid, norm_percent = (1-resid/fnorm(tnsr))*100))
}

###t-compress2 (Not Supported)
.t_compress2 <- function(tnsr,k1,k2){
	A = modeSum(tnsr,m=3,drop=TRUE)
	svdz <- svd(A@data,nu=k1,nv=k2)
	Util <- svdz$u
	Vtil <- svdz$v
	modes <- tnsr@modes
	n3 <- modes[3]
	core <- array(0,dim=c(k1,k2,n3))
	for(i in 1:n3){
	core[,,i]<-t(Util)%*%tnsr[,,i]@data%*%Vtil
	}
	est <- array(0,dim=modes)
	for(i in 1:k1){
		for (j in 1:k2){
			est = est + Util[,i] %o% Vtil[,j] %o% core[i,j,]
		}	
	}
	resid <- fnorm(tnsr - est)
	invisible(list(core = as.tensor(core), est=est, fnorm_resid = resid, norm_percent = (1-resid/fnorm(tnsr))*100))
}
