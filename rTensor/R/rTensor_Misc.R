###Functions that operate on Matrices and Arrays

#'List hadamard Product
#'
#'Returns the hadamard (element-wise) product from a list of matrices or vectors. Commonly used for n-mode products and various Tensor decompositions.
#'@name hadamard_list
#'@rdname hadamard_list
#'@aliases hadamard_list
#'@export
#'@param L list of matrices or vectors
#'@return matrix that is the hadamard product
#'@seealso \code{\link{kronecker_list}}, \code{\link{khatri_rao_list}}
#'@note The modes/dimensions of each element in the list must match.
#'@examples
#'lizt <- list('mat1' = matrix(runif(40),ncol=4), 
#' 'mat2' = matrix(runif(40),ncol=4),
#' 'mat3' = matrix(runif(40),ncol=4))
#'dim(hadamard_list(lizt))
hadamard_list <- function(L){
	isvecORmat <- function(x){is.matrix(x) || is.vector(x)}
	stopifnot(all(unlist(lapply(L,isvecORmat))))
	retmat <- L[[1]]
	for (i in 2:length(L)){
		retmat <- retmat*L[[i]]
	}
	retmat
}

#'List Kronecker Product
#'
#'Returns the Kronecker product from a list of matrices or vectors. Commonly used for n-mode products and various Tensor decompositions.
#'@name kronecker_list
#'@rdname kronecker_list
#'@aliases kronecker_list
#'@export
#'@param L list of matrices or vectors
#'@return matrix that is the Kronecker product
#'@seealso \code{\link{hadamard_list}}, \code{\link{khatri_rao_list}}, \code{\link{kronecker}}
#'@examples
#'smalllizt <- list('mat1' = matrix(runif(12),ncol=4), 
#' 'mat2' = matrix(runif(12),ncol=4),
#' 'mat3' = matrix(runif(12),ncol=4))
#'dim(kronecker_list(smalllizt))
kronecker_list <- function(L){
	isvecORmat <- function(x){is.matrix(x) || is.vector(x)}
	stopifnot(all(unlist(lapply(L,isvecORmat))))
	retmat <- L[[1]]
	for(i in 2:length(L)){
		retmat <- kronecker(retmat,L[[i]])
	}
	retmat
}

#'Khatri-Rao Product
#'
#'Returns the Khatri-Rao (column-wise Kronecker) product of two matrices. If the inputs are vectors then this is the same as the Kronecker product.
#'@name khatri_rao
#'@rdname khatri_rao
#'@aliases khatri_rao
#'@export
#'@param x first matrix
#'@param y second matrix
#'@return matrix that is the Khatri-Rao product
#'@seealso \code{\link{kronecker}}, \code{\link{khatri_rao_list}}
#'@note The number of columns must match in the two inputs.
#'@examples
#'dim(khatri_rao(matrix(runif(12),ncol=4),matrix(runif(12),ncol=4)))
khatri_rao <- function(x,y){
	if (!(is.matrix(x)&&is.matrix(y))) stop("Arguments must be matrices.")
	if (dim(x)[2]!=dim(y)[2]) stop("Arguments must have same number of columns.")
	retmat <- matrix(0,nrow=dim(x)[1]*dim(y)[1],ncol=dim(x)[2])
	for (j in 1:ncol(retmat)) retmat[,j] <- kronecker(x[,j],y[,j])
	retmat
}

#'List Khatri-Rao Product
#'
#'Returns the Khatri-Rao product from a list of matrices or vectors. Commonly used for n-mode products and various Tensor decompositions.
#'@name khatri_rao_list
#'@rdname khatri_rao_list
#'@aliases khatri_rao_list
#'@export
#'@param L list of matrices or vectors
#'@param reverse whether or not to reverse the order
#'@return matrix that is the Khatri-Rao product
#'@seealso \code{\link{khatri_rao}}
#'@note The number of columns must match in every element of the input list.
#'@examples
#'smalllizt <- list('mat1' = matrix(runif(12),ncol=4), 
#' 'mat2' = matrix(runif(12),ncol=4),
#' 'mat3' = matrix(runif(12),ncol=4))
#'dim(khatri_rao_list(smalllizt))
khatri_rao_list <- function(L,reverse=FALSE){
	stopifnot(all(unlist(lapply(L,is.matrix))))
	ncols <- unlist(lapply(L,ncol))
	stopifnot(length(unique(ncols))==1)
	ncols <- ncols[1]
	nrows <- unlist(lapply(L,nrow))
	retmat <- matrix(0,nrow=prod(nrows),ncol=ncols)
	if (reverse) L <- rev(L)
	for(j in 1:ncols){
			Lj <- lapply(L,function(x) x[,j])
			retmat[,j] <- kronecker_list(Lj)
	}
	retmat
}

#'Tensor Times Matrix (m-Mode Product)
#'
#'Contracted (m-Mode) product between a Tensor of arbitrary number of modes and a matrix. The result is folded back into Tensor.
#'@name ttm
#'@rdname ttm
#'@aliases ttm
#'@details By definition, \code{rs_unfold(ttm(tnsr,mat),m) = mat\%*\%rs_unfold(tnsr,m)}, so the number of columns in \code{mat} must match the \code{m}th mode of \code{tnsr}. For the math on the m-Mode Product, see Kolda and Bader (2009).
#'@export
#'@param tnsr Tensor object with K modes
#'@param mat input matrix with same number columns as the \code{m}th mode of \code{tnsr}
#'@param m the mode to contract on
#'@return a Tensor object with K modes
#'@seealso \code{\link{ttl}}, \code{\link{rs_unfold-methods}}
#'@note The \code{m}th mode of \code{tnsr} must match the number of columns in \code{mat}. By default, the returned Tensor does not drop any modes equal to 1.
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009.
#'@examples
#'tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
#'mat <- matrix(runif(50),ncol=5)
#'ttm(tnsr,mat,m=3)
ttm<-function(tnsr,mat,m=NULL){
	stopifnot(is.matrix(mat))
	if(is.null(m)) stop("m must be specified")
	mat_dims <- dim(mat)
	modes_in <- tnsr@modes
	stopifnot(modes_in[m]==mat_dims[2])
	modes_out <- modes_in
	modes_out[m] <- mat_dims[1]
	tnsr_m <- rs_unfold(tnsr,m=m)@data
	retarr_m <- mat%*%tnsr_m
	rs_fold(retarr_m,m=m,modes=modes_out)
}

#'Tensor Times List
#'
#'Contracted (m-Mode) product between a Tensor of arbitrary number of modes and a list of matrices. The result is folded back into Tensor.
#'@name ttl
#'@rdname ttl
#'@aliases ttl
#'@details Performs \code{ttm} repeated for a single Tensor and a list of matrices on multiple modes. For instance, suppose we want to do multiply a Tensor object \code{tnsr} with three matrices \code{mat1}, \code{mat2}, \code{mat3} on modes 1, 2, and 3. We could do \code{ttm(ttm(ttm(tnsr,mat1,1),mat2,2),3)}, or we could do \code{ttl(tnsr,list(mat1,mat2,mat3),c(1,2,3))}. The order of the matrices in the list should obviously match the order of the modes. This is a common operation for various Tensor decompositions such as CP and Tucker. For the math on the m-Mode Product, see Kolda and Bader (2009).
#'@export
#'@param tnsr Tensor object with K modes
#'@param list_mat a list of matrices
#'@param ms a vector of modes to contract on (order should match the order of \code{list_mat})
#'@return Tensor object with K modes
#'@seealso  \code{\link{ttm}}
#'@note The returned Tensor does not drop any modes equal to 1.
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009.
#'@examples
#'tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
#'lizt <- list('mat1' = matrix(runif(30),ncol=3), 
#' 'mat2' = matrix(runif(40),ncol=4),
#' 'mat3' = matrix(runif(50),ncol=5))
#'ttl(tnsr,lizt,ms=c(1,2,3))
ttl<-function(tnsr,list_mat,ms=NULL){
	if(is.null(ms)||!is.vector(ms)) stop ("m modes must be specified as a vector")
	if(length(ms)!=length(list_mat)) stop("m modes length does not match list_mat length")
	num_mats <- length(list_mat)
	if(length(unique(ms))!=num_mats) warning("consider pre-multiplying matrices for the same m for speed")
	mat_nrows <- vector("list", num_mats)
	mat_ncols <- vector("list", num_mats)
	for(i in 1:num_mats){
	mat <- list_mat[[i]]
	m <- ms[i]
	mat_dims <- dim(mat)
	modes_in <- tnsr@modes
	stopifnot(modes_in[m]==mat_dims[2])
	modes_out <- modes_in
	modes_out[m] <- mat_dims[1]
	tnsr_m <- rs_unfold(tnsr,m=m)@data
	retarr_m <- mat%*%tnsr_m
	tnsr <- rs_fold(retarr_m,m=m,modes=modes_out)
	}	
	tnsr	
}

#'Tensor Multiplication (T-MULT)
#'
#'Implements T-MULT based on block circulant matrices (Kilmer et al. 2013) for 3-tensors.
#'
#'@details Uses the Fast Fourier Transform (FFT) speed up suggested by Kilmer et al. 2013 instead of explicitly constructing the block circulant matrix. For the mathematical details of T-MULT, see Kilmer et al. (2013).
#'@export
#'@name t_mult
#'@rdname t_mult
#'@aliases t_mult
#'@param x a 3-tensor
#'@param y another 3-tensor
#'@return tensor product between \code{x} and \code{y}
#'@note This only works (so far) between 3-Tensors.
#'@references M. Kilmer, K. Braman, N. Hao, and R. Hoover, "Third-order tensors as operators on matrices: a theoretical and computational framework with applications in imaging". SIAM Journal on Matrix Analysis and Applications 2013.
#'@examples
#'tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
#'tnsr2 <- new("Tensor",3L,c(4L,3L,5L),data=runif(60))
#'t_mult(tnsr, tnsr2)
t_mult <- function(x,y){
	if((x@num_modes>3)||(y@num_modes>3)) stop("Tensor Multiplication currently only implemented for 3-Tensors")
	modes_x <- x@modes
	modes_y <- y@modes
	if(modes_x[2]!=modes_y[1]) stop("Mode 2 of x and Mode 1 of y must match")
	n3 <- modes_x[3]
	if(n3!=modes_y[3]) stop("Modes 3 of x and y must match")
	#fft's for x and y
	fft_x <- aperm(apply(x@data,MARGIN=1:2,fft),c(2,3,1))
	fft_y <- aperm(apply(y@data,MARGIN=1:2,fft),c(2,3,1))
	#multiply the faces (this is terribad! TO-DO: think of better way!)
	fft_ret <- array(0,dim=c(modes_x[1],modes_y[2],n3))
	for(i in 1:n3){
		first <- fft_x[,,i,drop=FALSE]
		second <- fft_y[,,i,drop=FALSE]
		fft_ret[,,i]<-matrix(first,nrow=dim(first)[1])%*%matrix(second,,nrow=dim(second)[1])
	}
	#ifft and return as Tensor
	#.ifft <- function(x){suppressWarnings(as.numeric(fft(x,inverse=TRUE))/length(x))}
	as.tensor(aperm(apply(fft_ret,MARGIN=1:2, .ifft),c(2,3,1)),drop=FALSE)
}

#####Special Tensors

#'Tensor with Random Entries
#'
#'Generate a Tensor with specified modes with iid normal(0,1) entries.
#'@export
#'@name rand_tensor
#'@rdname rand_tensor
#'@aliases rand_tensor
#'@param modes the modes of the output Tensor
#'@param drop whether or not modes equal to 1 should be dropped
#'@return a Tensor object with modes given by \code{modes}
#'@note Default \code{rand_tensor()} generates a 3-Tensor with modes \code{c(3,4,5)}.
#'@examples
#'rand_tensor()
#'rand_tensor(c(4,4,4))
#'rand_tensor(c(10,2,1),TRUE)
rand_tensor <- function(modes=c(3,4,5),drop=FALSE){
	as.tensor(array(rnorm(prod(modes)), dim=modes),drop=drop)
}

###Matrix Foldings

#'General Folding of Matrix
#'
#'General folding of a matrix into a Tensor. This is designed to be the inverse function to \code{\link{unfold-methods}}, with the same ordering of the indices. This amounts to following: if we were to unfold a Tensor using a set of \code{row_idx} and \code{col_idx}, then we can fold the resulting matrix back into the original Tensor using the same \code{row_idx} and \code{col_idx}.
#'@export
#'@details This function uses \code{aperm} as the primary workhorse.
#'@name fold
#'@rdname fold
#'@aliases fold
#'@param mat matrix to be folded into a Tensor
#'@param row_idx the indices of the modes that are mapped onto the row space
#'@param col_idx the indices of the modes that are mapped onto the column space
#'@param modes the modes of the output Tensor
#'@return Tensor object with modes given by \code{modes}
#'@seealso \code{\link{unfold-methods}}, \code{\link{k_fold}}, \code{\link{unmatvec}}
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009.
#'@examples
#'tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
#'matT3<-unfold(tnsr,row_idx=2,col_idx=c(3,1))
#'identical(fold(matT3,row_idx=2,col_idx=c(3,1),modes=c(3,4,5)),tnsr)
fold <- function(mat, row_idx = NULL, col_idx = NULL, modes=NULL){
	#checks
	rs <- row_idx
	cs <- col_idx
	if(is.null(rs)||is.null(cs)) stop("row space and col space indices must be specified")
	if(is.null(modes)) stop("Tensor modes must be specified")
	if(!is(mat,"Tensor")){
		if(!is.matrix(mat))  stop("mat must be of class 'matrix'")
		}else{
			stopifnot(mat@num_modes==2)
			mat <- mat@data			
			}
	num_modes <- length(modes)
	stopifnot(num_modes==length(rs)+length(cs))
	mat_modes <- dim(mat)
	if((mat_modes[1]!=prod(modes[rs])) || (mat_modes[2]!=prod(modes[cs]))) stop("matrix nrow/ncol does not match Tensor modes")
	#rearranges into array
	iperm <- match(1:num_modes,c(rs,cs))
	as.tensor(aperm(array(mat,dim=c(modes[rs],modes[cs])),iperm))
}

#'k-mode Folding of Matrix
#'
#'k-mode folding of a matrix into a Tensor. This is the inverse funtion to \code{k_unfold} in the m mode. In particular, \code{k_fold(k_unfold(tnsr, m),m,getModes(tnsr))} will result in the original Tensor.
#'@export
#'@details This is a wrapper function to \code{\link{fold}}.
#'@name k_fold
#'@rdname k_fold
#'@aliases k_fold
#'@param mat matrix to be folded into a Tensor
#'@param m the index of the mode that is mapped onto the row indices
#'@param modes the modes of the output Tensor
#'@return Tensor object with modes given by \code{modes}
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009.
#'@seealso \code{\link{k_unfold-methods}}, \code{\link{fold}}, \code{\link{unmatvec}}
#'@examples
#'tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
#'matT2<-k_unfold(tnsr,m=2)
#'identical(k_fold(matT2,m=2,modes=c(3,4,5)),tnsr)
k_fold <- function(mat,m=NULL,modes=NULL){
	if(is.null(m)) stop("mode m must be specified")
	if(is.null(modes)) stop("Tensor modes must be specified")
	num_modes <- length(modes)
	rs <- m
	cs <- (1:num_modes)[-m]
	fold(mat,row_idx=rs,col_idx=cs,modes=modes)
}

#'Unmatvec Folding of Matrix
#'
#'The inverse operation to \code{\link{matvec-methods}}, turning a matrix into a Tensor. For a full account of matrix folding/unfolding operations, consult Kolda and Bader (2009).
#'@export
#'@name unmatvec
#'@rdname unmatvec
#'@aliases unmatvec
#'@param mat matrix to be folded into a Tensor
#'@param modes the modes of the output Tensor
#'@return Tensor object with modes given by \code{modes}
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009.
#'@seealso \code{\link{matvec-methods}}, \code{\link{fold}}, \code{\link{k_fold}}
#'@examples
#'tnsr <- new("Tensor",3L,c(3L,4L,5L),data=runif(60))
#'matT1<-matvec(tnsr)
#'identical(unmatvec(matT1,modes=c(3,4,5)),tnsr)
unmatvec <- function(mat,modes=NULL){
	if(is.null(modes)) stop("Tensor modes must be specified")
	num_modes <- length(modes)
	cs <- 2
	rs <- (1:num_modes)[-2]
	fold(mat,row_idx=rs,col_idx=cs,modes=modes)	
}

#'Row Space Folding of Matrix
#'
#'DEPRECATED. Please see \code{\link{k_fold}}.
#'@export
#'@param mat matrix to be folded
#'@param m the mode corresponding to rs_unfold
#'@param modes the original modes of the tensor
#'@name rs_fold
#'@rdname rs_fold
#'@aliases rs_fold
rs_fold <- function(mat,m=NULL,modes=NULL){
	if(is.null(m)) stop("mode m must be specified")
	if(is.null(modes)) stop("Tensor modes must be specified")
	num_modes <- length(modes)
	rs <- m
	cs <- (1:num_modes)[-m]
	fold(mat,row_idx=rs,col_idx=cs,modes=modes)
}


#'Column Space Folding of Matrix
#'
#'DEPRECATED. Please see \code{\link{unmatvec}}
#'@export
#'@param mat matrix to be folded
#'@param m the mode corresponding to cs_unfold
#'@param modes the original modes of the tensor
#'@name cs_fold
#'@rdname cs_fold
#'@aliases cs_fold
cs_fold <- function(mat,m=NULL,modes=NULL){
	if(is.null(m)) stop("mode m must be specified")
	if(is.null(modes)) stop("Tensor modes must be specified")
	num_modes <- length(modes)
	cs <- m
	rs <- (1:num_modes)[-m]
	fold(mat,row_idx=rs,col_idx=cs,modes=modes)	
}

#'ORL Database of Faces
#'
#'A dataset containing pictures of 40 individuals under 10 different lightings. Each image has 92 x 112 pixels. Structured as a 4-tensor with modes 92 x 112 x 40 x 10.
#'@format A Tensor object with modes 92 x 112 x 40 x 10. The first two modes correspond to the image pixels, the third mode corresponds to the individual, and the last mode correpsonds to the lighting.
#'@source \url{http://www.cl.cam.ac.uk/research/dtg/attarchive/facedatabase.html}
#'@seealso \code{\link{plot_orl}}
"faces_tnsr"

#'Function to plot the ORL Database of Faces
#'
#'A wrapper function to image() to allow easy visualization of faces_tnsr, the ORL Face Dataset.
#'@export
#'@name plot_orl
#'@rdname plot_orl
#'@aliases plot_orl
#'@param subject which subject to plot (1-40)
#'@param condition which lighting condition (1-10)
#'@references AT&T Laboratories Cambridge. \url{http://www.cl.cam.ac.uk/research/dtg/attarchive/facedatabase.html}
#'@references F. Samaria, A. Harter, "Parameterisation of a Stochastic Model for Human Face Identification". IEEE Workshop on Applications of Computer Vision 1994.
#'@seealso \code{\link{faces_tnsr}}
#'@examples
#'plot_orl(subject=5,condition=4)
#'plot_orl(subject=2,condition=7)
plot_orl <- function(subject=1, condition=1){
	if (subject%in%seq(1,40)==FALSE) stop("subject must be between 1 and 40")
	if (condition%in%seq(1,10)==FALSE) stop("condition must be between 1 and 10")
	greyscale = grey(seq(0,1,length=256))
	faces_tnsr <- NULL
	rm(faces_tnsr)
	data(faces_tnsr, package='rTensor', envir=environment())
	image(faces_tnsr[,,subject,condition]@data,col=greyscale)
}


###Invisible Functions (undocumented)
#Wrapper to Inverse FFT
.ifft <- function(x){suppressWarnings(as.numeric(fft(x,inverse=TRUE))/length(x))}
#Creates a superdiagonal tensor
.superdiagonal_tensor <- function(num_modes,len,elements=1L){
	modes <- rep(len,num_modes)
	arr <- array(0, dim = modes)
	if(length(elements)==1) elements <- rep(elements,len)
	for (i in 1:len){
		txt <- paste("arr[",paste(rep("i", num_modes),collapse=","),"] <- ", elements[i],sep="")
		eval(parse(text=txt))
	}
	as.tensor(arr)
}
#3-Tensor Kilmer et. al (2013) identity
.identity_tensor3d <- function(modes){
	if(length(modes)!=3L) stop("identity tensor only implemented for 3d so far")
	n <- modes[1]
	stopifnot(n==modes[2])
	arr <- array(0,dim=modes)
	arr[,,1] <- diag(1,n,n)
	as.tensor(arr)
}
#Simple timing functions
.tic <- function (gcFirst = TRUE,overwrite=TRUE) {
   if(gcFirst) gc(FALSE)
   tic <- proc.time()
   ticExists <- ".tic"%in%ls(all.names=TRUE,envir=baseenv())
   if(overwrite||!ticExists){
   	assign(".tic", tic, envir=baseenv())
   	}
   	else{
   		stop("Another timing function running")
   		}
   invisible(tic)
}
.toc <- function (pr=FALSE) {
   toc <- proc.time()
   tic <- get(".tic", envir=baseenv())
   if(pr) print(toc - tic)
   invisible(toc - tic)
}
# #'MNIST Handwritten Digits Dataset in Tensor Format
# #'
# #'A dataset containing the MNIST training set, which contains 60,000 images (28 x 28 pixels) of handwritten digits (0-9).
# #'@format Organized into a List, with each element of the list being a 3-Tensor corresponding to a single digit. Each 3-Tensor is (N x 28 x 28), where N is the number of samples that correspond to this digit in the original training dataset.
# #'@source \url{http://yann.lecun.com/exdb/mnist/}
# #'@seealso \code{\link{MNIST_test}}, \code{\link{plot_MNIST}}
# "MNIST_train"

# #'MNIST Handwritten Digits Databse in Tensor Format
# #'
# #'A dataset containing the MNIST test set, which contains 10,000 images (28 x 28 pixels) of handwritten digits (0-9).
# #'@format Organized into a List, with each element of the list being a 3-Tensor corresponding to a single digit. Each 3-Tensor is (N x 28 x 28), where N is the number of samples that correspond to this digit in the original test dataset.
# #'@source \url{http://yann.lecun.com/exdb/mnist/}
# #'@seealso \code{\link{MNIST_train}}, \code{\link{plot_MNIST}}
# "MNIST_test"

# #'Function to plot the MNIST Dataset
# #'
# #'A wrapper function to image() to allow easy visualization of MNIST_train and MNIST_test, the training and testing datasets from MNIST.
# #'@export
# #'@name plot_MNIST
# #'@rdname plot_MNIST
# #'@aliases plot_MNIST
# #'@param train_or_test 'train' or 'test' dataset
# #'@param digit which digit (0-9)
# #'@param index which index for this digit (if NULL then a random index will be chosen)
# #'@references \url{http://yann.lecun.com/exdb/mnist/}
# #'@references Y. LeCun, L. Bottou, Y. Bengio, and P. Haffner, "Gradient-based Learning Applied to Document Recognition". Proceedings of the IEEE, 86, 2278-2324.
# #'@seealso \code{\link{MNIST_train}}, \code{\link{MNIST_test}}
# #'@examples
# #'plot_MNIST('train',digit=5,index=14)
# #'plot_MNIST('train',digit=9,index=14)
# #'plot_MNIST('test',digit=3)
# plot_MNIST <- function(train_or_test='train',digit=9,index=NULL){
	# if ((digit%in%seq(0,9))==FALSE) stop("digit must be an integer between 0 and 9")
	# greyscale = grey(seq(0,1,length=256))
	# MNIST_train <- NULL; MNIST_test <- NULL
	# rm(MNIST_train); rm(MNIST_test)
	# if (train_or_test=='train'){
		# data(MNIST_train, package='rTensor', envir=environment())
		# tmp_tnsr <- MNIST_train[[digit+1]]
		# num_length <- tmp_tnsr@modes[1]
		# if (is.null(index)) index <- sample(seq(1,num_length),1)
		# if ((index%in%seq(1,num_length))==FALSE) stop(paste("index must not exceed ",num_length," for this digit",sep=""))
		# image(tmp_tnsr[index,,]@data, col=greyscale) 
	# }else if(train_or_test=='test'){
		# data(MNIST_test, package='rTensor', envir=environment())
		# tmp_tnsr <- MNIST_test[[digit+1]]
		# num_length <- tmp_tnsr@modes[1]
		# if (is.null(index)) index <- sample(seq(1,num_length),1)
		# if ((index%in%seq(1,num_length))==FALSE) stop(paste("index must not exceed ",num_length," for this digit",sep=""))
		# image(tmp_tnsr[index,,]@data, col=greyscale) 
	# }else{
		# stop("train_or_test must be 'train' or 'test'")
	# }
# }
