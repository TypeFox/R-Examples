#'Kernel function
#'
#'@param x scalar
#'
#'@return density
#'
#'@export
ker <- function(x){
  1/sqrt(2*pi)*exp(-x^2/2);
}
#'Calculating estimated density value on some x with bandwidth h
#'
#'@param x scalar
#'@param h bandwidth
#'@param data data sample
#'
#'@return estimated density value
#'
#'@examples
#'fhat(0,0.2,rnorm(100))
#'
#'@export
fhat <- function(x,h,data){
  n <- length(data);
  return(1/(n*h)*sum(ker((x-data)/h)));
}
#'Get number order matrix which is used in md
#'
#'@param i integer
#'@param length the number of grids in domain
#'
#'@return matrix related to number order
#'
#'@export
nom <- function(i,length){
  m <- matrix(0,nrow=length,ncol=2);
  m[,1] <- i;
  m[,2] <- 1:length;
  return(m)
}
#'Get number order vector which is used in md
#'
#'@param j integer
#'
#'@return vector related to number order
#'
#'@export
xmm <- function(j){
  i <- j-1;
  m <- matrix(0,nrow=i,ncol=2);
  m[,1] <- 1:i;
  m[,2] <- j;
  xmmsub <- matrix(t(m),nrow=1,ncol=(i*2));
  return(xmmsub)
}
#'Get number order matrix which is used in md
#'
#'@param i integer
#'@param prod2 integer which is defined in md
#'
#'@return matrix related to number order
#'
#'@export
dmm <- function(i,prod2){
  m <- matrix(0,nrow=prod2,ncol=2);
  m[,1] <- i;
  m[,2] <- 1:prod2;
  return(m)
}
#'Get estimated values of kernel density estimator on domain
#'
#'@param ij number order vector
#'@param data_for_d data sample which is split to be used for
#'       kernel density estimator
#'@param h bandwidth
#'@param x scalar
#'
#'@return estimated values of kernel density estimator on domain
#'
#'@export
fhatboxm <- function(ij,data_for_d,h,x){
  i <- ij[1];
  j <- ij[2];
  return(fhat(x[j],h[i],data_for_d))
}
#'Calculating Scheffe sets
#'
#'@param ij number order vector
#'@param box estimated values of all kernel density estimators
#'
#'@return 0-1 vector
#'
#'@export
scheffe <- function(ij,box){
  s1 <- ij[1];
  s2 <- ij[2];
  a <- box[s1,]-box[s2,];
  asub <- abs(a);
  b <- ceiling(a/(2*max(asub)+5));
  a2 <- box[s2,]-box[s1,];
  b2 <- ceiling(a2/(2*max(asub)+5));
  return(c(b,b2));
}
#'Auxiliary function which is used in md
#'
#'@param l integer
#'@param data data sample
#'@param x scalar
#'
#'@return integer
#'
#'@export
subcounter <- function(l,data,x){
  e1 <- ceiling((x - data[l])/100);
  e2 <- which.max(e1);
  return(e2)
}
#'Get 0-1 vector which is used for calculating empirical measure
#'
#'@param e2 integer
#'@param box2 matrix which has 0-1 elements related to Scheffe set
#'
#'@return 0-1 vector
#'
#'@export
counter <- function(e2,box2){
  f1 <- box2[,(e2-1)]+box2[,e2];
  f2 <- floor(f1/2);
  return(f2)
}
#'Calculating delta
#'
#'@param ij number order vector
#'@param box matrix which has estimated values of all
#'       kernel density estimators
#'@param box2 matrix which has 0-1 elements related to Scheffe set
#'@param mu_box2 matrix which has values of all empirical measures
#'@param grid length of grid in domain
#'
#'@return delta value
#'
#'@export
deltaboxm <- function(ij,box,box2,mu_box2,grid){
  i <- ij[1];
  j <- ij[2];
  k <- abs(sum(box[i,]*box2[j,]*grid)-mu_box2[j]);
  return(k)
}
#'Get number order matrix which is used in md
#'
#'@param i integer
#'@param length the number of grids in domain
#'
#'@return matrix related to number order
#'
#'@export
xym <- function(i,length){
  m <- matrix(0,nrow=length,ncol=length);
  m[,1] <- i;
  m[,2] <- 1:length;
  return(m)
}
#'md selects bandwidth for kernel density estimator
#'with minimum distance method. Minimum distance method
#'directly selects optimal kernel density estimator
#'in countably infinite kernel density estimators
#'and indirectly selects optimal bandwidth.
#'md selects optimal bandwidth in countably finite
#'kernel density estimators.
#'
#'@param data data sample
#'@param hnumber the number of bandwidth which md can select.
#'               60 is enough. Of course, you can take it more.
#'@param ds rate of data split. Minimum distance method has to split
#'          data for constructing kernel density estimators and
#'          empirical measures.
#'
#'@return bandwidth
#'
#'@usage md(data,hnumber,ds)
#'
#'@examples
#'# select bandwidth
#'md(runif(100),20,0.6)
#'
#'# select bandwidth and plot
#'data <- rnorm(100)
#'bandwidth <- md(data,20,0.6)
#'x <- seq(min(data),max(data),length=100)
#'plot(x,sapply(x,fhat,bandwidth,data),type="l",ylab="density")
#'
#'@export
md <- function(data,hnumber,ds){
	xmin <- min(data);
	xmax <- max(data);
	xlength <- 25*(ceiling(abs(xmax-xmin))+2);
	x <- seq(xmin-1,xmax+1,length=xlength);
	xgrid <- abs(x[2]-x[1]);
	n <- length(data);
	number_kernel_data <- round(n*ds);
	prod <- xlength*hnumber;
	prod2 <- hnumber*(hnumber-1);
	kernel_data <- data[1:number_kernel_data];
	empirical_measure_data <- data[(number_kernel_data+1):n];
	hrange <- 0.79*IQR(empirical_measure_data)*(n-number_kernel_data)^(-1/5);
	h <- seq(0.001,2*hrange,length=hnumber);
	number_h <- 1:hnumber;
	sub1_number_matrix <- sapply(number_h,nom,xlength);
	sub2_number_matrix <- sub1_number_matrix[1:xlength,1:hnumber];
	dim(sub2_number_matrix) <- c(prod,1);
	sub3_number_matrix <- sub1_number_matrix[(xlength+1):(2*xlength),1:hnumber];
	dim(sub3_number_matrix) <- c(prod,1);
	number_matrix <- c(sub2_number_matrix,sub3_number_matrix);
	dim(number_matrix) <- c(prod,2);
	fhat_box <- apply(number_matrix,1,fhatboxm,kernel_data,h,x);
	fhat_box <- matrix(fhat_box,nrow=hnumber,ncol=xlength,byrow=T);
	number_sub1 <- 2:hnumber;
	sub4_number_matrix <- matrix(unlist(sapply(number_sub1,xmm)),nrow=(prod2/2),ncol=2,byrow=T);
	scheffe_matrix <- matrix(apply(sub4_number_matrix,1,scheffe,fhat_box),nrow=prod2,ncol=xlength,byrow=T);
	number_sub2 <- 1:(n-number_kernel_data);
	sub_mubox <- sapply(number_sub2,subcounter,empirical_measure_data,x);
	mubox <- sapply(sub_mubox,counter,scheffe_matrix);
	mubox2 <- apply(mubox,1,sum)/(n-number_kernel_data);
	sub1_delta_matrix <- sapply(number_h,dmm,prod2);
	sub2_delta_matrix <- sub1_delta_matrix[1:prod2,1:hnumber];
	dim(sub2_delta_matrix) <- c(prod2*hnumber,1);
	sub3_delta_matrix <- sub1_delta_matrix[(prod2+1):(2*prod2),1:hnumber];
	dim(sub3_delta_matrix) <- c(prod2*hnumber,1);
	sub4_delta_matrix <- c(sub2_delta_matrix,sub3_delta_matrix);
	dim(sub4_delta_matrix) <- c((prod2*hnumber),2);
	delta_matrix <- matrix(apply(sub4_delta_matrix,1,deltaboxm,fhat_box,scheffe_matrix,mubox2,xgrid),nrow=hnumber,ncol=prod2,byrow=T);
	deltamax <- apply(delta_matrix,1,max);
	return(h[which.min(deltamax)]);
}
