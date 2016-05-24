##' Discounting masses
##'
##' Discount masses using  given factors 
##'
##' @export
##' @param MassIn Matrix with \eqn{nb} columns and \eqn{2^n} rows. Parameter \eqn{n} is the number of elements in the discernment frame and \eqn{nb} is the number of experts. Each column is a bba. If there is only one bba, the input could be a vector with length \eqn{2^n}.
##' @param alpha Discounting factor. A number or a vector with length of \code{ncol(MassIn)}. If it is a number, all the bbas will be discounted using the same factor. If it is a vector with length \code{ncol{MassIn}}, the bbas will be discounted using the corresponding factor. 
##' @return Mass matrix with the same dimension as MassIn. The discounted masses, each column is a piece of mass. If the input is a vector, the output is also a vector.
##' @examples
##' ## The conflict table for two experts in a discernment frame with three elements
##' m1=c(0,0.4, 0.1, 0.2, 0.2, 0, 0, 0.1);
##' m2=c(0,0.2, 0.3, 0.1, 0.1, 0, 0.2, 0.1);
##' discounting(m1,0.95)
##' # if only one factor is given, all the masses are discounted using the same factor
##' discounting(cbind(m1,m2),0.95)
##' # if the factor vector is given, the masses are discounted using the corresponding factor
##' discounting(cbind(m1,m2),c(0.95,0.9))
discounting <- function(MassIn, alpha){

# MassIn: nb of masses = nb column
#         1 column=1 mass
    
	if(!is.matrix(MassIn)){
	  MassIn = matrix(MassIn, ,1)
	}
	# number of focal elements
    mm = nrow(MassIn);
	# number of masses
	nn = ncol(MassIn);
    
	natoms = round(log2(mm)); 		
	if(2^natoms != mm || mm == 2){
	  stop('The number of focal element should be 2^n (n>1), with n the number of elements in the discernment frame\n')
	}
    if(length(alpha) == 1){
	   alpha = rep(alpha, nn);
	}
    
    if(length(alpha) == nn){
		alpha_mat = t(matrix(alpha, nn, mm)) 
		Mass = alpha_mat * MassIn;
		Mass[mm, ] = 1 - colSums(matrix(Mass[-mm, ], ,nn))
		if(ncol(Mass) == 1) Mass = as.vector(Mass)
		return(Mass)
	}else{
		stop('Accident: in discounting the size of alpha is incorrect\n')
	}


}
