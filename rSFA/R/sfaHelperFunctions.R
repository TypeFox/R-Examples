###################################################################################
#' Helper Function of SFA.
#'
#' @param range Range
#'
#' @return numeric vector \code{int} 
#' @keywords internal
###################################################################################
# not exported
sfaGetIntRange <- function (range){
	if(length(range)==1){
		int=1:range;
	}else{
		int=range[1]:range[2]
	}
	return(int)
}

###################################################################################
#' Custom Repeater Function
#'
#' Faster than customRepmat in matlab package, in certain cases.
#'
#' @param a to be repeated
#' @param n repeat by
#'
#' @return Returns the repeated dataset
#' @keywords internal
###################################################################################
#not exported
customRep<-function(a,n){  #faster than customRepmat for dimension (x,1)
	res<-kronecker(array(1, c(n,1)), a)
	return(res)
}

###################################################################################
#' Custom repmat Function
#'
#' R version of the matlab function repmat, repeating a by m x n
#'
#' @param a to be repeated
#' @param n repeater parameter
#' @param m repeater parameter
#'
#'
#' @return Returns the repeated dataset
#' @export
#' @keywords internal
###################################################################################
customRepmat <- function(a,n,m) {
	res<-kronecker(matrix(1,n,m),a)
	return(res)
}

###################################################################################
#' Custom Size Function.
#'
#' custom R version of matlabs size function. Calls length for vectors, or else calls dim.
#'
#' @param x object to be checked for size
#' @param i 1, 2 or NULL. Defines if both or which size information should be returned.
#'
#' @return a vector if i is null, an integer if i is 1 or 2.
#' @export
#' @keywords internal
###################################################################################
customSize<-function(x,i=NULL){	
	if((class(x)=="numeric")||(class(x)=="vector")||(class(x)=="integer")){
		res=c(1,length(x))
	}
	else{
		res=dim(x)
	}
	if(is.null(i))return(res)
	if(i==1)return(res[1])
	if(i==2)return(res[2])
}

###################################################################################
#' Check Condition of a matrix for SFA
#'
#' Creates warnings with recommendations for different settings, if given matrix is ill-conditioned.
#'
#' @param matr matrix to be checked
#' @param datatype string to identify "input" or "expanded" data
#' @keywords internal
###################################################################################
#not exported
sfaCheckCondition <- function(matr, datatype){
	ev=eigen(matr)
	ev=sort(ev$values) #TODO: maybe allready sorted (?)
	cn = abs(ev[length(ev)]/ev[1]);        # Condition number of B					
	if (cn > 1e10){
		if (datatype=="input"){
			warning("It is recommended to use sfaList$ppType = PCA2. 
			Covariance matrix of ", datatype," data is ill-conditioned. 
			CN=",paste(cn," "))
		}
		if (datatype=="expanded"){
			warning("It is recommended to use method=SVDSFA in SFA-step.
			Covariance matrix of ", datatype," data is ill-conditioned. 
			CN=",paste(cn," "))
		}
	}
}

###################################################################################
#' Backslash operator.
#'
#' Reproduce what MATLAB's backslash operator can do, using qr() and qr.coef().
#'
#' @param X X matrix
#' @param Y Y vector
#'
#' @return Returns coefficients 
#' @export
#' @keywords internal
###################################################################################
sfaBSh<-function(X,Y){
	X<-qr(X)
	qr.coef(X,Y)
}
