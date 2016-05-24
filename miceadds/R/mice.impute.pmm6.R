mice.impute.pmm6 <- function (y, ry, x, donors=3 , noise = 10^5 , ridge = 10^(-5) , ...){
    x <- cbind(1, as.matrix(x))
	# dummy response indicator 
	ry01 <- 1*ry
	# sample regression coefficients
	coefu <- stats::rnorm( ncol(x) )
	# create donor samples
	donorsample <- base::sample( 1:(2*donors) , sum( ! ry) , replace=TRUE) - donors
	donorsample <- ifelse( donorsample <= 0 , donorsample - 1 , donorsample )
	# Rcpp code
	# SEXP ma_pmm6_C( SEXP y_, SEXP ry01_, SEXP x_, SEXP ridge_, 
	#          SEXP coefu_, SEXP donorsample_) ;	
	imp <- .Call("ma_pmm6_C", 
				y_=y ,  ry01_=ry01 ,  x_=x,  ridge_=ridge ,  coefu_=coefu ,  donorsample_=donorsample	,
				PACKAGE = "miceadds")	
    return(imp)
	}
