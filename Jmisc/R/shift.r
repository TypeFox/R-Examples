#' Repeat a vector by row 
#' @name shift
#' @aliases shift
#' @title shift a vector by \code{shift_by} unit
#' @param x a vector
#' @param shift_by number of shift 
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @examples    
#' d<-data.frame(x=1:15) 
#' #generate lead variable
#' d$df_lead2<-shift(d$x,2)
#' #generate lag variable
#' d$df_lag2<-shift(d$x,-2)
shift<-function(x,shift_by){
	stopifnot(is.numeric(shift_by))
	stopifnot(is.numeric(x))

	if (length(shift_by)>1)
		return(sapply(shift_by,shift, x=x))

	out<-NULL
	abs_shift_by=abs(shift_by)
	if (shift_by > 0 )
		out<-c(tail(x,-abs_shift_by),rep(NA,abs_shift_by))
	else if (shift_by < 0 )
		out<-c(rep(NA,abs_shift_by), head(x,-abs_shift_by))
	else 
		out<-x
	out
}

#test 
# library(testthat)
# expect_that(shift(1:10,2),is_identical_to(c(3:10,NA,NA)))
# expect_that(shift(1:10,-2), is_identical_to(c(NA,NA,1:8)))
# expect_that(shift(1:10,0), is_identical_to(1:10))
# expect_that(shift(1:10,0), is_identical_to(1:10))
# expect_that(shift(1:10,1:2), is_identical_to(cbind(c(2:10,NA),c(3:10,NA,NA))))

