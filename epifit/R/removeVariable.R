##' Remove variables containing NaNs or NAs from data.frames.
##'
##' In reading data, incomplete variables are sometimes included, and only NAs or NaNs are included in some variables.
##' This function removes such variables in the data.frame. In addition to NAs or NaNs, variables which contain specified keyword can also be removed. 
##' @param data a data.frame from which variables are removed.
##' @param string a flag character specifying variables to be removed.
##' @param na.prop a numeric value specifying NA threshold proportion for removing variables. The default threashold is 1.0, meaning that variables including only NAs are removed.
##' @param nan.prop a numeric value specifying NaN threshold proportion for removing variables. The default threashold is 1.0, meaning that variables including only NaNs are removed.
##' @return a data.frame from which some variables are removed.
##' @examples
##' df <- data.frame(imcomp=rep(NA, 10), imcomp2=rep(NaN, 10), cov1=rnorm(10), NO_USE=rnorm(10))
##' df2 <- removeVariable(df, string="NO_USE")
##' str(df)
##' str(df2)
##' @export
removeVariable <- function(data=NULL, string=NA, na.prop=1.0, nan.prop=1.0){

  if(is.null(data)||!is.data.frame(data))
    stop("data is not specified or not data.frame")

  n <- nrow(data)
  varname <- colnames(data)
  idx <- rep(TRUE, length(varname))

  for(i in 1:length(varname)){
    if(length(grep(string, varname[i])) > 0 ||
       sum(as.integer(is.na(data[,i]))) >= na.prop*n ||
       sum(as.integer(is.nan(data[,i]))) >= nan.prop*n){
      idx[i] = FALSE
    }
  }
  return(data[,idx])
}
