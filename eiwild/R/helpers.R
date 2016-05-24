###############################################################################

# str-function without attributes


str0 <- function(...)
  str(..., give.attr=FALSE)

###############################################################################
#' alphaCheck function
#' 
#' Given starting values for the alpha values of the dirichlet distribution of 
#' the second level of the ecological inference model \code{alphaCheck} checks 
#' if parameters are correctly specified
#' 
#' @param mat matrix with alpha values
#' @param r number of rows in the RxC-table
#' @param c number of cols in the RxC-table
#' 
#' @return 
#' alphamatrix if everything is okay. \code{Error}-message if something fails.
#' 
#' @examples
#' \dontrun{
#' 
#' # right alpha matrix
#' alphaRight <- matrix(1:9, 3,3)
#' alphaCheck(alphaRight, 3, 3)
#' # return is right alpha matrix
#' 
#' # wrong alpha matrix
#' alphaWrong <- matrix(0:8, 3,3)
#' alphaCheck(alphaWrong, 3,3)
#' alphaCheck(alphaRight, 3,4)
#' }
#' 
#' @export

alphaCheck <- function(mat, r,c){
  if(ncol(mat)!=c | nrow(mat)!=r)
    stop("starting values of alpha must have ", c, " columns and ", r, " rows!", call.=FALSE)
  if(any(mat<=0))
    stop("starting values of alpha have to be greater than 0!", call.=FALSE)
  return(mat)
}

###############################################################################
#' betaCheck function
#' 
#' Given starting values for the beta values on the first level of the ecological
#' inference model, \code{betaCheck} checks if parameters are correctly specified
#' 
#' @param arr array with beta values. 1st dimension: rows, 2nd dimension columns,
#'          3rd dimension precincts
#' @param r number of rows in the RxC-table
#' @param c number of cols in the RxC-table
#' @param prec number of precincts
#' 
#' @return 
#' betaarry if everything is okay. \code{Error}-message if something fails.
#' 
#' @examples
#' \dontrun{
#' # right beta array
#' beta1 <- rep(c(0,0.25,0.75), each=3)
#' beta2 <- rep(beta1, 3)
#' betaRight <- array(beta2, dim=c(3,3,3))
#' betaCheck(betaRight, 3,3,3)
#' 
#' # wrong beta array
#' betaWrong <- array(1:27, dim=c(3,3,3))
#' betaCheck(betaWrong)
#' }
#' 
#' @export

betaCheck <- function(arr, r,c,prec){
  if(!all(dim(arr) == c(r,c,prec)))
    stop("starting values of beta be an array with dimensions ",
         paste( c(r,c,prec),collapse=", "),"!", call.=FALSE)
  
  for(rr in 1:r)
    if(!isTRUE(all.equal( colSums(arr[rr,,]), rep(1,prec) )))
      stop("starting values of beta must have sum of 1 in every row! \n row: ",rr, call.=FALSE)
  
  if(any(arr<0))
    stop("starting values of beta have to be positiv numbers! \n districts: ",
         paste(ceiling(which(arr <0)/(r*c)),collapse=", "), call.=FALSE)
  return(arr)
}


###############################################################################
#' calculate mfrow
#' 
#' @param x length of certain vector
#' 

calcmfrow<- function(x){
   temp <- sqrt(x)
   
   if((temp %% floor(temp))==0){
      return(c(temp,temp))
   } else if((temp %% floor(temp))<0.5){
      return(c(floor(temp),ceiling(temp)))
   } else {
      return(c(ceiling(temp),ceiling(temp)))
   }  
}


###############################################################################
#' calculate rolling mean up to each iteration
#' 
#' @param x vector
#' 

rollMean<- function(x){
   sapply(1:length(x), function(j) sum(x[1:j])/length((x[1:j])) )
}


###############################################################################
# printing eiwild object


#' @S3method print eiwild
print.eiwild <- function(x, ...){
   str0(x)
}

###############################################################################
#' calculating cell number
#' 
#' @description 
#' Calculates Cell numbers for whichCell-parameters in eiwild-plot function
#' 
#' @param rc integer vector with 2 elements giving the dimension of the ecological inference table
#' @param whichRow integer naming the specified row number
#' @param whichCol integer naming the specified column number
#' 
#' @return
#' returns vector of innerCells
#' 
#' @export
#' 

calcWhichCell <- function(rc, whichRow=NULL, whichCol=NULL){
   if(is.null(whichRow) & is.null(whichCol))
      stop("\"whichRow\" or \"whichCol\" have to be specified!", call.=FALSE)
   if(!is.null(whichRow) & !is.null(whichCol))
      stop("Only one of parameters \"whichRow\" or \"whichCol\" must be specified!", call.=FALSE)
   
   if(!is.null(whichRow)){
      if(whichRow > rc[1])
         stop("\"whichRow\" is bigger than dimensions in \"rc\"!", call.=FALSE)
      return(seq(whichRow, rc[1]*rc[2], by=rc[1]))
   }
   
   if(!is.null(whichCol)){
      if(whichCol > rc[2])
         stop("\"whichCol\" is bigger than dimensions in \"rc\"!", call.=FALSE)
      return(c(1:rc[1]) + rc[1]*(whichCol-1))
   }
}








