#' @export
getArgs <- function(MLCall){
  # Replace the function by caller
  bracketLocation <- gregexpr(pattern ='\\(',MLCall)[[1]][1]
  arguments       <- substr(MLCall,bracketLocation,nchar(MLCall))
  functionCall    <- paste("caller",arguments,sep="")
  
  # Get all arguments
  args <- eval(parse(text=functionCall))
  return(args)
}

#' @export
caller <- function(...){
  get_args()
}

#' @export
get_args <- function(){
  as.list( sys.call(sys.parent()) )[-1]
}