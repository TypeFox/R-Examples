convertToProp <- function(y, T0=NULL, Ctrl=NULL){
  # If nothing specified: T0 is min(y), Ctrl is max(y)
  if( is.null(T0) ){
    T0 <- min(y, na.rm=TRUE)
  }
  if( is.null(Ctrl) ){
    Ctrl <- max(y, na.rm=TRUE)
  }
  
  ## CHECK T0 AND Ctrl TO ENSURE THEY ARE OF THE CORRECT FORMAT
  if( !(is.numeric(T0) & is.numeric(Ctrl)) ){
    stop("T0 and Ctrl must each be numeric")
  }
  if( length(T0) != 1 | length(Ctrl) != 1 ){
    stop("T0 and Ctrl must each be of length 1")
  }
  if( is.na(T0) | is.na(Ctrl) ){
    stop("Please provide non-NA T0 and Ctrl values.")
  }
  
  return((y-T0)/(Ctrl-T0))
}
