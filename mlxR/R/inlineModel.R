#' inline model
#' 
#' @param str model
#' @param filename where to write the temporary model
#' @export
inlineModel <- function(str,filename=NULL){
  if(is.null(filename))
  {    
    # uuidd<-uuid()
    # filename <-paste0("temp_model-",uuidd,".txt")
    filename <-paste0("temp_model.txt")
  }
  write(str,filename)
  return(filename)
}
