noba <-
function(x, batch) {

  if(any(is.na(x)))
	stop("Data contains missing values.")
  if(!(is.factor(batch) & all(levels(batch)==(1:length(levels(batch))))))
    stop("'batch' has to be of class 'factor' with levels '1','2',...")  
  if(!is.matrix(x))
    stop("'x' has to be of class 'matrix'.") 
  
  params <- list(xadj=x)
  params$nbatches <- length(unique(batch))
  params$batch <- batch  
  
  class(params) <- "noba"
   
  return(params)

}
