meancenter <-
function(x, batch) {

  if(any(is.na(x)))
	stop("Data contains missing values.")
  if(!(is.factor(batch) & all(levels(batch)==(1:length(levels(batch))))))
    stop("'batch' has to be of class 'factor' with levels '1','2',...")  
  if(!is.matrix(x))
    stop("'x' has to be of class 'matrix'.") 

   batches = levels(batch)
   nbatches = length(batches)
   
   if(nbatches > 1) {
     # Adjust for additive batch effects (Reference batch?!):
     adjustmentmod = lm(x~batch)
     design = model.matrix(~batch)
     adjustmentcoef = coef(adjustmentmod)
     xadj = x-design%*%adjustmentcoef
   }
   else {
	 adjustmentcoef <- colMeans(x)
     xadj = scale(x, center=adjustmentcoef, scale=FALSE)	 
   }
   
   params <- list(xadj=xadj)
   params$nbatches <- nbatches
   params$batch <- batch
   
   class(params) <- "meancenter"     
   
   return(params)
}
