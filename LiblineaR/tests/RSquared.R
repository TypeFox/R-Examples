RSquared<-function(pred,true){
	# Pierre Gramme
	# Jan 2015
	# Compute the R-squared coefficient of determination
	
	ss_tot = sum((true-mean(true)) ^ 2)
	ss_res = sum((true-pred) ^ 2)
  return(1 - ss_res/ss_tot)
}