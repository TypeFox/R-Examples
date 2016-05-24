#caSupplementalElementsPreProcessing <- function(SUP.DATA,hellinger=FALSE){
caSupplementalElementsPreProcessing <- function(SUP.DATA){

	#if(!hellinger){
	#	#sup.transform <- SUP.DATA/((rowSums(SUP.DATA))%*%matrix(1,1,ncol(SUP.DATA)))
		#return(SUP.DATA/((rowSums(SUP.DATA))%*%matrix(1,1,ncol(SUP.DATA))))
		return(rowNorms(SUP.DATA,'ca'))
	#}else{
		#sup.transform <- (SUP.DATA/repmat(rowSums(SUP.DATA),1,ncol(SUP.DATA)))^(1/2)
	#	return((SUP.DATA/repmat(rowSums(SUP.DATA),1,ncol(SUP.DATA)))^(1/2))
	#}
	#return(sup.transform)		
}