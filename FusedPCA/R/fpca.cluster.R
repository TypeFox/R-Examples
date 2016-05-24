fpca.cluster <- function(obj, K = 2, score = F){
	
	obj = obj$final.whole
	if(score == F){
		
		if(class(obj) == 'FPCA-RoE'){
			warning('This object is designed for score = T.')
			cluster.list = fpca.score.cluster(obj, K = K)	
		} 
		
		if(class(obj) == 'FPCA') cluster.list = fpca.nonscore.cluster(obj, K = K)
	
	}
	
	if(score == T){
		
		if(class(obj) == 'FPCA-RoE') cluster.list = fpca.score.cluster(obj, K = K)
		
		if(class(obj) == 'FPCA') stop('This object is designed for score = F.')

	}
		
	return(cluster.list)
	
}
