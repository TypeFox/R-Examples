MAP <-
function(mat, tie=c("random", "standard")){
	tie<- match.arg(tie)
	maxmat<- apply(mat, 1, max)
	M<-1:ncol(mat)
	map1<- rep(0, length(maxmat))
	
	if(tie=="random"){
	dummy<- mat == maxmat
	sel<- rowSums(dummy)>1
	if(any(sel==TRUE)) print("Some observations mapped to more than one group - these will be chosen at random.")	
	
	for(ind in (1: nrow(mat))[!sel]) map1[ind]<-M[mat[ind,]==maxmat[ind]]
	
	for(ind in (1: nrow(mat))[sel]) map1[ind]<- sample(M, 1, prob= dummy[ind,])
	}else { 
		for(ind in 1: nrow(mat)) map1[ind]<-M[mat[ind,]==maxmat[ind]]
 			}
		map1
	}
