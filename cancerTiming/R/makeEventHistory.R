#returns matrix of the possible matrices

makeEventHistory<-function(type=c("gain","LOH"), copies=NULL,totalCopy=sum(copies),onlyIdentifiable=TRUE){
	type<-match.arg(type)
	
	 doCopies<-if(!is.null(copies)) TRUE else FALSE
	if(doCopies){
		if(length(copies)!=2) stop("copies should be a vector of integers of length 2")
		copies<-round(copies)
		copies<-sort(copies)
	}
	if(totalCopy>5 & !onlyIdentifiable) stop("Only totalCopy<=5 is currently implemented.")
	if(type=="gain" & totalCopy<=2) stop("Cannot have 'gain' and totalCopy<3")
	if(!onlyIdentifiable & type=="LOH" & totalCopy>2) stop("LOH is only an option for totalCopy <3.")
	if(onlyIdentifiable & type=="LOH" & totalCopy!=2) stop("Identifiable LOH is only an option for totalCopy=2.")
	if(onlyIdentifiable & totalCopy==1) stop("There are no identifiable histories for totalCopy=1.")
	makeMatrix<-function(x){matrix(x,ncol=totalCopy-1,nrow=totalCopy-1,byrow=T)}

	if(totalCopy>5){
		if(!onlyIdentifiable) stop("for totalCopy>5 can only return the single identifiable event history")
		permMat<-diag(totalCopy-1)[(totalCopy-1):1,] 
		rank1Adjust<-makeMatrix(0)
		rank1Adjust[1,]<-1:(totalCopy-1)
		AMat<-list(permMat+rank1Adjust)
	}
	else{
		#for now, manually program them in
		AMat<-list()
		# if(totalCopy==1){
		# 	if(type=="LOH") AMat<-
		# }
		if(totalCopy==2){
			if(type=="LOH") AMat<-c(AMat,"1MLineage"=list(matrix(c(0,2,1,0),ncol=2,nrow=2,byrow=T)))
		}
		if(totalCopy==3){
			if(type=="gain") AMat<-c(AMat,"1MLineage"=list(makeMatrix(c(1,3,1,0))))
		}
		if(totalCopy==4){
			if((doCopies & all(copies==c(1,3))) | !doCopies) AMat<-c(AMat,list("1MLineage x 2"=makeMatrix(c(1,2,4,0,1,0,1,0,0))))
			if(!onlyIdentifiable){
				if((doCopies & all(copies==c(2,2))) | !doCopies) AMat<-c(AMat,list("1MLineage,1PLineage"=makeMatrix(c(0, 2, 4,2, 1, 0,0,0,0))))
			}
		}
		if(totalCopy==5){
			if((doCopies & all(copies==c(1,4))) | !doCopies){
				AMat<-c(AMat,list("1MLineage x 3"=makeMatrix(c(1, 2, 3, 5,0,0,1,0,0,1,0,0,1,0,0,0))))
				if(!onlyIdentifiable){
					AMat<-c(AMat,list("1MLineage x 2,2MLineage"=makeMatrix(c(1, 1, 3, 5,0,2,1,0,0,0,0,0,1,0,0,0))))			
				}
			} 
			if((doCopies & all(copies==c(2,3))) | !doCopies){
				AMat<-c(AMat,list("1MLineage x 2,1PLineage"=makeMatrix(c(0 ,1, 3, 5,1,2,1,0,1,0,0,0,0,0,0,0))))			
				AMat<-c(AMat,list("1MLineage,1PLineage,1MLineage"=makeMatrix(c(0,1,3,5,1,2,1,0,1,0,0,0,0,0,0,0))))			
				AMat<-c(AMat,list("1PLineage,1MLineage x 2"=makeMatrix(c(0, 2, 3, 5,1,0,1,0,1,1,0,0,0,0,0,0))))			
			}

		}
	}

	return(AMat)
}