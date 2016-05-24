wavesum <-
function(x, populations,popA=NA, popB=NA, ml=NULL, type = "la8",t.factor=1, fullWT =FALSE){

	#check input arguments
  	ctmp <- class(x)
  	if (is.null(ctmp))  stop("object has no class")
 	if (ctmp != "adsig")  stop("object is not of class adsig")
	#ancestral populations
      if(sum(is.na(populations[[popA]]))>0 | sum(is.na(populations[[popB]]))>0) stop("Check pop1 pop2. Exiting....")

	signals <- x$signals
	n.ind <- dim(signals)[2]
	T <- dim(signals)[1]

	popnames <- names(populations)
	n.group <- length(popnames)


	result <- list()
	result$n.ind <- n.ind
	result$n.group <-n.group
	if(is.null(ml)){ml <- floor(log2(dim(signals)[1]))}

	#setting objects for store
	result$rv.ind <- array( dim=c(n.ind, ml),data=NA); row.names(result$rv.ind) <- colnames(signals)
	result$rv.group <- array( dim=c(n.group,ml),data=NA);row.names(result$rv.group) <- popnames
	result$threshold  <- array( dim=c(ml),data=NA)
	result$iv.group <- result$rv.group
	result$iv.ind <- result$rv.ind
	result$abs.ind <- array( dim=c(n.ind),data=NA);names(result$abs.ind)<- colnames(signals) #Individual level centers
	result$abs.group <- array( dim=c(n.group),data=NA);names(result$abs.group)<- popnames #population level centers
	result$pws.ind <- result$abs.ind
	result$pws.group <- result$abs.group

	if(fullWT==TRUE){
	result$wtmatrix <- array( dim=c(T,n.ind, ml),data=NA);
	dimnames(result$wtmatrix)<- list(c(1:T),colnames(signals),c(1:ml))
	result$wtmatrix.group <- array( dim=c(T,n.group, ml),data=NA)
	dimnames(result$wtmatrix.group)<- list(c(1:T),popnames,c(1:ml))}



	#USING WAVESLIM
	for(i in 1:n.ind){
	wd1 <- modwt(t(signals)[i,], wf = type, n.levels = ml, boundary = "periodic")
	#extract coefficients to a matrix
 	 for(j in 1:ml){
    		result$rv.ind[i,j] <- mean(abs(wd1[[j]])) 
		if(fullWT ==TRUE){result$wtmatrix[,i,j] <- (abs(wd1[[j]])^2)}}
		}

	#population averages
	for(g in 1:n.group){
	result$rv.group[g,] <- colMeans(result$rv.ind[populations[[popnames[g]]],])
	if(fullWT ==TRUE){result$wtmatrix.group[,g,] <- apply(result$wtmatrix[,populations[[popnames[g]]],],c(1,3),FUN=mean)}
	}	

	result$threshold <- pmax(result$rv.group[popA,],result$rv.group[popB,])*t.factor
	#apply threshold
		for(g in 1:n.group){
	result$iv.group[g,] <-  pmax((result$rv.group[g,] - result$threshold) , rep(0,ml))
	result$abs.group[g] <- sum(result$iv.group[g,]*c(1:ml))/sum(result$iv.group[g,])
	result$pws.group[g]<- mean((1:ml)[result$iv.group[g,] == max(result$iv.group[g,])])

	if(sum(result$iv.group[g,]) == 0){result$abs.group[g] <- 0
								result$pws.group[g] <- 0} }

	#individual level
	for(i in 1:n.ind){
	result$iv.ind[i,] <- pmax((result$rv.ind[i,]- result$thresh), rep(0,ml))
	result$abs.ind[i] <- sum(result$iv.ind[i,]*(1:ml))/sum(result$iv.ind[i,])
	result$pws.ind[i] <- mean((1:ml)[result$iv.ind[i,] == max(result$iv.ind[i,])])
	if(sum(result$iv.ind[i,]) == 0){result$abs.ind[i] <- 0;result$pws.ind[i] <- 0}		}
	
	return(result)
	}
