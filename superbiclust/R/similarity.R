similarity <- function (x, index="jaccard",type="rows")
{
	options <- c("jaccard", "sorensen", "ochiai","kulczynski","sensitivity","specificity")
	if(any(!index %in% options)) {
		stop("index must be one of jaccard, sorensen, ochiai,kulczynski,sensitivity, or specificity")
	}		
	output <- switch(index,
			jaccard = jaccardMat(x, x,type), 
			sorensen = sorensenMat(x,x,type), 
			ochiai = ochiaiMat(x,x,type),
			kulczynski = kulczynskiMat(x,x,type),
			sensitivity = sensitivityMat(x,x,type),
			specificity = specificityMat(x,x,type))
	
	class(output)<- "similarity"
	output
}

jaccardMat <- function(x,y,type=c("rows","cols","both")){
	pA = x@Number
	pB = y@Number
	bicArows = t(x@GenesMembership)
	bicAcols =x@ColumnMembership
	
	bicBrows = y@GenesMembership
	bicBcols = y@ColumnMembership
	
	jamat <- matrix(0,pA,pB)	
	switch(type,rows={
				#option 1 :=rows
				bicArows = x@GenesMembership
				for(i in 1:pB){
					overlapABi <- bicArows+bicBrows[,i]
					AuB <- apply(overlapABi,2, function(x) sum(x>0)) #union	
					AnB <- apply(overlapABi,2, function(x) sum(x>1)) #intersection
					jamat[,i] <- AnB/AuB 
				}	},
			cols={
				#option 2 :=cols
				for(i in 1:pB){
					overlapABi <- t(bicAcols)+bicBcols[i,]
					AuB <- apply(overlapABi,2, function(x) sum(x>0)) #union	
					AnB <- apply(overlapABi,2, function(x) sum(x>1)) #intersection
					jamat[,i] <- AnB/AuB 
				}},
			both={
				#option 3 :=both
				n <- length(bicArows[1,]) #number of rows in data 
				l <- length(bicAcols[1,]) #number of columns in data
				u <- n*l
				for (i in 1:pA) {
					biclA <- tcrossprod(bicArows[i,],bicAcols[i,])
					apos <- biclA > 0
					sizeBiclA <- sum(apos) #size of bicluster A
					for (j in 1:pB) {
						biclB <- tcrossprod(bicBrows[,j],bicBcols[j,])
						bpos <- biclB > 0
						sizeBiclB <- sum(bpos)
						
						biclAB <- biclA + biclB
						abpos <- biclAB > 0
						sab <- sum(abpos)
						
						if (sab>0) {
							jamat[i,j] <- (sizeBiclA + sizeBiclB)/sab - 1
						} else {
							jamat[i,j] <- 0
						}
						
					}
				}})
	return(jamat)
}

sorensenMat <- function(x,y,type=c("rows","cols","both")){
	pA = x@Number
	pB = y@Number
	bicArows = t(x@GenesMembership)
	bicAcols =x@ColumnMembership
	
	bicBrows = y@GenesMembership
	bicBcols = y@ColumnMembership
	
	n <- length(bicArows[1,]) #number of rows in data 
	l <- length(bicAcols[1,]) #number of columns in data
	u <- n*l
	somat <- matrix(0,pA,pB)	
	switch(type,rows={
				#option 1 :=rows
				bicArows = x@GenesMembership
				NumberA <- colSums(bicArows)
				NumberB <- colSums(bicBrows)
				for(i in 1:pB){
					overlapABi <- bicArows+bicBrows[,i]
					AnB <- apply(overlapABi,2, function(x) sum(x>1)) #intersection
					somat[,i] <- 2*AnB/(NumberA+NumberB[i]) 
				}	},
			cols={
				#option 2 :=cols
				NumberA <- rowSums(bicAcols)
				NumberB <- rowSums(bicBcols)
				for(i in 1:pB){
					overlapABi <- t(bicAcols)+bicBcols[i,]
					AnB <- apply(overlapABi,2, function(x) sum(x>1)) #intersection
					somat[,i] <- 2*AnB/(NumberA+NumberB[i]) 
				}},
			both={
				#option 3 :=both
				
				for (i in 1:pA) {
					biclA <- tcrossprod(bicArows[i,],bicAcols[i,])
					apos <- biclA > 0
					sizeBiclA <- sum(apos) #size of bicluster A
					for (j in 1:pB) {
						biclB <- tcrossprod(bicBrows[,j],bicBcols[j,])
						bpos <- biclB > 0
						sizeBiclB <- sum(bpos)
						
						biclAB <- biclA + biclB
						abpos <- biclAB > 0
						sab <- sum(abpos)
						
						if (sab>0) {
							somat[i,j] <- 2.0-2.0*sab/(sizeBiclA+sizeBiclB)
						} else {
							somat[i,j] <- 0
						}
						
					}
				}})
	return(somat)	
}

kulczynskiMat <- function(x,y,type=c("rows","cols","both")){
	pA = x@Number
	pB = y@Number
	bicArows = t(x@GenesMembership)
	bicAcols =x@ColumnMembership
	
	bicBrows = y@GenesMembership
	bicBcols = y@ColumnMembership
	
	n <- length(bicArows[1,]) #number of rows in data 
	l <- length(bicAcols[1,]) #number of columns in data
	u <- n*l
	kumat <- matrix(0,pA,pB)	
	switch(type,rows={
				#option 1 :=rows
				bicArows = x@GenesMembership
				NumberA <- colSums(bicArows)
				NumberB <- colSums(bicBrows)
				for(i in 1:pB){
					overlapABi <- bicArows+bicBrows[,i]
					AnB <- apply(overlapABi,2, function(x) sum(x>1)) #intersection
					kumat[,i] <- 0.5*AnB*(NumberA+NumberB[i])/(NumberA*NumberB[i])
				}	},
			cols={
				#option 2 :=cols
				NumberA <- rowSums(bicAcols)
				NumberB <- rowSums(bicBcols)
				for(i in 1:pB){
					overlapABi <- t(bicAcols)+bicBcols[i,]
					AnB <- apply(overlapABi,2, function(x) sum(x>1)) #intersection
					kumat[,i] <- 0.5*AnB*(NumberA+NumberB[i])/(NumberA*NumberB[i]) 
				}},
			both={
				#option 3 :=both
				for (i in 1:pA) {
					biclA <- tcrossprod(bicArows[i,],bicAcols[i,])
					apos <- biclA > 0
					sizeBiclA <- sum(apos) #size of bicluster A
					for (j in 1:pB) {
						biclB <- tcrossprod(bicBrows[,j],bicBcols[j,])
						bpos <- biclB > 0
						sizeBiclB <- sum(bpos)
						
						biclAB <- biclA + biclB
						abpos <- biclAB > 0
						sab <- sum(abpos)
						
						if ((sizeBiclA>0)&&(sizeBiclB>0))
						{
							kumat[i,j] <- 1.0+0.5*( (sizeBiclA-sab)/sizeBiclB + (sizeBiclB -sab)/sizeBiclA )				
						}else {
							kumat[i,j] <- 0
						}			
					}
				}})
	return(kumat)	
}

ochiaiMat <- function(x,y,type=c("rows","cols","both")){
	pA = x@Number
	pB = y@Number
	bicArows = t(x@GenesMembership)
	bicAcols =x@ColumnMembership
	
	bicBrows = y@GenesMembership
	bicBcols = y@ColumnMembership
	
	n <- length(bicArows[1,]) #number of rows in data 
	l <- length(bicAcols[1,]) #number of columns in data
	u <- n*l
	ocmat <- matrix(0,pA,pB)
	switch(type,rows={
				#option 1 :=rows
				bicArows = x@GenesMembership
				NumberA <- colSums(bicArows)
				NumberB <- colSums(bicBrows)
				for(i in 1:pB){
					overlapABi <- bicArows+bicBrows[,i]
					AnB <- apply(overlapABi,2, function(x) sum(x>1)) #intersection
					ocmat[,i] <- AnB/sqrt(NumberA*NumberB[i])
				}	},
			cols={
				#option 2 :=cols
				NumberA <- rowSums(bicAcols)
				NumberB <- rowSums(bicBcols)
				for(i in 1:pB){
					overlapABi <- t(bicAcols)+bicBcols[i,]
					AnB <- apply(overlapABi,2, function(x) sum(x>1)) #intersection
					ocmat[,i] <- AnB/sqrt(NumberA*NumberB[i]) 
				}},
			both={
				#option 3 :=both
				for (i in 1:pA) {
					biclA <- tcrossprod(bicArows[i,],bicAcols[i,])
					apos <- biclA > 0
					sizeBiclA <- sum(apos) #size of bicluster A
					for (j in 1:pB) {
						biclB <- tcrossprod(bicBrows[,j],bicBcols[j,])
						bpos <- biclB > 0
						sizeBiclB <- sum(bpos)
						
						biclAB <- biclA + biclB
						abpos <- biclAB > 0
						sizeBiclAB <- sum(abpos)
						
						if ((sizeBiclA>0)&&(sizeBiclB>0))
						{
							ocmat[i,j] <- (sizeBiclA+sizeBiclB-sizeBiclAB)/sqrt(sizeBiclB*sizeBiclA)				
						}else {
							ocmat[i,j] <- 0
						}			
					}
				}})
	return(ocmat)	
}

sensitivityMat <- function(x,y,type=c("rows","cols","both")){
	pA = x@Number
	pB = y@Number
	bicArows = t(x@GenesMembership)
	bicAcols =x@ColumnMembership
	
	bicBrows = y@GenesMembership
	bicBcols = y@ColumnMembership
	
	n <- length(bicArows[1,]) #number of rows in data 
	l <- length(bicAcols[1,]) #number of columns in data
	u <- n*l
	senmat <- matrix(0,pA,pB)
	switch(type,rows={
				#option 1 :=rows
				bicArows = x@GenesMembership
				NumberA <- colSums(bicArows)
				NumberB <- colSums(bicBrows)
				for(i in 1:pB){
					overlapABi <- bicArows+bicBrows[,i]
					AnB <- apply(overlapABi,2, function(x) sum(x>1)) #intersection
					senmat[,i] <- AnB/NumberA
				}	},
			cols={
				#option 2 :=cols
				NumberA <- rowSums(bicAcols)
				NumberB <- rowSums(bicBcols)
				for(i in 1:pB){
					overlapABi <- t(bicAcols)+bicBcols[i,]
					AnB <- apply(overlapABi,2, function(x) sum(x>1)) #intersection
					senmat[,i] <- AnB/NumberA
				}},
			both={
				#option 3 :=both
				for (i in 1:pA) {
					biclA <- tcrossprod(bicArows[i,],bicAcols[i,])
					apos <- biclA > 0
					sizeBiclA <- sum(apos) #size of bicluster A
					for (j in 1:pB) {
						biclB <- tcrossprod(bicBrows[,j],bicBcols[j,])
						bpos <- biclB > 0
						sizeBiclB <- sum(bpos)
						
						biclAB <- biclA + biclB
						abpos <- biclAB > 0
						sizeBiclAB <- sum(abpos)
						
						if ((sizeBiclA>0)&&(sizeBiclB>0))
						{
							senmat[i,j] <- 1+ (sizeBiclB-sizeBiclAB)/sizeBiclA
						}else {
							senmat[i,j] <- 0
						}			
					}
				}})
	return(senmat)	
}

specificityMat <- function(x,y,type=c("rows","cols","both")){
	pA = x@Number
	pB = y@Number
	bicArows = t(x@GenesMembership)
	bicAcols =x@ColumnMembership
	
	bicBrows = y@GenesMembership
	bicBcols = y@ColumnMembership
	
	n <- length(bicArows[1,]) #number of rows in data 
	l <- length(bicAcols[1,]) #number of columns in data
	u <- n*l
	spemat <- matrix(0,pA,pB)
	switch(type,rows={
				#option 1 :=rows
				bicArows = x@GenesMembership
				NumberA <- colSums(bicArows)
				NumberB <- colSums(bicBrows)
				for(i in 1:pB){
					overlapABi <- bicArows+bicBrows[,i]
					AnB <- apply(overlapABi,2, function(x) sum(x>1)) #intersection
					spemat[,i] <- AnB/NumberB[i]
				}	},
			cols={
				#option 2 :=cols
				NumberA <- rowSums(bicAcols)
				NumberB <- rowSums(bicBcols)
				for(i in 1:pB){
					overlapABi <- t(bicAcols)+bicBcols[i,]
					AnB <- apply(overlapABi,2, function(x) sum(x>1)) #intersection
					spemat[,i] <- AnB/NumberB[i]
				}},
			both={
				for (i in 1:pA) {
					biclA <- tcrossprod(bicArows[i,],bicAcols[i,])
					apos <- biclA > 0
					sizeBiclA <- sum(apos) #size of bicluster A
					for (j in 1:pB) {
						biclB <- tcrossprod(bicBrows[,j],bicBcols[j,])
						bpos <- biclB > 0
						sizeBiclB <- sum(bpos)
						
						biclAB <- biclA + biclB
						abpos <- biclAB > 0
						sizeBiclAB <- sum(abpos)
						
						if ((sizeBiclA>0)&&(sizeBiclB>0))
						{
							spemat[i,j] <- 1+ (sizeBiclA-sizeBiclAB)/sizeBiclB
						}else {
							spemat[i,j] <- 0
						}			
					}
				}})
	return(spemat)	
}