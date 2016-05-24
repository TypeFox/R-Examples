XandYmatrix.endogMNP <- function(selX, selY, outX, outY, base1=NULL, 
base2=NULL, extra=FALSE, verbose=verbose){
	selXnames <- colnames(selX)
	outXnames <- colnames(outX)
	rownames(selX) <- NULL
	rownames(outX) <- NULL
	selX <- data.matrix(selX)
	

	
	outX <- data.matrix(outX)
	Y1 <- as.factor(selY)
	lev1 <- levels(Y1)
	if(!is.null(base1))
	if (base1 %in% lev1){
		Y1 <- relevel(Y1, ref=base1)
		lev1 <- levels(Y1)}
		else{
			stop(paste("Error: 'base' does not exist in the response variable."))	}
	base1 <- lev1[1]
	counts1 <- table(Y1)
	if(any(counts1==0)){
		warning(paste("group(s)", paste(lev1[counts1 == 0], collapse = " "), "are empty"))
		Y1 <- factor(Y1, levels=lev1[counts1 > 0])
		lev1 <- lev1[counts1 > 0]
		}
	Y2 <- as.factor(outY)
	lev2 <- levels(Y2)
	if(!is.null(base2))
	if (base2 %in% lev2){
		Y2 <- relevel(Y2, ref=base2)
		lev2 <- levels(Y2)}
		else{
			stop(paste("Error: `base' does not exist in the response variable."))	}
	base2 <- lev2[1]
	counts2 <- table(Y2)	
	if(any(counts2==0)){
		warning(paste("group(s)", paste(lev2[counts2 == 0], collapse = " "), "are empty"))
		Y2 <- factor(Y2, levels=lev2[counts2 >0])
		lev2 <- lev2[counts2 > 0]
	}	
	p1 <- length(lev1)
	p2 <- length(lev2)
	Y1 <- as.matrix(unclass(Y1)-1, nc=1) 
	Y2 <- as.matrix(unclass(Y2)-1)

 naMat <- matrix(NA, nc=p1, nr=length(Y1))
	Y <- matrix(cbind(Y1, naMat), nc=p1+1)
whichNotObs <- NULL	
notObsCat <- FALSE

	for(i in 1:length(Y1)){
		if(!is.na(Y1[i])){
		   Y[i,(Y1[i]+2)] <- Y2[i]}}
	
	
 if(sum(apply(is.na(Y), 2, prod))>0){
	if(sum(apply(is.na(Y), 2, prod))>1){
	stop(paste("Error: Only one selection category can have no observed outcomes\n"))}
	else {whichNotObs <- which.max(colSums(is.na(Y))) -1
		Y <- Y[,colSums(!is.na(Y)>0)>0]
		notObsCat <- TRUE
		
	}}
		
			
	
	Y <- matrix(Y, nc=1, byrow=TRUE)	
	if(notObsCat){
		xMatrix <- matrix(0, nr=((p1-1+(p1-1)*(p2-1))*length(Y1)), nc=((length(selX[1,]))*(p1-1)+
																   (length(outX[1,]))*(p1-1)*(p2-1)))
	}
	else{
xMatrix <- matrix(0, nr=((p1-1+p1*(p2-1))*length(Y1)), nc=((length(selX[1,]))*(p1-1)+
														   (length(outX[1,]))*p1*(p2-1)))
	}

for(i in 1:length(Y1)){

		selVec <- matrix(selX[i,],nr=1)
		outVec <- matrix(outX[i,], nr=1)
	

	selPart <- kronecker(diag(1,(p1-1)), selVec)
	if(notObsCat) outPart <- kronecker(diag(1,(p1-1)*(p2-1)), outVec)
	
	else outPart <- kronecker(diag(1, p1*(p2-1)), outVec)
	
	zeroPadUpRight <- matrix(0, nr=dim(selPart)[1], nc=dim(outPart)[2])
	upperRows <- cbind(selPart,zeroPadUpRight)
	zeroPadLowLeft <- matrix(0, nr=dim(outPart)[1], nc=dim(selPart)[2])
	lowRows <- cbind(zeroPadLowLeft, outPart)
	Xi <- rbind(upperRows, lowRows)
	xMatrix[(1+(i-1)*dim(Xi)[1]):(i*dim(Xi)[1]),] <- Xi
	}
	
#	print(head(xMatrix))	
selCovNum = length(selVec)
outCovNum = length(outVec)	
xMatrix <- matrix(xMatrix, nc=1, byrow=TRUE)		
if(extra){
	return(list(Y=Y, X=xMatrix, lev1=lev1, p1=p1, base1=base1, lev2=lev2, p2=p2, base2=base2, 
				selCovNum=selCovNum, outCovNum=outCovNum, selXnames=selXnames, outXnames=outXnames,
		   notObsCat=notObsCat, whichNotObs=whichNotObs))}
	else{
		return(list(Y=Y, X=xMatrix))}
	}			