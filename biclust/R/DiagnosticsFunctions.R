#' Bootstrapping of Columns Outside the Bicluster
#' The objective of the bootstrapping is to compare the structure within the bicluster 
#' and outside of the  bicluster
#' @param x original expression matrix
#' @param bicResult output of the class biclust
#' @param number number of a bicluster to diagnose
#' @param nResamplings number of resampling iterations (either bootstrap or permutation, cf. 'replace') 
#' @param replace should one use a bootstrapping (TRUE; default) or permutation (FALSE) based empirical distribution
#'   outputs the bootstrapping distributions of F scores and p-values from ANOVA table
#' @return 
#' @export

diagnoseColRow <- function(x, bicResult, number, nResamplings, replace = TRUE){
	
	coreBiclusterSamples <- which(bicResult@NumberxCol[number,])
	coreBiclusterGenes <- which(bicResult@RowxNumber[,number]) 
	
	B <- nResamplings
	fval.row <- fval.col <- rep(NA, nResamplings)
	N <- ncol(x)
	
	outCols <- c(1:N)[-coreBiclusterSamples] # columns outside bicluster
	xmat <- x[coreBiclusterGenes, ]
	xA <- as.vector((xmat[, coreBiclusterSamples]))
	
	roweff <- rep(c(1:length(coreBiclusterGenes)), length(coreBiclusterSamples))
	coleff <- sort(rep(c(1:length(coreBiclusterSamples)), length(coreBiclusterGenes) ))
	
	xx1 <- anova(aov(xA ~ as.factor(roweff) + as.factor(coleff)))
	
	fval.row.obs <- xx1[1, "F value"]
	fval.col.obs <- xx1[2, "F value"]
	
	for (b in 1:B){
		
		index <- sample(outCols, size = length(coreBiclusterSamples), replace=replace) #sampling without replacement from original data
		xB <- as.vector((xmat[, index])) # converting data from matrix to the vector form
		
		xx2 <- anova(aov(xB ~ as.factor(roweff) + as.factor(coleff)))
		fval.row[b] <- xx2[1,4]
		fval.col[b] <- xx2[2,4]
	}
	pval.row <- sum(fval.row > fval.row.obs)/(nResamplings + 1)
	pval.col <- sum(fval.col > fval.col.obs)/(nResamplings + 1)
	retval <- cbind(fval.row, fval.col)

	retval <- list(bootstrapFstats = retval, observedFstatRow = fval.row.obs, observedFstatCol = fval.col.obs,
					bootstrapPvalueRow = pval.row, bootstrapPvalueCol = pval.col)
	
	return(retval)
}
#' Compute Observed F Statistics and Asymptotic P Values for Gene And Sample Effect and Interaction within Bicluster 
#' @param x original expression matrix
#' @param bicResult - biclust object, output from biclustering algorithm
#' @param number - number of bicluster in an output
#' @return dataframe with two variables - F-values and p-values for the three effects 
#' @export
computeObservedFstat<- function(x, bicResult, number){ 
	coreBiclusterSamples <- which(bicResult@NumberxCol[number,])
	coreBiclusterGenes <- which(bicResult@RowxNumber[,number]) 
	xmat <- x[coreBiclusterGenes, ]
	xA <- as.vector((xmat[, coreBiclusterSamples]))
	
	roweff<-rep(c(1:length(coreBiclusterGenes)), length(coreBiclusterSamples))
	coleff<-sort(rep(c(1:length(coreBiclusterSamples)), length(coreBiclusterGenes) ))
	
	xx1 <- anova(aov(xA ~ as.factor(roweff) + as.factor(coleff)))
	
	fval.row.obs <- xx1[1, "F value"]
	fval.col.obs <- xx1[2, "F value"]
	
	pval.row.obs <- xx1[1, "Pr(>F)"]
	pval.col.obs <- xx1[2, "Pr(>F)"]
	
	xmat2 <- x[coreBiclusterGenes, coreBiclusterSamples]
	
	y.bar.j<-colMeans(xmat2)
	y.bar.i<-rowMeans(xmat2)
	x.vec<-as.vector(xmat2)
	y.bar <- mean(x.vec)
	B <- sum((y.bar.i - y.bar)^2)
	C <- sum((y.bar.j - y.bar)^2)
	A <- 0
	for(i in 1: length(coreBiclusterGenes )){
		for(j in 1: length(coreBiclusterSamples)){
			A <- A + xmat2[i,j]*(y.bar.i[i]-y.bar)*(y.bar.j[j]-y.bar)
		}
	}
	F.interaction <- (A^2)/(B*C)
	pval.Tukey.obs <- 1 - pf(q=F.interaction, df1=1, df2=(length(coreBiclusterGenes)*length(coreBiclusterSamples) -(length(coreBiclusterGenes)+length(coreBiclusterSamples))))
	retval <- data.frame(Fstat = c(fval.row.obs, fval.col.obs, F.interaction), 
			PValue = c(pval.row.obs, pval.col.obs, pval.Tukey.obs))
	rownames(retval) <- c("Row Effect", "Column Effect", "Tukey test")
	return(retval)
	
}
#' Plot result of Bootstrap distribution of Scores 
#' @param bootstrapOutput object as returned by the diagnoseColRow function 
#' @return no return value; a plot is drawn to the current device 
#' @export
diagnosticPlot <- function(bootstrapOutput){
  
  fval.row.obs <- bootstrapOutput[["observedFstatRow"]]
  fval.col.obs <- bootstrapOutput[["observedFstatCol"]]
  fval.row <- bootstrapOutput[["bootstrapFstats"]][,"fval.row"]
  fval.col <- bootstrapOutput[["bootstrapFstats"]][,"fval.col"]
  
  op <- par(mfrow = c(1, 2))
  dx1<-density(fval.row)
  hist(fval.row, nclass = 30, probability = TRUE, xlab = "F(A)", main = "row scores")
  lines(dx1$x, dx1$y)
  lines(rep(fval.row.obs, 2), c(0, 5), col="green", lwd = 4)
  
  
  dx1<-density(fval.col)
  hist(fval.col, nclass = 30, probability = TRUE, xlab = "F(B)", main = "column scores")
  lines(dx1$x, dx1$y)
  lines(rep(fval.col.obs, 2), c(0, 5), col = "green", lwd = 4)
  
  par(op)
  
}
#====Chia and Karuturi scores=====#
#' function computing scores as described in the paper of Chia and Karuturi (2010)
#' @param bicResult  - biclust object from biclust package
#' @param number - number of bicluster in the output for computing the scores
#' @return dataframe with 6 slots: T, B scores for within and outside bicluster, SB and TS scores 
ChiaKaruturi <- function(x, bicResult, number){
	
	biclCols <- which(bicResult@NumberxCol[number,])
	biclRows <- which(bicResult@RowxNumber[,number]) 
	xmat <- x[biclRows, biclCols]
	
	y.bar.j<-colMeans(xmat)
	y.bar.i<-rowMeans(xmat)
	x.vec<-as.vector(xmat)
	y.bar <- mean(x.vec)
	E <- 0
	for(i in 1: length(biclRows)){
		for(j in 1: length(biclCols)){
			E <- E + (xmat[i,j] - y.bar.i[i] - y.bar.j[j] + y.bar)^2
		}
	}
	E <- E/((ncol(xmat) - 1) * (nrow(xmat) - 1))
	
	T1 <- sum((y.bar.i)^2)/nrow(xmat) - E/ncol(xmat)
	B1 <- sum((y.bar.j)^2)/ncol(xmat) - E/nrow(xmat)
	
	
	xmat2 <- x[biclRows, -biclCols]
	
	y.bar.j2 <- colMeans(xmat2)
	y.bar.i2 <- rowMeans(xmat2)
	x.vec2 <- as.vector(xmat2)
	y.bar2 <- mean(x.vec2)
	E2 <- 0
	for(i in 1: nrow(xmat2) ){
		for(j in 1: ncol(xmat2) ){
			E2 <- E2 + (xmat2[i,j] - y.bar.i2[i] - y.bar.j2[j] + y.bar2)^2
		}
	}
	E2 <- E2/((ncol(xmat2) - 1) * (nrow(xmat2) - 1))
	T2 <- sum((y.bar.i2)^2)/nrow(xmat2) - E/ncol(xmat2)
	B2 <- sum((y.bar.j2)^2)/ncol(xmat2) - E/nrow(xmat2)
	
	a <- 0.01
	
	SB <- log(max((T1+a), (B1+a)) /max((T2+a), (B2+a) ))
	if(SB > 0 ) {TS <- log((T1+a)/(B1+a))} else  TS <- log((T2+a)/(B2+a))
	
	retvalue <- data.frame(Tscore1 = T1, Bscore1 = B1, Tscore2 = T2, Bscore2 = B2,
							SBscore = SB, TSscore = TS)
	return(retvalue)
}

