low.anova <-
function(x,y,index,trace=NULL,...){
  if(is.null(trace)) trace <- FALSE

  n <- nrow(index)
  c <- ncol(index)
  if(c<2)stop("The order of interaction should be greater than 2")
  Max.R2 <- numeric()

  for(i in 1:n){
	if(trace)cat("Stage 2, Step:",i,"/",n,"\n")

	X1 <- x[,index[i,]]
	ok <- complete.cases(X1) ## Missing 3/3/2014	
	colnames(X1)= NULL
	X1 <- data.frame(X=X1)
	X1 <- X1[ok,]
	y1 <- matrix(y,,1)
	y1 <- y1[ok,]
        
	nR <- nrow(X1)
	pR <- 2
	max.R2 <- -1

	for(j in 1:c){
		x1 <- X1[,j]
		x2 <- X1[,-j]
		x3 <- x2[,1]
		for(k in 1:(c-2))x3<-paste(x3,x2[,k+1],sep=":")
		x3 <- factor(x3)
		x1 <- factor(x1)

		reg <- lm(y1~x1+x3)
		B <- summary(reg)
		R2 <- B$r.squared
		R2 <- 1-(1-R2)*(nR-1)/(nR-pR-1)

		if(R2 > max.R2) max.R2 <- R2
	}
	Max.R2[i] <- max.R2

  }
  
  res <- cbind(Max.R2,index)									
  return(res)

}
