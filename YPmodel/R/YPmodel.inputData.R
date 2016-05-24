YPmodel.inputData <-
function(data=c(), ...)
{
	if(is.null(data)){
		stop(paste(fun.errorMessage('DataSet')))
	}

	if(class(data)=="character"){
	  data <- read.table(data)
	}

	(o <- order(data$V1,data$V2,data$V3))
	data <- data[o,]

	X <- as.numeric(data[,1])
	Delta <- as.numeric(data[,2])
	Z <- as.numeric(data[,3])
	n <- length(X)

	#Group data
#	temp <- matrix(c(1:n), nrow = n, ncol = 1)
#	Z1 <- temp[Z==0]
#	Z2 <- temp[Z==1]
#	X1 <- X[Z1]
#	X2 <- X[Z2]
#	n1 <- length(X1)
	n2 <- sum(Z)
	n1 <- n-n2
	glength <- list(Num1=n1,Num2=n2)
	GroupLength <- list(length=glength)

	X <- matrix(as.numeric(X), nrow = n, ncol = 1)
	Delta <- matrix(as.numeric(Delta), nrow = n, ncol = 1)
	Z <- matrix(as.numeric(Z), nrow = n, ncol = 1)
	Data <- list(X=X, Delta=Delta, Z=Z, length=n, GroupLength=GroupLength)

	Data

}
