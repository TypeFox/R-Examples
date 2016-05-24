AUCRF <-
function(formula,data,k0=1,pdel=0.2,ranking=c("MDG","MDA"),...){

AUC.randomForest <-
function(rf,clase=1){
r <- rank(rf$votes[,as.character(clase)])
rd <- mean(r[rf$y==clase])
nd <- sum(rf$y==clase)
nnd <- length(rf$y)-nd
return((rd-nd/2-0.5)/nnd)
}

MDGRanking <-
function(formula,data,...){
fitRF <- randomForest(formula,data=data,...)
  mdgRanking <- sort(fitRF$importance[,"MeanDecreaseGini"],decreasing=TRUE)
  return(mdgRanking)
}

MDARanking <-
function(formula,data,...){
fitRF <- randomForest(formula,data=data,importance=TRUE,...)
  mdaRanking <- sort(fitRF$importance[,"MeanDecreaseAccuracy"],decreasing=TRUE)
  return(mdaRanking)
}

t <- 0	
cl <- match.call()
mf <- match("formula", names(cl), 0L)
y <- eval(eval(cl[[mf]])[[2]],data)
if(!is.factor(y) && length(levels(y))!=2) stop("Outcome must be a factor with two levels")
if(pdel<0 || pdel>=1) stop("pdel must be in the interval [0,1)")

ranking <- match.arg(ranking)
switch(ranking,
         MDG = {ranking <- MDGRanking(formula,data,...); ImpMes <- "MDG"},
         MDA = {ranking <- MDARanking(formula,data,...); ImpMes <- "MDA"},
         stop("Not valid ranking")
       )
        

mf <- match("formula", names(cl), 0L)
yname <- as.character(eval(cl[[mf]])[[2]])
  vars <- names(ranking)
  AUCcurve <- data.frame()
  auxThres <- 0
  auxMaxAUC <- 0
  k <- length(vars)
while(k>=k0){
	fitRF <- randomForest(formula,data=data[,c(yname,vars[1:k])],...)
	getAUC <- AUC.randomForest(fitRF)
	if(getAUC>=auxMaxAUC){
		auxMaxAUC <- getAUC
		auxThres <- auxMaxAUC-t
	}
	if(getAUC>=auxThres) RFopt <- fitRF
	AUCcurve <- rbind(c(k,getAUC),AUCcurve)
	k <- k-as.integer(pdel*k)-1
}
 
colnames(AUCcurve) <- c("k","AUC")
maxAUC <- max(AUCcurve$AUC)
opthreshold <- maxAUC-t
optimal <- AUCcurve[AUCcurve$AUC>=opthreshold,][1,]

objectList <- list()
objectList$call <- cl
objectList$data <- data
objectList$ranking <- ranking
objectList$Xopt <- names(ranking)[1:(optimal$k)]
objectList$"OOB-AUCopt" <- optimal$AUC
  objectList$Kopt <- optimal$k
objectList$AUCcurve <- AUCcurve
#objectList$AUCthres <- opthreshold
objectList$RFopt <- RFopt
objectList$ImpMeasure <- ImpMes
  class(objectList) <- "AUCRF" 
  return(objectList)
}

