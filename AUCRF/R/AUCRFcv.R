AUCRFcv <-
function(x,nCV=5,M=20){

AUC.votes <-
function(votes,y=NULL,clase=1){
if(missing(y) || is.null(y)) y <- votes$y 
r <- rank(votes[,as.character(clase)])
rd <- mean(r[y==clase])
nd <- sum(y==clase)
nnd <- length(y)-nd
return((rd-nd/2-0.5)/nnd)
}

cl <- match.call()
switch(class(x),
  "AUCRF" = { 
      callRF <- x$call
      data <- x$data
      callRF$data <- as.name("newData")
      yname <- as.character(eval(x$call$formula)[[2]])
  },
stop("x must be a AUCRF object.")
)

cvAUC <- NULL
varnames <- colnames(data)[colnames(data)!=yname]
nSelect <- rep(0,length(varnames))
names(nSelect) <- varnames
for(m in 1:M){
  CV <- list()
  mpredict <- NULL
  indPermuted <- matrix(c(sample(rownames(data)),rep(NA,nCV-nrow(data)%%nCV)),ncol=nCV,byrow=TRUE) 
  for(k in 1:nCV){
  indTest <- indPermuted[,k]
  indTest <- indTest[!is.na(indTest)]
  indTrain <- rownames(data)[!(rownames(data) %in% indTest)]
  newData <- data[indTrain,]
    kaucRF <- eval(callRF)
      mpredict <- rbind(mpredict, predict(kaucRF$RFopt,newdata=data[indTest,],type="vote"))
      nSelect[kaucRF$Xopt] <-  nSelect[kaucRF$Xopt]+1
   }
   mvotes <- data.frame(y=data[,yname],mpredict[rownames(data),])
   class(mvotes) <- c("votes","data.frame")
   colnames(mvotes) <- c("y","0","1")
    cvAUC <- c(cvAUC, AUC.votes(mvotes))
  }
  
  if(class(x)=="AUCRF")
  objectList <- x
  else
    objectList <- list()
  
  objectList$cvAUC <- mean(cvAUC)
  objectList$Psel <- nSelect/(M*nCV)
  objectList$callcv <- cl
  class(objectList) <- c("AUCRFcv","AUCRF")
  return(objectList)
}

