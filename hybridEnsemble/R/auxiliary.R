.tuneMember <- function(call,tuning,xtest,ytest, predicttype=NULL,probability=TRUE){
  
  grid <- expand.grid(tuning)
  
  perf <- numeric()
  for (i in 1:nrow(grid)){
  Call <- c(as.list(call), grid[i,])
  model <-  eval(as.call(Call))
 
  predictions <- predict(model,xtest,type=predicttype, probability=probability)    
  
  if (class(model)[2] == "svm") predictions <- attr(predictions,"probabilities")[,2]
  
     
  if (is.matrix(predictions)) if (ncol(predictions) == 2 ) predictions <- predictions[,2]
  perf[i] <- performance(prediction(predictions,ytest),"auc")@y.values[[1]]
  }
  perf <- data.frame(grid, auc=perf)
  perf[which.max(perf$auc),]
}


################################################################################################################################


#currently limited 
#-to binary problems only, in the future it will work for more classes too, with a loop.
#-repeated cross validation (e.g., 5 times 2 fold, 1 times 2 fold)

.partition <- function(y,p=0.5,times=1) {
  
  #STEP 1: split up 0 and 1
  class1_ind <- which(y==as.integer(levels(y)[1]))
  class2_ind <- which(y==as.integer(levels(y)[2]))
  
  l <- list()
  for (i in 1:times){
  
  #STEP 2: take subsamples for both 0 and 1
  class1_ind_train <- sample(class1_ind, floor(p*table(y)[1]),replace=FALSE)
  class2_ind_train <- sample(class2_ind, floor(p*table(y)[2]),replace=FALSE)

  class1_ind_test <- class1_ind[!class1_ind %in% class1_ind_train]
  class2_ind_test <- class2_ind[!class2_ind %in% class2_ind_train]

  #STEP 3: combine 0 and 1 for both train and test
  
  l[[i]] <- list(train=c(class1_ind_train,class2_ind_train),test=c(class1_ind_test,class2_ind_test))
  }
  l
}



################################################################################################################################




#maak equal frequency bins
#average probabilities en proportion labels
#op die binned dataset een regressie (bv., glm, random Forest) voor proportions
#gebruik dat model via predict om nieuwe data (i.e., test data predictions) te calibreren
#rule base cre?ren die ervoor zorgt dat de calibrated probabilities isotonic zijn



# x : numeric vector of predicted probabilities from a classifier on a validation set
# y:  factor of observed labels on a validation set
.calibrate <- function(x,y) {
  
  trainIND <- .partition(y,p=0.8)[[1]]$train
  
  xTRAIN <- x[trainIND]
  yTRAIN <- y[trainIND]
  
  xVALIDATE <- x[-trainIND]
  yVALIDATE <- y[-trainIND]
  
    #DETERMINE OPTIMAL NUMBER OF BREAKS 
    AUC <- data.frame(matrix(ncol=2))
    i <- 0
    for (nbreaks in 2: if (length(xTRAIN) > 1000) 500 else length(xTRAIN))  {
       
      #create equal frequency bins
      x_bin <- cut(xTRAIN,breaks=nbreaks ,labels=FALSE)
      x_mean <- data.frame(aggregate(xTRAIN,by=list(x_bin),mean)$x)
      names(x_mean) <- "x_mean"
      y_prop <- aggregate(as.integer(as.character(yTRAIN)),by=list(x_bin),mean)$x
      y_prop <- cummax(y_prop)
     
      if (length(unique(y_prop)) > 5) {
    
      
        names(x_mean) <- "x"
        rf <- randomForest(y=y_prop,x=x_mean)
    
        xVALIDATE <- data.frame(xVALIDATE)
        names(xVALIDATE) <- "x"
        predrfCAL <- predict(object=rf, newdata=xVALIDATE)
    
        i <- i + 1
        AUC[i,1:2] <- c(performance(prediction(predrfCAL,yVALIDATE), measure="auc")@y.values[[1]],nbreaks)
      }
      names(AUC) <- c("AUC","nbreaks")
    }  
  
  nbreaks <- AUC$nbreaks[which.max(AUC$AUC)]
  
  if (length(nbreaks) != 0) {
  
      #ESTIMATE FINAL MODEL USING OPTIMAL NUMBER OF BREAKS
      #create equal frequency bins
      x_bin <- cut(x,breaks=nbreaks ,labels=FALSE)
      x_mean <- data.frame(aggregate(x,by=list(x_bin),mean)$x)
      names(x_mean) <- "x_mean"
      y_prop <- aggregate(as.integer(as.character(y)),by=list(x_bin),mean)$x
      y_prop <- cummax(y_prop)
    
      result <- list(probmapper=randomForest(y=y_prop,x=x_mean),performance=AUC,nbrbreaks=nbreaks)
  } else {
      result <- list(rawprobs=x)
  }
  
  class(result) <- "calibrate"
  result
} 





#object: calibrate object
#newdata: test data

.predict.calibrate <- function(object,newdata) {
  if (length(object)!=1) {
    newdata <- data.frame(newdata)
    names(newdata) <- "x_mean"
    pr <- predict(object=object[[1]], newdata=newdata, type="response")
  } else {
    pr <- newdata
  }
  pr  
#   df <- data.frame(rn=as.integer(rownames(newdata)),newdata,pr)
#   df <- df[order(df[,2]),]
#   cummax(df$pr)[order(df$rn)]
#  }

}



#binary, unbinary, bit and maxBit are from the compositions package. We copied them here because there were collisions with their S3methods

# A function to format a number in binary digits
# x the number
# mb the number of binary digits
.binary <- function(x,mb=max(.maxBit(x,g)),g=2) {
  if( is.character(x) ) x<- .unbinary(x)
  if( g==2 )
    do.call(paste,c(sep="",lapply(mb:0,function(i) ifelse(.bit(x,i,g=2),"1","0"))))
  else{
    .toDigit <- function(x) c(0:9,LETTERS)[x+1]
    do.call(paste,c(sep="",lapply(mb:0,function(i) .toDigit(.bit(x,i,g=g)))))
  }
}
# Converts a binary character string to a number
.unbinary <- function(x,g=2) {
  if( is.numeric(x) )
    return(x)
  nc =nchar(x)
  D = max(max(nchar(x)),2)
  .asDigit <- (if(g==2) function(x) as.logical((match(x,c("0","1","F","T"))-1)%%2) else function(x) as.integer((match(x,as.character(c(0:9,LETTERS,10:19,letters)[c(1:g,37:(36+g))]))-1)%%g))
  c(sapply(1:D,function(i) ifelse(i<=nc,.asDigit(substring(x,i,i))*g^(nc-i),rep(0,length(x))))%*%rep(1,D)) 
}

# a function to extract a bit from a binary number
# either given as number or as character string
# x the number or string (may be vectors)
# b the bit to be extracted (may be a vector)
.bit <- function(x,b,g=2) UseMethod(".bit")                       
.bit.numeric   <- function(x,b=0:.maxBit(x,g),g=2)  {
  erg <- sapply(b,function(b) (x%/% (g^b) %% g ))
  structure((if(g==2) as.logical else as.integer)(erg),dim=dim(erg))
}
.bit.character <- function(x,b=0:.maxBit(x,g),g=2)  {
  nc = nchar(x)
  .asDigit <- (if(g==2) function(x) as.logical((match(x,c("0","1","F","T"))-1)%%2) else function(x) as.integer((match(x,as.character(c(0:9,LETTERS,10:19,letters)[c(1:g,37:(36+g))]))-1)%%g))
  erg <- sapply(b,function(b) ifelse(b<nc,substring(x,nc-b,nc-b),"0"))
  structure(.asDigit(erg),dim=dim(erg))
}


.maxBit <- function(x,g=2) UseMethod(".maxBit")
.maxBit.numeric <- function(x,g=2) ceiling(log(x+1,g))-1
.maxBit.character <- function(x,g=2) max(nchar(x))-1
