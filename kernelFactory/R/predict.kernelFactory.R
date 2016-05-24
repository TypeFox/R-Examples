#' Predict method for kernelFactory objects
#'
#' Prediction of new data using kernelFactory.
#' 
#' @param object An object of class \code{kernelFactory}, as created by the function \code{kernelFactory}
#' @param newdata A data frame with the same predictors as in the training data.
#' @param predict.all TRUE or FALSE. If TRUE and rp and cp are 1 then the individual predictions of the random forest are returned. If TRUE and any of rp and cp or bigger than 1 then the predictions of all the members are returned.
#' @param ... Not used currently.
#' 
#' @examples
#' #Credit Approval data available at UCI Machine Learning Repository
#' data(Credit)
#' #take subset (for the purpose of a quick example) and train and test
#' Credit <- Credit[1:100,]
#' train.ind <- sample(nrow(Credit),round(0.5*nrow(Credit)))
#' 
#' #Train Kernel Factory on training data
#' kFmodel <- kernelFactory(x=Credit[train.ind,names(Credit)!= "Response"], 
#'           y=Credit[train.ind,"Response"], method=random)
#'  
#' #Deploy Kernel Factory to predict response for test data
#' predictedresponse <- predict(kFmodel, newdata=Credit[-train.ind,names(Credit)!= "Response"])
#' @references Ballings, M. and Van den Poel, D. (2013), Kernel Factory: An Ensemble of Kernel Machines. Expert Systems With Applications, 40(8), 2904-2913.
#' @seealso \code{\link{kernelFactory}}
#' @return A vector containing the response probabilities.
#' @author Authors: Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
#' @method predict kernelFactory
#' @keywords classification
predict.kernelFactory <-
function(object, 
   newdata=NULL, predict.all=FALSE, ...  )
{
  
  newdata <- newdata[,!object$constants]
  
  #ERROR HANDLING
  if (!inherits(object, "kernelFactory")) stop("object not of class kernelFactory")
  
  if (is.null(object)) stop("You must provide the name of the model through object.")
  if (is.null(newdata)) {
  stop("Data cannot be NULL.")
  }else if (any(is.na(newdata))){
  stop("NAs not permitted.")
  }
  if(all(colnames(object$trn[,names(object$trn)!= "Y"]) %in% colnames(newdata))==FALSE ) stop("Column names have to be equal in data and new data.")


 
#STEP1 ### SCALING 

  numericcolumns <- newdata[,object$nmDtrn]
  numericcolumns <- data.frame(base::scale(numericcolumns,center=FALSE,scale=object$rngs))
  
  
  #add factor variables
  facID <- sapply(newdata,is.factor)
  if (sum(facID) >= 1){
    
    newdata <- data.frame(numericcolumns,newdata[,facID])
   
    # in case there is only 1 factor: set colnames
    colnames(newdata) <-  c(colnames(numericcolumns),names(facID)[facID])
  } else {
        
		newdata <- data.frame(numericcolumns)
		
  }
  
  
  rm(numericcolumns)


## ORDER VARIABLES ACCORDING TO TRAINING SET
newdata <- newdata[,colnames(object$trn[,names(object$trn)!= "Y"])]


#STEP3 ### CREATE PARTITIONS 


testlist <- list()
colparts <- object$cpr
rowparts <- object$rpr

startcol <- 1
counter <- 0
#HANDLE UNEVEN NUMBER OF VARIABLES

if (ncol(newdata)%%colparts > 0) {
numbercols <- ncol(newdata)-(ncol(newdata)%%colparts )
} else {
numbercols <-  ncol(newdata)  }



for (i in 1:colparts ){

  
endcol <-(i/colparts )*numbercols 
if (i==colparts ) endcol <- ncol(newdata)

startrow <- 1
for (j in 1:rowparts) {
counter = counter + 1

testlist[[counter]] <- newdata[,c(startcol:endcol)]



}

startcol <- endcol + 1
}

#STEP4 ### SCORING


#STEP 5.1 PREPARING 
predicted <- list()
resultsKIRFScored <- list()
for (i in 1:object$cntr) {                


if (is.null(testlist[[i]][,sapply(testlist[[i]],is.numeric)]) == FALSE) {


#Extract split variables for all trees: these are observation numbers from the trainingdata: split observations
trainobs <- object$trnlst[[i]][,sapply(object$trnlst[[i]],is.numeric)]



           
#compute dot product per split observation from training data with all validation observations: This results in 1 column per split observation
resultsKIRFScored[[i]] <- data.frame(data.frame(kernelMatrix(object$rbfstre[[i]], as.matrix(testlist[[i]][,sapply(testlist[[i]],is.numeric)]),as.matrix(trainobs))),
data.frame(testlist[[i]]))
                 
          
colnames(resultsKIRFScored[[i]]) <- colnames(object$rbfmtrX[[i]])

#STEP 5.2 ACTUAL PREDICTION
    if (predict.all && object$cpr==1 && object$rpr==1 ) { # predict.all
    predicted[[i]]  <- data.frame(sapply(data.frame(predict(object$rsltsKF[[i]], 
                                                        resultsKIRFScored[[i]],type="prob", 
                                                        predict.all=TRUE)$individual,
                                                stringsAsFactors=FALSE),
                                     as.integer,
                                     simplify=FALSE)
                              )
  } else { #combined prediction
    predicted[[i]] <- predict(object$rsltsKF[[i]],resultsKIRFScored[[i]],type="prob")[,2]
  }


} else {
  
  if (predict.all && object$cpr==1 && object$rpr==1 ) { #predict.all
    predicted[[i]]  <- data.frame(sapply(data.frame(predict(object$rsltsKF[[i]], 
                                                        testlist[[i]],type="prob", 
                                                        predict.all=TRUE)$individual,
                                                stringsAsFactors=FALSE),
                                     as.integer,
                                     simplify=FALSE)
                              )
  } else { #combined prediction
    predicted[[i]] <- predict(object$rsltsKF[[i]],testlist[[i]],type="prob")[,2]
  }

}
}

if (predict.all && object$cpr==1 && object$rpr==1  ) {
  return(predicted[[1]])
} else {
  
    #STEP 5.3 COLLECTING PREDICTIONS
    final <- data.frame(matrix(nrow=nrow(newdata),ncol=(object$cntr)))
    for (i in 1:object$cntr) {
    final[,i] <- predicted[[i]]
    }

    
    #STEP 5.4 APPLYING WEIGHTS  
    res <- t(as.numeric(object$wghts) * t(final))

    if (predict.all) {
      return(data.frame(final))
    } else {
      return(rowSums(res))
    }
}
  
}
