#' Binary classification with Kernel Factory
#'
#' \code{kernelFactory} implements an ensemble method for kernel machines (Ballings and Van den Poel, 2013).
#' 
#' @param x A data frame of predictors (numeric, integer or factor). Categorical variables need to be factors. Indicator values should not be too imbalanced because this might produce constants in the subsetting process.
#' @param y A factor containing the response vector. Only \{0,1\} is allowed.
#' @param cp The number of column partitions.
#' @param rp The number of row partitions.
#' @param method Can be one of the following: POLynomial kernel function (\code{pol}), LINear kernel function (\code{lin}), Radial Basis kernel Function \code{rbf}), random choice (random={pol, lin, rbf}) (\code{random}), burn- in choice of best function (burn={pol, lin, rbf }) (\code{burn}). Use \code{random} or \code{burn} if you don't know in advance which kernel function is best.
#' @param ntree Number of trees in the Random Forest base classifiers.
#' @param filter either NULL (deactivate) or a percentage denoting the minimum class size of dummy predictors. This parameter is used to remove near constants. For example if nrow(xTRAIN)=100, and filter=0.01 then all dummy predictors with any class size equal to 1 will be removed. Set this higher (e.g., 0.05 or 0.10) in case of errors.
#' @param popSize Population size of the genetic algorithm.
#' @param iters  Number of generations of the genetic algorithm.
#' @param mutationChance Mutationchance of the genetic algorithm.
#' @param elitism  Elitism parameter of the genetic algorithm.
#' @param oversample  Oversample the smallest class. This helps avoid problems related to the subsetting procedure (e.g., if rp is too high).
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
#' #predictedresponse <- predict(kFmodel, newdata=Credit[-train.ind,names(Credit)!= "Response"])
#' @references Ballings, M. and Van den Poel, D. (2013), Kernel Factory: An Ensemble of Kernel Machines. Expert Systems With Applications, 40(8), 2904-2913.
#' @seealso \code{\link{predict.kernelFactory}}
#' @return An object of class \code{kernelFactory}, which is a list with the following elements:
#'   \item{trn}{Training data set.}
#'   \item{trnlst}{List of training partitions.}
#'   \item{rbfstre}{List of used kernel functions.}
#'   \item{rbfmtrX}{List of augmented kernel matrices.}
#'   \item{rsltsKF}{List of models.}
#'   \item{cpr}{Number of column partitions.}
#'   \item{rpr}{Number of row partitions.}
#'   \item{cntr}{Number of partitions.}
#'   \item{wghts}{Weights of the ensemble members.}
#'   \item{nmDtrn}{Vector indicating the numeric (and integer) features.}
#'   \item{rngs}{Ranges of numeric predictors.}
#'   \item{constants}{To exclude from newdata.}
#' @author Authors: Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
#' @keywords classification
kernelFactory <-
function(x=NULL, 
	   y=NULL, 
	   cp=1, 
	   rp=round(log(nrow(x),10)), 
	   method="burn" ,
     ntree=500,
     filter= 0.01,
     popSize = rp*cp*7,
     iters = 80, 
     mutationChance = 1/ (rp*cp),
     elitism= max(1, round((rp*cp)*0.05)),
     oversample=TRUE){
  
  
  #ERROR HANDLING
  if (!is.data.frame(x)) stop("x must be a data frame")
  
  
	if (is.null(x) || is.null(y)) {
		stop("x or y cannot be NULL.")
	}else if (any(is.na(x)) || any(is.na(y))){
		stop("NAs not permitted.")
	}

	if (!is.factor(y)) stop("y must be a factor. Support for regression will be added later.")
	
	
	if (any(table(y) == 0)) stop("Cannot have empty classes in y.")
	
	if (length(unique(y)) != 2) stop("Must have 2 classes. Support for more classes will be added later.")

	if (length(y) != nrow(x)) stop("x and dependent variable have to be of equal length.")


  #OVERSAMPLING
  tab <- table(y)
  if (!all(tab==tab[1])){
      if (oversample) {

        #oversample instances from the smallest class
        whichmin <- which(y==as.integer(names(which.min(tab))))
        indmin <- sample(whichmin,max(tab),replace=TRUE)
        indmin <- c(whichmin,indmin)[1:max(tab)]
        #take all the instances of the dominant class
        indmax <- which(y==as.integer(names(which.max(tab))))
        x <- x[c(indmin,indmax),]
        y <- y[c(indmin,indmax)]    
      }
  }
  
  #STEP0 SEPERATE DATA INTO TRAINING AND VALIDATION

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
  
  train.ind <- .partition(y,0.75)[[1]]$train
	
  xtrain <- x[train.ind,]
  ytrain <- y[train.ind]
  
  xtest <-  x[-train.ind,]
  ytest <-  y[-train.ind]
  
	constants <- sapply(xtrain,function(x){all(as.numeric(x[1])==as.numeric(x))})
  if (!is.null(filter)) constants <- sapply(xtrain,function(x) (length(unique(x))==2 && any(table(x) <= round(nrow(xtrain)*filter))) || all(as.numeric(x[1])==as.numeric(x))   )
  xtrain <- xtrain[,!constants]
  xtest <- xtest[,!constants]


  #STEP1 ### SCALING 
	
  #select only numerics and integers
	numIDtrain <- sapply(xtrain,is.numeric)
  
	numericcolumnsTRAIN <- data.frame(xtrain[,numIDtrain ])
	numericcolumnsVAL <- data.frame(xtest[,numIDtrain ])
  colnames(numericcolumnsTRAIN) <- colnames(numericcolumnsVAL) <- colnames(xtrain)[numIDtrain]

  #scale
  ranges <- sapply(numericcolumnsTRAIN ,function(x) max(x)-min(x))
  numericcolumnsTRAIN <- data.frame(base::scale(numericcolumnsTRAIN,center=FALSE,scale=ranges))
  numericcolumnsVAL   <- data.frame(base::scale(numericcolumnsVAL,center=FALSE,scale=ranges))
  
  
  #add factor variables and response
  facID <- sapply(xtrain,is.factor)
  if (sum(facID) >= 1){
  
  
    train <- data.frame(numericcolumnsTRAIN,xtrain[,facID],Y=ytrain)
    test <- data.frame(numericcolumnsVAL,xtest[,facID],Y=ytest)
    
    # in case there is only 1 factor: set colnames
    colnames(train) <- colnames(test) <- c(colnames(numericcolumnsTRAIN),names(facID)[facID],"Y")
  } else {
      	
		train <- data.frame(numericcolumnsTRAIN, Y=ytrain)
		test  <- data.frame(numericcolumnsVAL, Y=ytest)
  }
  
  rm(xtrain,xtest,ytrain,ytest,numericcolumnsTRAIN,numericcolumnsVAL)
  

#STEP2 ### RANDOMLY ORDER ROWS AND COLUMNS
	train <- data.frame(train[order(runif(nrow(train))),order(runif(ncol(train)))])
	train <- data.frame(train[,names(train)!= "Y"], Y=train$Y)
	
	test <- test[,colnames(train)]

  

#STEP3 ### CREATE PARTITIONS 

  #STEP 3.1 ### PARTITIONS FOR TRAINING SET	
	trainlist <- list()
	colparts <- cp
	rowparts <- rp

	startcol <- 1
	counter <- 0
	
  #HANDLE UNEVEN NUMBER OF VARIABLES	
	if (ncol(train[,names(train)!= "Y"])%%colparts > 0) {
	numbercols <- ncol(train[,names(train)!= "Y"])-(ncol(train[,names(train)!= "Y"])%%colparts )
	} else {
	numbercols <- ncol(train[,names(train)!= "Y"]) }

	 #HANDLE UNEVEN NUMBER OF ROWS
	if (nrow(train)%%rowparts > 0) {
	numberrows <- nrow(train)-(nrow(train)%%rowparts )
	} else {
	numberrows <- nrow(train)}


	for (i in 1:colparts ){

		 
		endcol <-(i/colparts )*numbercols 
		if (i==colparts ) endcol <- ncol(train[,names(train)!= "Y"])
		
		startrow <- 1
		for (j in 1:rowparts) {
			counter = counter + 1
			endrow <-(j/rowparts )*numberrows
			if (j==rowparts) endrow <- nrow(train)


			trainlist[[counter]] <- train[startrow:endrow ,c(startcol:endcol,which(colnames(train)=="Y"))]
	
			
			if (any(table(trainlist[[counter]]$Y) == 0)) stop("Cannot have empty classes in y. Make sure number of rp is not too high.")
			
			if (length(unique(trainlist[[counter]]$Y)) != 2) stop("Must have 2 classes. Make sure number of rp is not too high.")

			startrow <- endrow + 1
		}

		startcol <- endcol + 1
	}



	#STEP 3.2 ### PARTITIONS FOR TEST SET
	testlist <- list()
	colparts <- colparts
	rowparts <- rowparts 

	startcol <- 1
	counter <- 0
	#HANDLE UNEVEN NUMBER OF VARIABLES
	if (ncol(test[,names(test)!= "Y"])%%colparts > 0) {
	numbercols <- ncol(test[,names(test)!= "Y"])-(ncol(test[,names(test)!= "Y"])%%colparts )
	} else {
	numbercols <- ncol(test[,names(test)!= "Y"]) }



	for (i in 1:colparts ){

		 
			endcol <-(i/colparts )*numbercols 
			if (i==colparts ) endcol <- ncol(test[,names(test)!= "Y"])
		
			startrow <- 1
			for (j in 1:rowparts) {
				counter = counter + 1
			
				testlist[[counter]] <- test[,c(startcol:endcol,which(colnames(test)=="Y"))]
				#Checks: no empty classes in dependent
				if (any(table(testlist[[counter]]$Y) == 0)) stop("Cannot have empty classes in y. Make sure number of rp is not too high.")
				#Can only have 2 classes (not less, not more)
				if (length(unique(testlist[[counter]]$Y)) != 2) stop("Must have 2 classes. Make sure number of rp is not too high.")

	
			
			}

			startcol <- endcol + 1
	}





#STEP4 ### BUILD MODELS 	
	rbfstore <- list()
	rbfmatrX <- data.frame()
	resultsKF <- data.frame()


	#STEP 4.0 APPLY METHOD	
	
	
	if (as.character(substitute(method))=="rbf") {
		for (ii in 1:counter) {
			rbfstore[[ii ]]<- rbfdot(sigma = 1)
		}
	} else if (as.character(substitute(method))=="pol") {
		for (ii in 1:counter) {
			rbfstore[[ii ]]<- polydot(degree=2,scale=1)
		}
	} else if (as.character(substitute(method))=="lin") {
		for (ii in 1:counter) {
			rbfstore[[ii ]]<- vanilladot()
		}
	} else if (as.character(substitute(method))=="random") {
		for (ii in 1:counter) {
			randomnumber <- runif(1,min=0, max=1)
			if(randomnumber <= 0.33) {
				rbfstore[[ii ]]<- rbfdot(sigma = 1)
			} else if (randomnumber > 0.33 & randomnumber <= 0.66 ) {
				rbfstore[[ii ]]<- polydot(degree=2,scale=1)
			} else if (randomnumber > 0.66 ) {
				rbfstore[[ii ]]<- vanilladot()
			}
		}
	} else if (as.character(substitute(method))=="burn") {
		rbfstore <- list(rbfdot(sigma = 1),polydot(degree=2,scale=1),vanilladot())
		burnperf <- list()
		for (ii in 1:3) {
			
			numID <- sapply(trainlist[[1]],is.numeric)
			numericcolumns <- trainlist[[1]][,numID]
			
			if (is.null(trainlist[[1]][,sapply(trainlist[[1]],is.numeric)]) == FALSE) {
								
											
				rbfdt<-as.matrix(numericcolumns)
				rbfmatr<-kernelMatrix(rbfstore[[ii]], rbfdt)
				rbfmatr <-  data.frame(rbfmatr)


				rbfmatrX <- data.frame(rbfmatr, trainlist[[1]][,names(trainlist[[1]])!= "Y"])
				rbfmatrY <- trainlist[[1]][,which(colnames(trainlist[[1]])=="Y")]
                   
				rm(numID)
				rm(numericcolumns)
        								
            
		
				
				resultsKF  <-  	 randomForest(x=rbfmatrX,y=as.factor(rbfmatrY),  ntree=ntree, proximity=FALSE, importance=FALSE, na.action=na.omit)
				
				
				#Extract split variables for all trees: these are observation numbers from the trainingdata: split observations
				trainobs <- trainlist[[1]][,sapply(trainlist[[1]],is.numeric)]
								
																
				#compute dot product per split observation from training data with all validation observations: This results in 1 column per split observation
           
								
				resultsKFScored <- data.frame(data.frame(kernelMatrix(rbfstore[[ii]], as.matrix(testlist[[1]][,sapply(testlist[[1]],is.numeric)]),as.matrix(trainobs))),
													data.frame(testlist[[1]][,names(testlist[[1]])!= "Y"]))
                 
          			
				colnames(resultsKFScored) <- colnames(rbfmatrX)
					
				predicted <- predict(resultsKF,resultsKFScored,type="prob")[,2]
			
				burnperf[[ii]] <- AUC::auc(roc(predicted,testlist[[1]]$Y))
                          
			} else {
				
		 		resultsKF <-  randomForest(x=trainlist[[1]][,names(trainlist[[1]])!= "Y"],y=as.factor(trainlist[[1]]$Y),  ntree=ntree, proximity=FALSE, importance=FALSE, na.action=na.omit )
				
				
				predicted <- predict(resultsKF,testlist[[1]],type="prob")[,2]

				burnperf[[ii]] <- AUC::auc(roc(predicted,testlist[[1]]$Y))
                          
  			}

		}
		
		bestkernelfunction <- rbfstore[[which.max(burnperf)]]
		rbfstore <- list()
		for (ii in 1:counter) {
			rbfstore[[ii]] <- bestkernelfunction
		}
				
	}


	rbfmatrX <- list()
	resultsKF <- list()


	for (i in 1:counter) {

			#STEP 4.1 ### MODELING
								
			numID <- sapply(trainlist[[i]],is.numeric)
			numericcolumns <- trainlist[[i]][,numID]
				
			
			if (is.null(trainlist[[i]][,sapply(trainlist[[i]],is.numeric)]) == FALSE) {
								
				
				#STEP 4.1 ### COMPUTE KERNEL MATRIX
				rbfdt<-as.matrix(numericcolumns)
				rbfmatr<-kernelMatrix(rbfstore[[i]], rbfdt)
				rbfmatr <-  data.frame(rbfmatr)


				rbfmatrX[[i]] <- data.frame(rbfmatr, trainlist[[i]][,names(trainlist[[i]])!= "Y"])
				rbfmatrY <- trainlist[[i]][,which(colnames(trainlist[[i]])=="Y")]
                   
				rm(numID)
				rm(numericcolumns)
        								


				
        #STEP 4.3 BUILDING KERNEL MODELS    
		
				
				resultsKF[[i]] <-   randomForest(x=rbfmatrX[[i]],y=as.factor(rbfmatrY),  ntree=ntree, proximity=FALSE, importance=FALSE, na.action=na.omit)
				

			} else {
				
		 		resultsKF[[i]]  <-   randomForest(x=trainlist[[i]][,names(trainlist[[i]])!= "Y"],y=as.factor(trainlist[[i]]$Y),  ntree=ntree, proximity=FALSE, importance=FALSE, na.action=na.omit )
				

			}

			
	}



#STEP5 ### OPTIMIZE WEIGHTS

  #STEP 5.1 PREPARING 
	
	predicted <- list()
	resultsKFScored <- list()




	for (i in 1:counter) {                
	
			if (is.null(trainlist[[i]][,sapply(trainlist[[i]],is.numeric)]) == FALSE) {
	
				#Extract split variables for all trees: these are observation numbers from the trainingdata: split observations
				trainobs <- trainlist[[i]][,sapply(trainlist[[i]],is.numeric)]
								
																
				
        #compute dot product per split observation from training data with all validation observations: This results in 1 column per split observation   
								
				resultsKFScored[[i]] <- data.frame(data.frame(kernelMatrix(rbfstore[[i]], as.matrix(testlist[[i]][,sapply(testlist[[i]],is.numeric)]),as.matrix(trainobs))),
													data.frame(testlist[[i]][,names(testlist[[i]])!= "Y"]))
                 
          			
				colnames(resultsKFScored[[i]]) <- colnames(rbfmatrX[[i]])
	
				#STEP 5.2 ACTUAL PREDICTION
				predicted[[i]] <- predict(resultsKF[[i]],resultsKFScored[[i]],type="prob")[,2]
				


	
			} else {
				predicted[[i]] <- predict(resultsKF[[i]],testlist[[i]],type="prob")[,2]
				

			
			}
	}
	



	#STEP 5.3 COLLECTING PREDICTIONS
	final <- data.frame(matrix(nrow=nrow(test),ncol=(counter)))
	for (i in 1:counter) {
	final[,i] <- as.numeric(predicted[[i]])
	}


	#STEP 5.4 Genetic algorithm
	

	evaluate <- function(string=c()) {
 
     -AUC::auc(roc(as.numeric(rowSums(t(as.numeric(as.numeric(string)/sum(as.numeric(string))) * t(final)))),testlist[[i]]$Y ))
   	   
	}


	rbga.results = rbga(rep(0,counter), rep(1,counter), 
	                    suggestions=t( as.data.frame( c(rep((1/counter),counter)))), 
                      popSize=popSize, 
                      iters=iters,  
                      mutationChance=mutationChance, 
                      elitism=elitism,
                      evalFunc=evaluate)

	

	weights <- rbga.results$population[which.min(rbga.results$evaluations),]

	weights <- as.numeric(weights)/sum(as.numeric(weights))


result <- list(trn=train, 						
		   trnlst=trainlist,  					
		   rbfstre=rbfstore, 					
		   rbfmtrX=rbfmatrX, 					
		   rsltsKF=resultsKF, 					
		   cpr=cp,  						
		   rpr=rp, 	 						
		   cntr=counter,  					
		   wghts=weights,  					
		   nmDtrn=numIDtrain ,  					
		   rngs=ranges,
       constants=constants  ) 

class(result) <- "kernelFactory"					

return(result)
}
