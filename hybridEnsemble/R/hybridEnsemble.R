#' Binary classification with Hybrid Ensemble
#'
#' \code{hybridEnsemble} builds an ensemble consisting of seven different sub-ensembles: Bagged Logistic Regressions, Random Forest, Stochastic AdaBoost, Kernel Factory, Bagged Neural Networks, Bagged Support Vector Machines, and Rotation Forest.
#' 
#' @param x A data frame of predictors. Categorical variables need to be transformed to binary (dummy) factors.
#' @param y A factor of observed class labels (responses) with the only allowed values \{0,1\}.,
#' @param combine Additional methods for combining the sub-ensembles. The simple mean, authority-based weighting and the single best are automatically provided since they are very efficient.  Possible additional methods: Genetic Algorithm: "rbga", Differential Evolutionary Algorithm: "DEopt", Generalized Simulated Annealing: "GenSA", Memetic Algorithm with Local Search Chains: "malschains", Particle Swarm Optimization: "psoptim", Self-Organising Migrating Algorithm: "soma", Tabu Search Algorithm: "tabu", Non-negative binomial likelihood: "NNloglik", Goldfarb-Idnani Non-negative least squares: "GINNLS", Lawson-Hanson Non-negative least squares: "LHNNLS".
#' @param eval.measure Evaluation measure for the following combination methods: authority-based method, single best, "rbga", "DEopt", "GenSA", "malschains", "psoptim", "soma", "tabu". Default is the area under the receiver operator characteristic curve 'auc'. The area under the sensitivity curve ('sens') and the area under the specificity curve ('spec') are also supported.
#' @param verbose TRUE or FALSE. Should information be printed to the screen while estimating the Hybrid Ensemble.
#' @param oversample TRUE or FALSE. Should oversampling be used? Setting oversample to TRUE helps avoid computational problems related to the subsetting process.
#' @param filter either NULL (deactivate) or a percentage denoting the minimum class size of dummy predictors. This parameter is used to remove near constants. For example if nrow(xTRAIN)=100, and filter=0.01 then all dummy predictors with any class size equal to 1 will be removed. Set this higher (e.g., 0.05 or 0.10) in case of errors.
#' @param LR.size Logistic Regression parameter. Ensemble size of the bagged logistic regression sub-ensemble.
#' @param RF.ntree Random Forest parameter. Number of trees to grow.
#' @param AB.iter Stochastic AdaBoost parameter. Number of boosting iterations to perform.
#' @param AB.maxdepth Stochastic AdaBoost parameter. The maximum depth of any node of the final tree, with the root node counted as depth 0.
#' @param KF.cp Kernel Factory parameter. The number of column partitions.
#' @param KF.rp Kernel Factory parameter. The number of row partitions.
#' @param KF.ntree Kernel Factory parameter. Number of trees to grow.
#' @param NN.rang Neural Network parameter. Initial random weights on [-rang, rang].
#' @param NN.maxit Neural Network parameter. Maximum number of iterations. 
#' @param NN.size Neural Network parameter. Number of units in the single hidden layer. Can be mutiple values that need to be optimized.
#' @param NN.decay Neural Network parameter. Weight decay. Can be mutiple values that need to be optimized.
#' @param NN.ens.size Neural Network parameter. Ensemble size of the neural network sub-ensemble.
#' @param SV.gamma Support Vector Machines parameter. Width of the Guassian for radial basis and sigmoid kernel. Can be mutiple values that need to be optimized.
#' @param SV.cost Support Vector Machines parameter. Penalty (soft margin constant). Can be mutiple values that need to be optimized.
#' @param SV.degree Support Vector Machines parameter. Degree of the polynomial kernel. Can be mutiple values that need to be optimized.
#' @param SV.kernel Support Vector Machines parameter. Kernels to try. Can be one or more of: 'radial','sigmoid','linear','polynomial'. Can be mutiple values that need to be optimized.
#' @param SV.size Support Vector Machines parameter. Ensemble size of the SVM sub-ensemble.
#' @param RoF.L Rotation Forest parameter. Number of trees to grow.
#' @param KNN.K K-Nearest Neigbhors parameter. Number of nearest neighbors to try. For example c(10,20,30). The optimal K will be selected. If larger than nrow(xTRAIN) the maximum K will be reset to 50\% of nrow(xTRAIN). Can be mutiple values that need to be optimized.
#' @param KNN.size K-Nearest Neigbhors parameter. Ensemble size of the K-nearest neighbor sub-ensemble.
#' @param rbga.popSize Genetic Algorithm parameter. Population size.
#' @param rbga.iters Genetic Algorithm parameter.  Number of iterations.
#' @param rbga.mutationChance Genetic Algorithm parameter. The chance that a gene in the chromosome mutates.
#' @param rbga.elitism Genetic Algorithm parameter. Number of chromosomes that are kept into the next generation.
#' @param DEopt.nP Differential Evolutionary Algorithm parameter. Population size.
#' @param DEopt.nG Differential Evolutionary Algorithm parameter. Number of generations.
#' @param DEopt.F Differential Evolutionary Algorithm parameter. Step size.
#' @param DEopt.CR Differential Evolutionary Algorithm parameter. Probability of crossover.
#' @param GenSA.maxit Generalized Simulated Annealing. Maximum number of iterations.
#' @param GenSA.temperature Generalized Simulated Annealing. Initial value for temperature.
#' @param GenSA.visiting.param Generalized Simulated Annealing. Parameter for visiting distribution.
#' @param GenSA.acceptance.param Generalized Simulated Annealing. Parameter for acceptance distribution.
#' @param GenSA.max.call Generalized Simulated Annealing. Maximum number of calls of the objective function.
#' @param malschains.popsize Memetic Algorithm with Local Search Chains parameter. Population size. 
#' @param malschains.ls Memetic Algorithm with Local Search Chains parameter. Local search method.
#' @param malschains.istep Memetic Algorithm with Local Search Chains parameter. Number of iterations of the local search.
#' @param malschains.effort Memetic Algorithm with Local Search Chains parameter. Value between 0 and 1. The ratio between the number of evaluations for the local search and for the evolutionary algorithm. A higher effort means more evaluations for the evolutionary algorithm. 
#' @param malschains.alpha Memetic Algorithm with Local Search Chains parameter. Crossover BLX-alpha. Lower values (<0.3) reduce diversity and a higher value increases diversity.
#' @param malschains.threshold Memetic Algorithm with Local Search Chains parameter. Threshold that defines how much improvement in the local search is considered to be no improvement.
#' @param malschains.maxEvals Memetic Algorithm with Local Search Chains parameter. Maximum number of evaluations.
#' @param psoptim.maxit Particle Swarm Optimization parameter. Maximum number of iterations.
#' @param psoptim.maxf Particle Swarm Optimization parameter. Maximum number of function evaluations.
#' @param psoptim.abstol Particle Swarm Optimization parameter. Absolute convergence tolerance.
#' @param psoptim.reltol Particle Swarm Optimization parameter. Tolerance for restarting.
#' @param psoptim.s Particle Swarm Optimization parameter. Swarm size.
#' @param psoptim.k Particle Swarm Optimization parameter. Exponent for calculating number of informants.
#' @param psoptim.p Particle Swarm Optimization parameter. Average percentage of informants for each particle.
#' @param psoptim.w Particle Swarm Optimization parameter. Exploitation constant.
#' @param psoptim.c.p Particle Swarm Optimization parameter. Local exploration constant.
#' @param psoptim.c.g Particle Swarm Optimization parameter. Global exploration constant.
#' @param soma.pathLength Self-Organising Migrating Algorithm parameter. Distance (towards the leader) that individuals may migrate.
#' @param soma.stepLength Self-Organising Migrating Algorithm parameter. Granularity at which potential steps are evaluated.
#' @param soma.perturbationChance Self-Organising Migrating Algorithm parameter. Probability that individual parameters are changed on any given step.
#' @param soma.minAbsoluteSep Self-Organising Migrating Algorithm parameter. Smallest absolute difference between maximum and minimum cost function values. Below this minimum the algorithm will terminate.
#' @param soma.minRelativeSep Self-Organising Migrating Algorithm parameter. Smallest relative difference between maximum and minimum cost function values. Below this minimum the algorithm will terminate.
#' @param soma.nMigrations Self-Organising Migrating Algorithm parameter. Maximum number of migrations to complete.
#' @param soma.populationSize Self-Organising Migrating Algorithm parameter. Population size.
#' @param tabu.iters Number of iterations in the preliminary search of the algorithm.
#' @param tabu.listSize Tabu list size.
#' @examples
#' 
#' data(Credit)
#' 
#' \dontrun{
#' hE <-hybridEnsemble(x=Credit[1:100,names(Credit) != 'Response'],
#'                     y=Credit$Response[1:100],
#'                     RF.ntree=50,
#'                     AB.iter=50,
#'                     NN.size=5,
#'                     NN.decay=0,
#'                     SV.gamma = 2^-15,
#'                     SV.cost = 2^-5,
#'                     SV.degree=2,
#'                     SV.kernel='radial')
#' }                     
#' @references Ballings, M., Vercamer, D., Van den Poel, D., Hybrid Ensemble: Many Ensembles is Better Than One, Forthcoming.
#' @seealso \code{\link{predict.hybridEnsemble}}, \code{\link{importance.hybridEnsemble}}, \code{\link{CVhybridEnsemble}}, \code{\link{plot.CVhybridEnsemble}}, \code{\link{summary.CVhybridEnsemble}}
#' @return A list of class \code{hybridEnsemble} containing the following elements:
#' \item{LR}{Bagged Logistic Regression model}
#' \item{LR.lambda}{Shrinkage parameter}
#' \item{RF}{Random Forest model}
#' \item{AB}{Stochastic AdaBoost model}
#' \item{KF}{Kernel Factory model}
#' \item{NN}{Neural Network model}
#' \item{SV}{Bagged Support Vector Machines model}
#' \item{SB}{A label denoting which sub-ensemble was the single best}
#' \item{KNN.K}{Optimal number of nearest neighbors}
#' \item{x_KNN}{The full data set for finding the nearest neighbors in the deployment phase}
#' \item{y_KNN}{The full response vector to compute the response of the nearest neigbhors}
#' \item{KNN.size}{Size of the nearest neigbhor sub-ensemble}
#' \item{weightsAUTHORITY}{The weights for the authority-based weighting method}
#' \item{combine}{Combination methods used}
#' \item{constants}{A vector denoting which predictors are constants}
#' \item{minima}{Minimum values of the predictors required for preprocessing the data for the Neural Network}
#' \item{maxima}{Maximum values of the predictors required for preprocessing the data for the Neural Network}
#' \item{calibratorLR}{The calibrator for the Bagged Logistic Regression model}
#' \item{calibratorRF}{The calibrator for the Random Forest model}
#' \item{calibratorAB}{The calibrator for the Stochastic AdaBoost model}
#' \item{calibratorKF}{The calibrator for the Kernel Factory model}
#' \item{calibratorNN}{The calibrator for the Neural Network model}
#' \item{calibratorSV}{The calibrator for the Bagged Support Vector Machines model}
#' \item{calibratorKNN}{The calibrator for the Bagged Nearest Neigbhors}
#' \item{xVALIDATE}{Predictors of the validation sample}
#' \item{predictions}{The seperate predictions by the sub-ensembles}
#' \item{yVALIDATE}{Response variable of the validation sample}
#' \item{eval.measure}{The evaluation measure that was used}
#' @author Michel Ballings, Dauwe Vercamer, and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
hybridEnsemble <- function(  x=NULL,
                             y=NULL,
                             combine=NULL,
                             eval.measure='auc',
                             verbose=FALSE,
                             oversample=TRUE,
                             filter= 0.01,
                             LR.size=10,
                             RF.ntree=500,
                             AB.iter=500,
                             AB.maxdepth=3,
                             KF.cp=1,
                             KF.rp=round(log(nrow(x),10)),
                             KF.ntree=500,
                             NN.rang=0.1,
                             NN.maxit=10000,
                             NN.size=c(5,10,20),
                             NN.decay=c(0,0.001,0.01,0.1),
                             NN.ens.size=10,  
                             SV.gamma = 2^(-15:3),
                             SV.cost = 2^(-5:13),
                             SV.degree=c(2,3),
                             SV.kernel=c('radial','sigmoid','linear','polynomial'),
                             SV.size=10,
                             RoF.L=10,
                             KNN.K=c(1:150),
                             KNN.size=10,
                             rbga.popSize = 42,
                             rbga.iters = 300, 
                             rbga.mutationChance = 1/ rbga.popSize,
                             rbga.elitism= max(1, round(rbga.popSize*0.05)) ,
                             DEopt.nP=20,
                             DEopt.nG=300,
                             DEopt.F=0.9314,
                             DEopt.CR=0.6938,
                             GenSA.maxit=300,
                             GenSA.temperature=0.5,
                             GenSA.visiting.param=2.7,
                             GenSA.acceptance.param=-5,
                             GenSA.max.call=1e7,
                             malschains.popsize=60,
                             malschains.ls="cmaes",
                             malschains.istep=300,
                             malschains.effort=0.5,
                             malschains.alpha=0.5,
                             malschains.threshold=1e-08,
                             malschains.maxEvals=300,
                             psoptim.maxit=300,
                             psoptim.maxf=Inf,
                             psoptim.abstol=-Inf,
                             psoptim.reltol=0,
                             psoptim.s=40,
                             psoptim.k=3,
                             psoptim.p=1-(1-1/psoptim.s)^psoptim.k,
                             psoptim.w=1/(2*log(2)),
                             psoptim.c.p=.5+log(2),
                             psoptim.c.g=.5+log(2),
                             soma.pathLength=3,
                             soma.stepLength=0.11,
                             soma.perturbationChance=0.1,
                             soma.minAbsoluteSep=0,
                             soma.minRelativeSep= 0.001,
                             soma.nMigrations=300,
                             soma.populationSize=10,
                             tabu.iters=300,
                             tabu.listSize=c(5:12)
                             ){
  
  
  if (!is.null(combine) && !tolower(combine) %in% tolower(c("rbga","DEopt","GenSA","malschains","psoptim","soma",
                                                           "NNloglik","GINNLS","LHNNLS",'tabu'))) {
    stop("Please check spelling of combine parameter")
  }

  

  
  #ERROR HANDLING
  if (!is.data.frame(x)) stop("x must be a data frame")
  
  if (is.null(x) || is.null(y)) {
		stop("x or y cannot be NULL.")
	}else if (any(is.na(x)) || any(is.na(y))){
		stop("NAs not permitted.")
	}

	if (!is.factor(y)) stop("y must be a factor")
		
	if (any(table(y) == 0)) stop("Cannot have empty classes in y.")
	
	if (length(unique(y)) != 2) stop("Must have 2 classes.")

	if (length(y) != nrow(x)) stop("x and y have to be of equal length.")

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
  options(warn=-1)
  #STEP0 SEPERATE DATA INTO TRAINING AND VALIDATION

    
    trainIND <- .partition(y,p=0.75)[[1]]$train
      
    xTRAIN <- x[trainIND,]
    yTRAIN <- y[trainIND]
    
    xVALIDATE <- x[-trainIND,]
    yVALIDATE <- y[-trainIND]
  
  constants <- sapply(xTRAIN,function(x){all(as.numeric(x[1])==as.numeric(x))})
  if (!is.null(filter)) constants <- sapply(xTRAIN,function(x) length(unique(x))<=2 && any(table(x) <= round(nrow(xTRAIN)*filter)))
  xTRAIN <- xTRAIN[,!constants]
  xVALIDATE <- xVALIDATE[,!constants]
  x <- x[,!constants]

  
  
  
  if (verbose==TRUE) cat('Create base classifiers \n')
  

  
  #BAGGED logistic regression
  if (verbose==TRUE) cat('   Bagged Logistic Regression \n')
  
  LR <-  glmnet(x=data.matrix(xTRAIN), y=yTRAIN, family="binomial")

  
  #cross validate lambda
  aucstore <- numeric()
  for (i in 1:length(LR$lambda) ) {
    predglmnet <- predict(LR,newx=data.matrix(xVALIDATE),type="response",s=LR$lambda[i])
    aucstore[i] <- AUC::auc(roc(as.numeric(predglmnet),yVALIDATE))
    
  }
  
  LR.lambda <- LR$lambda[which.max(aucstore)]
  
  
  LR <- list()
  predLR <- data.frame(matrix(nrow=nrow(xVALIDATE),ncol=LR.size))

  noconst <- FALSE
  for (i in 1:LR.size) {
        #make sure we don't make constants in bootstrapping
        while(noconst==FALSE){
            ind <- sample(nrow(xTRAIN),size=round(1*nrow(xTRAIN)), replace=TRUE)
            const <- sapply(xTRAIN,function(x){all(as.numeric(x[1])==as.numeric(x))})
         	  if (!is.null(filter)) const <- sapply(xTRAIN,function(x) length(unique(x))<=2 && any(table(x) <= round(nrow(xTRAIN)*filter)))
      		  noconst <- sum(const) == 0
        } 
          
        LR[[i]] <-  glmnet(x=data.matrix(xTRAIN[ind,]), y=yTRAIN[ind], family="binomial")
        
                
        #Compute the predictions for weight optimization
        predLR[,i] <- as.numeric(predict(LR[[i]],newx=data.matrix(xVALIDATE),type="response",s=LR.lambda))
  }

  predLR <- as.numeric(rowMeans(predLR))
  calibratorLR <-  .calibrate(x=predLR,y=yVALIDATE)
  predLR <- .predict.calibrate(object=calibratorLR, newdata=predLR)
        

  if (tolower(eval.measure)=='spec') {
          
      evalLR <-  AUC::auc(specificity(as.numeric(predLR),yVALIDATE))
          
  } else if (tolower(eval.measure)=='sens') {
        
      evalLR <-  AUC::auc(sensitivity(as.numeric(predLR),yVALIDATE))
          
  } else  if (tolower(eval.measure)=='auc') { 
  
      evalLR <- AUC::auc(roc(as.numeric(predLR),yVALIDATE))
  }
        
  
  
  #Create final model
  noconst <- FALSE
  for (i in 1:LR.size) {
    #make sure we don't make constants in bootstrapping
    while(noconst==FALSE){
      ind <- sample(nrow(x),size=round(1*nrow(x)), replace=TRUE)
      const <- sapply(x,function(x){all(as.numeric(x[1])==as.numeric(x))})
      if (!is.null(filter)) const <- sapply(x,function(z) length(unique(z))<=2 && any(table(z) <= round(nrow(x)*filter)))
    	noconst <- sum(const) == 0
    }
  
    LR[[i]] <-  glmnet(x=data.matrix(x[ind,]), y=y[ind], family="binomial")
  }
  
  
  ####################################################################################################
  #random forest
  if (verbose==TRUE) cat('   Random Forest \n')
  RF <- randomForest(xTRAIN,as.factor(yTRAIN),  ntree=RF.ntree, importance=FALSE, na.action=na.omit )
  
  #Compute the predictions for weight optimization
  predRF <- as.numeric(predict(RF,xVALIDATE,type="prob")[,2])
  calibratorRF <-  .calibrate(x=predRF,y=yVALIDATE)
  predRF <- .predict.calibrate(object=calibratorRF, newdata=predRF)
  
  if (tolower(eval.measure)=='spec') {
    
    evalRF <- AUC::auc(specificity(as.numeric(predRF),yVALIDATE))
    
  } else if (tolower(eval.measure)=='sens') {
    
    evalRF <- AUC::auc(sensitivity(as.numeric(predRF),yVALIDATE))
    
  } else  if (tolower(eval.measure)=='auc') { 
    
    evalRF <- AUC::auc(roc(as.numeric(predRF),yVALIDATE))
  }
  
 
  
  #Create final model
  RF <- randomForest(x,as.factor(y),  ntree=RF.ntree, importance=FALSE, na.action=na.omit )
  
  ####################################################################################################
  #ada boost
  if (verbose==TRUE) cat('   AdaBoost \n')
  AB <- ada(xTRAIN,as.factor(yTRAIN),iter=AB.iter, control=rpart.control(maxdepth=AB.maxdepth))
  
  #Compute the predictions for weight optimization
  predAB <- as.numeric(predict(AB,xVALIDATE,type="probs")[,2])
  calibratorAB <-  .calibrate(x=predAB,y=yVALIDATE)
  predAB <- .predict.calibrate(object=calibratorAB, newdata=predAB)

  if (tolower(eval.measure)=='spec') {
    
    evalAB <-  AUC::auc(specificity(as.numeric(predAB),yVALIDATE))
    
  } else if (tolower(eval.measure)=='sens') {
    
    evalAB <-  AUC::auc(sensitivity(as.numeric(predAB),yVALIDATE))
    
  } else  if (tolower(eval.measure)=='auc') { 
    
    evalAB <-  AUC::auc(roc(as.numeric(predAB),yVALIDATE))
  }
  
  
  
  #Create final model
  AB <- ada(x,as.factor(y),iter=AB.iter, control=rpart.control(maxdepth=AB.maxdepth))
  
  ####################################################################################################
  #kernelFactory
  if (verbose==TRUE) cat('   Kernel Factory \n')
  KF <- kernelFactory(xTRAIN,as.factor(yTRAIN), rp=KF.rp, cp=KF.cp, ntree=KF.ntree)
  
  #Compute the predictions for weight optimization
  predKF <- as.numeric(predict(KF,xVALIDATE))
  calibratorKF <-  .calibrate(x=predKF,y=yVALIDATE)
  predKF <- .predict.calibrate(object=calibratorKF, newdata=predKF)
  
  if (tolower(eval.measure)=='spec') {
    
    
    evalKF <- AUC::auc(specificity(as.numeric(predKF),yVALIDATE))
    
  } else if (tolower(eval.measure)=='sens') {
    
    evalKF <- AUC::auc(sensitivity(as.numeric(predKF),yVALIDATE))
 
    } else  if (tolower(eval.measure)=='auc') { 
    
    evalKF <- AUC::auc(roc(as.numeric(predKF),yVALIDATE))
  }
  
  
  
  #Create final model
  KF <- kernelFactory(x,as.factor(y), rp=KF.rp, cp=KF.cp, ntree=KF.ntree)
  
  ####################################################################################################
  #BAGGED neural network (version: only tune once)
  if (verbose==TRUE) cat('   Bagged Neural Network \n')
  
  xTRAINscaled <- data.frame(sapply(xTRAIN, as.numeric))
  xTRAINscaled <- data.frame(sapply(xTRAINscaled, function(x) if(length(unique(x))==2 && min(x)==1) x-1 else x))
  
  minima <- sapply(xTRAINscaled,min)
  maxima <- sapply(xTRAINscaled,max)
  
  xTRAINscaled <- data.frame(t((t(xTRAINscaled) - ((minima + maxima)/2))/((maxima-minima)/2)))

  
  call <- call("nnet", formula = yTRAIN ~ ., data=xTRAINscaled,  rang=NN.rang, maxit=NN.maxit, trace=FALSE, MaxNWts= Inf)
  tuning <- list(size=NN.size, decay=NN.decay)

    #tune nnet

  xVALIDATEscaled <- data.frame(sapply(xVALIDATE, as.numeric))
  xVALIDATEscaled <- data.frame(sapply(xVALIDATEscaled, function(x) if(length(unique(x))==2 && min(x)==1) x-1 else x))
  
  xVALIDATEscaled <- data.frame(t((t(xVALIDATEscaled) - ((minima + maxima)/2))/((maxima-minima)/2)))

  result <- .tuneMember(call=call,
                         tuning=tuning,
                         xtest=xVALIDATEscaled,
                         ytest=yVALIDATE,
                         predicttype="raw")
                     
  
  predNN <- data.frame(matrix(nrow=nrow(xVALIDATE),ncol=NN.ens.size))

  nozerorange <- FALSE
  for (i in 1:NN.ens.size) {

     #make sure we don't divide by 0
        while(nozerorange==FALSE){
          ind <- sample(nrow(xTRAIN),size=round(1*nrow(xTRAIN)), replace=TRUE)
      
          xTRAINscaled <- data.frame(sapply(xTRAIN[ind,], as.numeric))
          xTRAINscaled <- data.frame(sapply(xTRAINscaled, function(x) if(length(unique(x))==2 && min(x)==1) x-1 else x))
         
          nozerorange <- all(sapply(xTRAINscaled[ind,],function(x) (max(x)-min(x))!=0))

        }
        
        minima <- sapply(xTRAINscaled,min)
        maxima <- sapply(xTRAINscaled,max)
        
        xTRAINscaled <- data.frame(t((t(xTRAINscaled) - ((minima + maxima)/2))/((maxima-minima)/2)))

             
        #use the optimal parameters to train final model
        NN <- nnet(yTRAIN[ind] ~ ., xTRAINscaled, size = result$size,  rang = NN.rang, decay = result$decay, maxit = NN.maxit, trace=FALSE, MaxNWts= Inf)
        
        #Compute the predictions for weight optimization
        predNN[,i] <- as.numeric(predict(NN,xVALIDATEscaled,type="raw"))
  
  }
  
  
  predNN <- rowMeans(predNN)
  calibratorNN <-  .calibrate(x=predNN,y=yVALIDATE)
  predNN <- .predict.calibrate(object=calibratorNN, newdata=predNN)
  
  if (tolower(eval.measure)=='spec') {
    
    evalNN <- AUC::auc(specificity(as.numeric(predNN),yVALIDATE))
    
  } else if (tolower(eval.measure)=='sens') {
    
    evalNN <- AUC::auc(sensitivity(as.numeric(predNN),yVALIDATE))
    
  } else  if (tolower(eval.measure)=='auc') { 
    
    evalNN <- AUC::auc(roc(as.numeric(predNN),yVALIDATE))
  }
  

  
  #Create final model
  

  minima <- list()
  maxima <- list()
  NN <- list()
  nozerorange <- FALSE
  for (i in 1:NN.ens.size) {
      #make sure we don't divide by 0
      while(nozerorange==FALSE){
        ind <- sample(nrow(x),size=round(1*nrow(x)), replace=TRUE)
          
        xscaled <- data.frame(sapply(x[ind,], as.numeric))
        xscaled <- data.frame(sapply(xscaled, function(x) if(length(unique(x))==2 && min(x)==1) x-1 else x))
        
        nozerorange <- all(sapply(xscaled[ind,],function(x) (max(x)-min(x))!=0))

      }
      minima[[i]] <- sapply(xscaled,min)
      maxima[[i]] <- sapply(xscaled,max)
        
      xscaled <- data.frame(t((t(xscaled) - ((minima[[i]] + maxima[[i]])/2))/((maxima[[i]]-minima[[i]])/2)))

      
      NN[[i]] <- nnet(y[ind] ~ ., xscaled, size = result$size, rang = NN.rang,  decay = result$decay, maxit = NN.maxit, trace=FALSE, MaxNWts= Inf)
  }
  
  
  ####################################################################################################
 
  #BAGGED support vector machines (version: only tune once and not on every fold)
  if (verbose==TRUE) cat('   Bagged Support Vector Machine \n')
  
  result <- list()
  for (i in SV.kernel) {
          call <- call("svm", formula = as.factor(yTRAIN) ~ ., data=xTRAIN, type = "C-classification",probability=TRUE)
          
          if (i=='radial'  ) tuning <- list(gamma = SV.gamma, cost = SV.cost, kernel='radial')
          if (i=='sigmoid') tuning <- list(gamma = SV.gamma, cost = SV.cost, kernel='sigmoid')
          if (i=='linear') tuning <- list(cost = SV.cost, kernel='linear')
          if (i=='polynomial') tuning <- list(gamma = SV.gamma, cost = SV.cost, degree=SV.degree, kernel='polynomial')
          
            
          #tune svm
          result[[i]] <- .tuneMember(call=call,
                               tuning=tuning,
                               xtest=xVALIDATE,
                               ytest=yVALIDATE,
                               probability=TRUE)
          
  }
  auc <- numeric()
  for (i in 1:length(result)) auc[i] <- result[[i]]$auc

  result <- result[[which.max(auc)]]

  
  
  predSV <- data.frame(matrix(nrow=nrow(xVALIDATE),ncol=SV.size))

  noconst <- FALSE
  for (ii in 1:SV.size) {
        #make sure we don't make constants in bootstrapping
        while(noconst==FALSE){
           ind <- sample(nrow(xTRAIN),size=round(1*nrow(xTRAIN)), replace=TRUE)
           const <- sapply(xTRAIN,function(x){all(as.numeric(x[1])==as.numeric(x))})
           if (!is.null(filter)) const <- sapply(xTRAIN,function(x) length(unique(x))<=2 && any(table(x) <= round(nrow(xTRAIN)*filter)))
           noconst <- sum(const) == 0
        
        }   
        
        #use the optimal parameters to train final model
        SV <- svm(as.factor(yTRAIN[ind]) ~ ., data = xTRAIN[ind,],
                      type = "C-classification", kernel = as.character(result$kernel), degree= if (is.null(result$degree)) 3 else result$degree, 
                      cost = result$cost, gamma = if (is.null(result$gamma)) 1 / ncol(xTRAIN) else result$gamma , probability=TRUE)
            
        #Compute the predictions for weight optimization
        predSV[,ii] <- as.numeric(attr(predict(SV,xVALIDATE, probability=TRUE),"probabilities")[,2])
  }
  predSV <- rowMeans(predSV)
  calibratorSV <-  .calibrate(x=predSV,y=yVALIDATE)
  predSV <- .predict.calibrate(object=calibratorSV, newdata=predSV)
  
  if (tolower(eval.measure)=='spec') {
            
    evalSV <-  AUC::auc(specificity(as.numeric(predSV),yVALIDATE))
    
  } else if (tolower(eval.measure)=='sens') {
    
    evalSV <-  AUC::auc(sensitivity(as.numeric(predSV),yVALIDATE))
    
  } else  if (tolower(eval.measure)=='auc') { 
    
    evalSV <-  AUC::auc(roc(as.numeric(predSV),yVALIDATE))
  }
  
  
 
  SV <- list()
  noconst <- FALSE
  #Create final model
  for (i in 1:SV.size) {
      while(noconst==FALSE){
        ind <- sample(nrow(x),size=round(1*nrow(x)), replace=TRUE)
        const <- sapply(x,function(x){all(as.numeric(x[1])==as.numeric(x))})
        if (!is.null(filter)) const <- sapply(x,function(z) length(unique(z))<=2 && any(table(z) <= round(nrow(x)*filter)))
        noconst <- sum(const) == 0
      }
    
      SV[[i]] <- svm(as.factor(y[ind]) ~ ., data = x[ind,],
                    type = "C-classification", kernel = as.character(result$kernel), degree= if (is.null(result$degree)) 3 else result$degree,
                    cost = result$cost, gamma = if (is.null(result$gamma)) 1 / ncol(xTRAIN) else result$gamma, probability=TRUE)
  }        
  
  ###################################################################################################
  #rotation forest


  if (verbose==TRUE) cat('   Rotation Forest \n')

  RoF <- rotationForest(x=xTRAIN,y=as.factor(yTRAIN),L=RoF.L)
  #Compute the predictions for weight optimization
  predRoF <- as.numeric(predict(RoF,xVALIDATE))
  calibratorRoF <-  .calibrate(x=predRoF,y=yVALIDATE)
  predRoF <- .predict.calibrate(object=calibratorRoF, newdata=predRoF)
  
  if (tolower(eval.measure)=='spec') {
    
    evalRoF <- AUC::auc(specificity(as.numeric(predRoF),yVALIDATE))
    
  } else if (tolower(eval.measure)=='sens') {
    
    evalRoF <- AUC::auc(sensitivity(as.numeric(predRoF),yVALIDATE))
    
  } else  if (tolower(eval.measure)=='auc') { 
    
    evalRoF <- AUC::auc(roc(as.numeric(predRoF),yVALIDATE))
  }
  
 
  
  #Create final model
  RoF <- rotationForest(x=x,y=as.factor(y),L=RoF.L)

  ###################################################################################################
  #Bagged K-Nearest Neigbhors
  if (verbose==TRUE) cat('   Bagged K-Nearest Neigbhors \n')
  #the knnx function requires all indicators to be numeric so we first convert our data.
  
  xTRAIN_KNN <- data.frame(sapply(xTRAIN, as.numeric))
  xTRAIN_KNN <- data.frame(sapply(xTRAIN_KNN, function(x) if(length(unique(x))==2 && min(x)==1) x-1 else x))
         
  xVALIDATE_KNN <- data.frame(sapply(xVALIDATE, as.numeric))
  xVALIDATE_KNN <- data.frame(sapply(xVALIDATE_KNN, function(x) if(length(unique(x))==2 && min(x)==1) x-1 else x))
  
  x_KNN <- data.frame(sapply(x, as.numeric))
  x_KNN <- data.frame(sapply(x_KNN, function(x) if(length(unique(x))==2 && min(x)==1) x-1 else x))
  
  KNN.K <- KNN.K[KNN.K <= round(0.5*nrow(xTRAIN_KNN))]
  
  auc <- numeric(length(KNN.K))
  
  #determine optimal k
  kk <- 0
  for (k in KNN.K) {
    kk <- kk + 1
    #retrieve the indicators of the k nearest neighbors of the query data 
    indicatorsKNN <- as.integer(knnx.index(data=xTRAIN_KNN, query=xVALIDATE_KNN, k=k))
    #retrieve the actual y from the tarining set
    predKNN <- as.integer(as.character(yTRAIN[indicatorsKNN]))
    #if k > 1 then we take the proportion of 1s
    predKNN <- rowMeans(data.frame(matrix(data=predKNN,ncol=k,nrow=nrow(xVALIDATE_KNN))))
    
    #COMPUTE AUC 
    auc[kk] <- AUC::auc(roc(predKNN,yVALIDATE))
  
  }
  
  k <- which.max(auc)
  
  #create ensemble predictions using optimal k (needed for weight estimation)
  predKNN <- data.frame(matrix(nrow=nrow(xVALIDATE_KNN),ncol=KNN.size))
  for (i in 1:KNN.size){
    ind <- sample(1:nrow(xTRAIN_KNN),size=round(nrow(xTRAIN_KNN)), replace=TRUE)
    #retrieve the indicators of the k nearest neighbors of the query data 
    indicatorsKNN <- as.integer(knnx.index(data=xTRAIN_KNN, query=xVALIDATE_KNN, k=k))
    #retrieve the actual y from the tarining set
    predKNNoptimal <- as.integer(as.character(yTRAIN[indicatorsKNN]))
    #if k > 1 than we take the proportion of 1s
    predKNN[,i] <- rowMeans(data.frame(matrix(data=predKNNoptimal,ncol=k,nrow=nrow(xVALIDATE_KNN))))
  }
  
  predKNN <- rowMeans(predKNN)
  
  calibratorKNN <-  .calibrate(x=predKNN,y=yVALIDATE)
  predKNN <- .predict.calibrate(object=calibratorKNN, newdata=predKNN)
    
    if (tolower(eval.measure)=='spec') {
              
      evalKNN <-  AUC::auc(specificity(as.numeric(predKNN),yVALIDATE))
      
    } else if (tolower(eval.measure)=='sens') {
      
      evalKNN <-  AUC::auc(sensitivity(as.numeric(predKNN),yVALIDATE))
      
    } else  if (tolower(eval.measure)=='auc') { 
      
      evalKNN <-  AUC::auc(roc(as.numeric(predKNN),yVALIDATE))
    }
    


  
  ###################################################################################################
  #storing objects 

    
  predictions <- data.frame(predLR,predRF,predAB,predKF,predNN,predSV,predRoF,predKNN)
  
  performance <- c(evalLR,evalRF,evalAB,evalKF,evalNN,evalSV,evalRoF,evalKNN)
  performance <- (performance) / sum(performance)
  
  #select single best
  SB <- c("LR","RF","AB","KF","NN","SV","RoF","KN")[which.max(performance)]
  
  result <- list(LR=LR,LR.lambda=LR.lambda, RF=RF,AB=AB,KF=KF,NN=NN,SV=SV,RoF=RoF,SB=SB,KNN.K=k,x_KNN=x_KNN,y_KNN=y,KNN.size=KNN.size, weightsAUTHORITY=performance, 
                combine=combine,constants=constants,minima=minima,maxima=maxima,
                 calibratorLR=calibratorLR, calibratorRF=calibratorRF, calibratorAB=calibratorAB,calibratorKF=calibratorKF,
                 calibratorNN=calibratorNN, calibratorSV=calibratorSV, calibratorRoF=calibratorRoF,calibratorKNN=calibratorKNN, xVALIDATE=xVALIDATE, predictions=predictions, yVALIDATE=yVALIDATE,
                 eval.measure=eval.measure)

                 


  #compute weights for combine
  if (verbose==TRUE) cat('Optimize weights \n') 

  #objective function
  
  predictions <- data.matrix(predictions)

    

          
  evaluate <- function(string = c()) {
            
            stringRepaired <- as.numeric(string)/sum(as.numeric(string))
            
            weightedprediction <- as.numeric(rowSums(t(as.numeric(stringRepaired) * t(predictions))))
            
            

            if (tolower(eval.measure)=='spec') {
              returnVal <- -AUC::auc(specificity(weightedprediction,yVALIDATE))
            } else if (tolower(eval.measure)=='sens') {
              returnVal <- -AUC::auc(sensitivity(weightedprediction,yVALIDATE))
            } else  if (tolower(eval.measure)=='auc') { 
              returnVal <- -AUC::auc(roc(weightedprediction,yVALIDATE))
            }
            returnVal
   }
   


  ########################################################################################################################################
  #genalg package: genetic algorithm
  if (tolower('rbga') %in% tolower(combine)) {
    
    if (verbose==TRUE) cat('   Genetic Algorithm \n') 
  
    
    
    tuning <- list(popSize=rbga.popSize,iters=rbga.iters,mutationChance=rbga.mutationChance,elitism=rbga.elitism)  
    grid <- expand.grid(tuning)
    
    perf <- data.frame(matrix(nrow=nrow(grid),ncol=9+ncol(grid)))
    #note: everywhere in perf, the auc is not only for auc but also for sens and spec
    names(perf) <- c("weightLR","weightRF","weightAB","weightKF","weightNN","weightSV","weightRoF","weightKN","auc")
    colnames(perf)[(sum(!is.na(colnames(perf))) + 1) : (sum(!is.na(colnames(perf))) + length(colnames(grid)))] <- colnames(grid)
    
    for (i in 1:nrow(grid)){
      
      
      
      rbga.results <- rbga(stringMin=rep(0, ncol(predictions)), 
                           stringMax=rep(1, ncol(predictions)), 
                           suggestions = rbind(
                             performance,
                             c(rep((1/ncol(predictions)),ncol(predictions))),
                             c(0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.993),
                             c(0.001,0.001,0.001,0.001,0.001,0.001,0.993,0.001),
                             c(0.001,0.001,0.001,0.001,0.001,0.993,0.001,0.001),
                             c(0.001,0.001,0.001,0.001,0.993,0.001,0.001,0.001),
                             c(0.001,0.001,0.001,0.993,0.001,0.001,0.001,0.001),
                             c(0.001,0.001,0.993,0.001,0.001,0.001,0.001,0.001),
                             c(0.001,0.993,0.001,0.001,0.001,0.001,0.001,0.001),
                             c(0.993,0.001,0.001,0.001,0.001,0.001,0.001,0.001)), 
                           popSize = grid[i,]$popSize, 
                           iters = grid[i,]$iters, 
                           mutationChance = grid[i,]$mutationChance,
                           elitism=grid[i,]$elitism,
                           evalFunc = evaluate)
      
      weightsRBGA <- rbga.results$population[which.min(rbga.results$evaluations),]
      weightsRBGA <- as.numeric(weightsRBGA)/sum(as.numeric(weightsRBGA))
      
      perf[i,] <- c(weightsRBGA,-rbga.results$evaluations[which.min(rbga.results$evaluations)],grid[i,])   
    }
    
    result$weightsRBGA  <- perf[which.max(perf$auc),c("weightLR","weightRF","weightAB","weightKF","weightNN","weightSV","weightRoF","weightKN")]
    
  }
  
  
  ########################################################################################################################################
  #NMOF package: differential evolutionary algorithm

  if (tolower('DEopt') %in% tolower(combine)) {
    
    if (verbose==TRUE) cat('   Differential Evolutionary Algorithm \n')     

    
    initial  <- cbind(t(rbind(c(performance),
                              rev(performance),
                              c(rep((1/ncol(predictions)),ncol(predictions))),
                             c(0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.993),
                             c(0.001,0.001,0.001,0.001,0.001,0.001,0.993,0.001),
                             c(0.001,0.001,0.001,0.001,0.001,0.993,0.001,0.001),
                             c(0.001,0.001,0.001,0.001,0.993,0.001,0.001,0.001),
                             c(0.001,0.001,0.001,0.993,0.001,0.001,0.001,0.001),
                             c(0.001,0.001,0.993,0.001,0.001,0.001,0.001,0.001),
                             c(0.001,0.993,0.001,0.001,0.001,0.001,0.001,0.001),
                             c(0.993,0.001,0.001,0.001,0.001,0.001,0.001,0.001))),
                        array(runif(8 * DEopt.nP), dim = c(8, DEopt.nP-11)))
    
    
    initial <- t(t(initial)/colSums(initial))
    
    
    tuning <- list(nP=DEopt.nP,nG=DEopt.nG,F=DEopt.F,CR=DEopt.CR) 
    grid <- expand.grid(tuning)
    
    perf <- data.frame(matrix(nrow=nrow(grid),ncol=9+ncol(grid)))
    #note: everywhere in perf, the auc is not only for auc but also for sens and spec
    names(perf) <- c("weightLR","weightRF","weightAB","weightKF","weightNN","weightSV","weightRoF","weightKN","auc")
    colnames(perf)[(sum(!is.na(colnames(perf))) + 1) : (sum(!is.na(colnames(perf))) + length(colnames(grid)))] <- colnames(grid)
    
    for (i in 1:nrow(grid)){
     
         algo <- list( nP = grid[i,]$nP, ## population size
                nG = grid[i,]$nG, ## number of generations
                F = grid[i,]$F, ## step size
                CR = grid[i,]$CR,
                min = rep(0, ncol(predictions)),
                max = rep(1, ncol(predictions)),
                initP=initial,
                minmaxConstr=TRUE,
                repair = NULL,
                pen = NULL,
                printBar = FALSE,
                printDetail = FALSE,
                loopOF = TRUE, ## do not vectorise
                loopPen = TRUE, ## do not vectorise
                loopRepair = TRUE, ## do not vectorise
                storeF=FALSE) 

          DEopt.results <- DEopt(OF = evaluate,algo = algo)

          weightsDEOPT <- as.numeric(DEopt.results$xbest)/sum(as.numeric(DEopt.results$xbest))
          perf[i,] <- c(weightsDEOPT,-DEopt.results$OFvalue,grid[i,])
         
  }
  
  result$weightsDEOPT  <- perf[which.max(perf$auc),c("weightLR","weightRF","weightAB","weightKF","weightNN","weightSV","weightRoF","weightKN")]
    
  }
  
  ########################################################################################################################################
  #GenSA package: generalized simulated annealing

  if (tolower('GenSA') %in% tolower(combine)) {
  
  if (verbose==TRUE) cat('   Generalized Simulated Annealing \n')     
 
  
  tuning <- list(maxit=GenSA.maxit,temperature=GenSA.temperature,max.call=GenSA.max.call,visiting.param=GenSA.visiting.param, acceptance.param= GenSA.acceptance.param) 
  grid <- expand.grid(tuning)
  
  perf <- data.frame(matrix(nrow=nrow(grid),ncol=9+ncol(grid)))
  #note: everywhere in perf, the auc is not only for auc but also for sens and spec
  names(perf) <- c("weightLR","weightRF","weightAB","weightKF","weightNN","weightSV","weightRoF","weightKN","auc")
  colnames(perf)[(sum(!is.na(colnames(perf))) + 1) : (sum(!is.na(colnames(perf))) + length(colnames(grid)))] <- colnames(grid)
  
  for (i in 1:nrow(grid)){
    
    GenSA.results <- GenSA(par=rep((1/ncol(predictions)),ncol(predictions)),lower=rep(0, ncol(predictions)),upper=rep(1, ncol(predictions)), 
                           fn=evaluate, control=list(maxit=grid[i,]$maxit,
                                                     temperature=grid[i,]$temperature,
                                                     max.call=grid[i,]$max.call,
                                                     visiting.param=grid[i,]$visiting.param,
                                                     acceptance.param=grid[i,]$acceptance.param))
    
    weightsGENSA <- as.numeric(GenSA.results$par)/sum(as.numeric(GenSA.results$par))
        
    perf[i,] <- c(weightsGENSA,-GenSA.results$value,grid[i,])

  }
  
  result$weightsGENSA  <- perf[which.max(perf$auc),c("weightLR","weightRF","weightAB","weightKF","weightNN","weightSV","weightRoF","weightKN")]
  
  }
  
  ########################################################################################################################################
  #Rmalschains package: memetic algorithm with local search chains
  

  if (tolower('malschains') %in% tolower(combine)) {
  
  if (verbose==TRUE) cat('   Memetic Algorithm with Local Search Chains \n')     
  
  quiet <-function(f){
    return(function(...) {capture.output(w<-f(...));return(w);});
  }
  


  initial  <- rbind(t(cbind(performance,
                            c(rep((1/ncol(predictions)),ncol(predictions))),
                             c(0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.993),
                             c(0.001,0.001,0.001,0.001,0.001,0.001,0.993,0.001),
                             c(0.001,0.001,0.001,0.001,0.001,0.993,0.001,0.001),
                             c(0.001,0.001,0.001,0.001,0.993,0.001,0.001,0.001),
                             c(0.001,0.001,0.001,0.993,0.001,0.001,0.001,0.001),
                             c(0.001,0.001,0.993,0.001,0.001,0.001,0.001,0.001),
                             c(0.001,0.993,0.001,0.001,0.001,0.001,0.001,0.001),
                             c(0.993,0.001,0.001,0.001,0.001,0.001,0.001,0.001))),
                    array(runif(8 * malschains.popsize), dim = c(malschains.popsize-10,8)))
 
  initial <- initial/rowSums(initial)
 
  
  tuning <- list(maxEvals=malschains.maxEvals,
                 popsize=malschains.popsize,
                 ls=as.character(malschains.ls), 
                 istep= malschains.istep, 
                 effort=malschains.effort, 
                 alpha=malschains.alpha,
                 threshold= malschains.threshold,
                 optimum=-1.1)
  grid <- expand.grid(tuning)
  
  perf <- data.frame(matrix(nrow=nrow(grid),ncol=9+ncol(grid)))
  #note: everywhere in perf, the auc is not only for auc but also for sens and spec
  names(perf) <- c("weightLR","weightRF","weightAB","weightKF","weightNN","weightSV","weightRoF","weightKN","auc")
  colnames(perf)[(sum(!is.na(colnames(perf))) + 1) : (sum(!is.na(colnames(perf))) + length(colnames(grid)))] <- colnames(grid)
  
  for (i in 1:nrow(grid)){
    
    
 
   malschains.results <-quiet(malschains)(fn=evaluate, 
                                        lower=rep(0, ncol(predictions)),
                                        upper=rep(1, ncol(predictions)), 
                                        initialpop=initial,
                                        control=malschains.control(popsize=grid[i,]$popsize,
                                                                   ls=as.character(grid[i,]$ls), 
                                                                   istep= grid[i,]$istep, 
                                                                   effort=grid[i,]$effort, 
                                                                   alpha=grid[i,]$alpha,
                                                                   threshold= grid[i,]$threshold,
                                                                    optimum=-1),
                                        maxEvals=grid[i,]$maxEvals)

   
    weightsMALSCHAINS <- as.numeric(malschains.results$sol)/sum(as.numeric(malschains.results$sol))

    perf[i,] <- c(weightsMALSCHAINS,-malschains.results$fitness,grid[i,])
  
  }
  
  
    result$weightsMALSCHAINS <- perf[which.max(perf$auc),c("weightLR","weightRF","weightAB","weightKF","weightNN","weightSV","weightRoF","weightKN")]
    
 }

  ########################################################################################################################################
  #pso package: Particle Swarm Optimization
  
  if (tolower('psoptim') %in% tolower(combine)) {
 
  if (verbose==TRUE) cat('   Particle Swarm Optimization \n')     
    
    
    tuning <- list(maxit= psoptim.maxit, maxf=psoptim.maxf, abstol= psoptim.abstol,reltol=psoptim.reltol,s=psoptim.s,k=psoptim.k,p=psoptim.p,w=psoptim.w,c.p=psoptim.c.p,c.g= psoptim.c.g)
    grid <- expand.grid(tuning)
    
    perf <- data.frame(matrix(nrow=nrow(grid),ncol=9+ncol(grid)))
    #note: everywhere in perf, the auc is not only for auc but also for sens and spec
    names(perf) <- c("weightLR","weightRF","weightAB","weightKF","weightNN","weightSV","weightRoF","weightKN","auc")
    colnames(perf)[(sum(!is.na(colnames(perf))) + 1) : (sum(!is.na(colnames(perf))) + length(colnames(grid)))] <- colnames(grid)
    
    for (i in 1:nrow(grid)){
      
      
      
      psoptim.results <- psoptim(par=rep((1/ncol(predictions)),ncol(predictions)),
                                 fn=evaluate,
                                 lower=rep(0, ncol(predictions)),
                                 upper=rep(1, ncol(predictions)),
                                 control=list(maxit= grid[i,]$maxit, 
                                          maxf=grid[i,]$maxf, 
                                          abstol= grid[i,]$abstol,
                                          reltol=grid[i,]$reltol,
                                          s=grid[i,]$s,
                                          k=grid[i,]$k,
                                          p=grid[i,]$p,
                                          w=grid[i,]$w,
                                          c.p=grid[i,]$c.p,
                                          c.g= grid[i,]$c.g,
                                          type="SPSO2011")
                                 )
      
  
      weightsPSOPTIM <- as.numeric(psoptim.results$par)/sum(as.numeric(psoptim.results$par))
      perf[i,] <- c(weightsPSOPTIM,-psoptim.results$value,grid[i,])
      
    }
    
    result$weightsPSOPTIM  <- perf[which.max(perf$auc),c("weightLR","weightRF","weightAB","weightKF","weightNN","weightSV","weightRoF","weightKN")]
  
  
  }
  
  ########################################################################################################################################
  #soma: Self-Organising Migrating Algorithm

  if (tolower('soma') %in% tolower(combine)) {
  
  if (verbose==TRUE) cat('   Self-Organising Migrating Algorithm \n')    
  
  #suppress info
  setOutputLevel(OL$Warning)
    
  
  tuning <- list(  pathLength=soma.pathLength,
                   stepLength=soma.stepLength,
                   perturbationChance=soma.perturbationChance,
                   minAbsoluteSep=soma.minAbsoluteSep,
                   minRelativeSep=soma.minRelativeSep,
                   nMigrations=soma.nMigrations,
                   populationSize=soma.populationSize) 
  grid <- expand.grid(tuning)
  
  perf <- data.frame(matrix(nrow=nrow(grid),ncol=9+ncol(grid)))
  #note: everywhere in perf, the auc is not only for auc but also for sens and spec
  names(perf) <- c("weightLR","weightRF","weightAB","weightKF","weightNN","weightSV","weightRoF","weightKN","auc")
  colnames(perf)[(sum(!is.na(colnames(perf))) + 1) : (sum(!is.na(colnames(perf))) + length(colnames(grid)))] <- colnames(grid)
  
  for (i in 1:nrow(grid)){
    soma.results <- soma(evaluate, list(min=rep(0, ncol(predictions)),max=rep(1, ncol(predictions))),
                         options = list(  pathLength=grid[i,]$pathLength,
                                          stepLength=grid[i,]$stepLength,
                                          perturbationChance=grid[i,]$perturbationChance,
                                          minAbsoluteSep=grid[i,]$minAbsoluteSep,
                                          minRelativeSep=grid[i,]$minRelativeSep,
                                          nMigrations=grid[i,]$nMigrations,
                                          populationSize=grid[i,]$populationSize)
                         )
    
    weightsSOMA <- soma.results$population[,which.min(soma.results$cost)]
    weightsSOMA <- as.numeric(weightsSOMA)/sum(as.numeric(weightsSOMA))
    perf[i,] <- c(weightsSOMA,- soma.results$cost[soma.results$leader],grid[i,])
      
      
  }
  
  result$weightsSOMA  <- perf[which.max(perf$auc),c("weightLR","weightRF","weightAB","weightKF","weightNN","weightSV","weightRoF","weightKN")]
  
  
  }

  ########################################################################################################################################
  #tabu: Tabu Search Algorithm
  
  if (tolower('tabu') %in% tolower(combine)) {
    
    if (verbose==TRUE) cat('   Tabu Search Algorithm \n')  
    
            BinaryVectorToReal <- function(y) {
              ii=0
              y.paste <- as.character()
              for (i in seq(0,(length(y)-11),by=11) )  {ii <- ii + 1 ;  y.paste[ii] <- paste(y[(i+1):(i+11)],collapse='')}
              real <- .unbinary(y.paste)
              res <- real/sum(real)
              res
            }
              
            
            RealVectorToBinary <- function(x) {
              x <- x*100
              y <- .binary(x,mb=10)
              y <- as.integer(unlist(strsplit(as.character(y), "")))
              y
            }
            
    
            
    
            #tabusearch maximizes, hence we need to change the objective function

              
            evaluateTabu <- function(string = c()) {
                      
                        string <- as.integer(unlist(strsplit(as.character(string), "")))
                        x <- BinaryVectorToReal(y=string)
                          
                        returnVal <- NA
                          
                          
                        weightedprediction <- rowSums(t(as.numeric(x) * t(predictions)))
                    
                        if (tolower(eval.measure)=='spec') {
                          returnVal <- AUC::auc(specificity(weightedprediction,yVALIDATE))
                        } else if (tolower(eval.measure)=='sens') {
                          returnVal <- AUC::auc(sensitivity(weightedprediction,yVALIDATE))
                        } else  if (tolower(eval.measure)=='auc') { 
                          returnVal <- AUC::auc(roc(weightedprediction,yVALIDATE))  
                        }
                        returnVal
             }
            
    
    

            config <- RealVectorToBinary(performance)
            
            
            tuning <- list(  iters=tabu.iters,
                             listSize=tabu.listSize
                          ) 
            
            grid <- expand.grid(tuning)
            
            perf <- data.frame(matrix(nrow=nrow(grid),ncol=9+ncol(grid)))
            #note: everywhere in perf, the auc is not only for auc but also for sens and spec
            names(perf) <- c("weightLR","weightRF","weightAB","weightKF","weightNN","weightSV","weightRoF","weightKN","auc")
            colnames(perf)[(sum(!is.na(colnames(perf))) + 1) : (sum(!is.na(colnames(perf))) + length(colnames(grid)))] <- colnames(grid)
            
            for (i in 1:nrow(grid)){
          
                tabu.results <- tabuSearch(size = length(config), 
                                         iters = grid[i,]$iters,
                                         listSize=grid[i,]$listSize, 
                                         objFunc = evaluateTabu, 
                                         config=config) 
          
                weightsTABU <- tabu.results$configKeep[which.max(tabu.results$eUtilityKeep),]
                    
                weightsTABU <-  BinaryVectorToReal(y=weightsTABU)
                    
                perf[i,] <- c(weightsTABU, tabu.results$eUtilityKeep[which.max(tabu.results$eUtilityKeep)],grid[i,])
                    
                    
            }
            
            result$weightsTABU  <- perf[which.max(perf$auc),c("weightLR","weightRF","weightAB","weightKF","weightNN","weightSV","weightRoF","weightKN")]
    
    
    
  }
  

  ########################################################################################################################################
  #LHNNLS (Code adapted from the SuperLearner package)
  
  if (tolower('LHNNLS') %in% tolower(combine)) {
    #adapted from Superlearner package  
    if (verbose==TRUE) cat('   Lawson-Hanson Non-negative least squares \n')    
    
    # compute coef
    fit.nnls <- nnls(as.matrix(predictions), yVALIDATE)
    initCoef <- coef(fit.nnls)
    initCoef[is.na(initCoef)] <- 0.0
    
    # normalize so sum(coef) = 1 if possible
    result$weightsLHNNLS <- initCoef/sum(initCoef)
      
  
  
  }
  
  ########################################################################################################################################
  #GINNLS (Code adapted from the SuperLearner package)
  
  if (tolower('GINNLS') %in% tolower(combine)) {
    #adapted from Superlearner package
    if (verbose==TRUE) cat('   Goldfarb-Idnani Non-negative least squares \n')    


    # compute coef
    .NNLS <- function(x, y) {
      
      D <- t(x) %*% x
      d <- t(t(y) %*% x)
      A <- diag(ncol(x))
      b <- rep(0, ncol(x))
      fit <- solve.QP(Dmat = D, dvec = d, Amat = t(A), bvec = b, meq=0)
      invisible(fit)
    }
    fit.nnls <- .NNLS(x = as.matrix(predictions), y = as.integer(as.character(yVALIDATE)))
    initCoef <- fit.nnls$solution
    initCoef[initCoef < 0] <- 0.0
    initCoef[is.na(initCoef)] <- 0.0
    
    result$weightsGINNLS <- initCoef/sum(initCoef)
    
  
  

  
  }
  ########################################################################################################################################
  #NNloglik (Code adapted from the SuperLearner package)
  
  
  if (tolower('NNloglik') %in% tolower(combine)) {
    #adapted from Superlearner package
    if (verbose==TRUE) cat('   Non-negative binomial likelihood  \n')    
    
  
    trimLogit <- function(x, trim=0.00001) {
      x[x < trim] <- trim
      x[x > (1-trim)] <- (1-trim)
      foo <- log(x/(1-x))
      return(foo)
    }
  
  
    
    .NNloglik <- function(x, y,  start = rep(0, ncol(x))) {
      # adapted from MASS pg 445
      fmin <- function(beta, X, y) {
        p <- plogis(crossprod(t(X), beta))
        -sum(2 *  ifelse(y, log(p), log(1-p)))
      }
      gmin <- function(beta, X, y) {
        eta <- X %*% beta
        p <- plogis(eta)
        -2 * t( dlogis(eta) * ifelse(y, 1/p, -1/(1-p))) %*% X
      }
      fit <- optim(start, fmin, gmin, X = x, y = y,  method = "L-BFGS-B", lower = 0)
      invisible(fit)
    }
    
    
    tempZ <- trimLogit(predictions)
    fit.nnloglik <- .NNloglik(x = as.matrix(tempZ), y = as.integer(as.character(yVALIDATE)))
    
    initCoef <- fit.nnloglik$par
    initCoef[initCoef < 0] <- 0.0
    initCoef[is.na(initCoef)] <- 0.0
    
    result$weightsNNloglik  <- initCoef/sum(initCoef)
    
 
  
  

    
  }
  ########################################################################################################################################  
  class(result) <- "hybridEnsemble"
  
  options(warn=0)
  result
}