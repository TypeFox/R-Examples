#FinalDrugClust


#' DrugClustKSeeds
#'
#' Function Implementing metrics calculation DrugClust
#'
#'
#' @name DrugClustKSeeds
#' @param num_folds number of folds
#' @param num_clusters number of clusters
#' @param num_iterations number of iterations
#' @param features features matrix
#' @param side_effects side_effects matrix
#' @return (list(AUCFinal,AUPRFinal)) first value is the mean AUC on the various folders, second value is the mean AUPR on the various folders
#' @examples
#' # num_folds=3
#' # num_clusters=4
#' # num_iterations= 5
#' # features is the features matrix (see InitFeatures function)
#' # side effects is the matrix containing side effects (see InitSideEffects function)
#' #result<-DrugClustKSeeds(num_folds,num_clusters,num_iterations,features,side_effects)
#' @export




DrugClustKSeeds<-function(num_folds,num_clusters,num_iterations,features,side_effects){
#vectors where the means of the AUC will be stored
vectorMedieAUC<-numeric()

#vectors where the means of the AUPR will be stored
vectorMedieAUPR<-numeric()

#various Iterations of the k-folds cross validation procedure
for(j in 1:num_iterations){
  i=0 # variable to trace that is the first time of the k-fold cross validation procedure
  vectorAUC<-numeric() #instantiate a vecotr of AUC for first iteration of the k-fold cross validation procedure
  vectorAUPR<-numeric() # vector of AUPR for first iteration of the k-folds cross validation procedure

  folds<-CreateFolds(features,num_folds)
  train = features[folds != i,]
  trainpharmat = side_effects[folds != i,]
  test = features[folds == i,]
  testpharmat = side_effects[folds == i,]

  #KSeeds clustering method application
  s<-RandomSeedGenerator(num_clusters,nrow(train)) #function that generates randomly the numbers that will be the seeds of the cluster
  Seed<-SeedSelection(train,num_clusters,s)
  #Return the list of clusters
  clusters<-KSeedsClusters (train,num_clusters,Seed,s)
  #Return the matrix A of KSeeds scores
  A<-KSeedsScores(train,trainpharmat,num_clusters,Seed,s,clusters)

  #KSeeds predictions for drugs in the test set
  predizioni<-PredictionKSeeds(test,Seed,num_clusters,A)

  #Function for obtaining AUC
  vectorAUC<-AUC(predizioni,testpharmat,vectorAUC,"KSeeds")

  #Function for obtaining AUPR
  vectorAUPR<-AUPR(predizioni,testpharmat,vectorAUPR,"KSeeds")

  for(i in 1:(num_folds-1)){
    train = features[folds != i,]
    trainpharmat = side_effects[folds != i,]
    test = features[folds == i,]
    testpharmat = side_effects[folds == i,]

    #KSeeds clustering
    s<-RandomSeedGenerator(num_clusters,nrow(train)) #function that generates randomly the numbers that will be the seeds of the cluster
    Seed<-SeedSelection(train,num_clusters,s)

    #Return the list of clusters
    clusters<-KSeedsClusters (train,num_clusters,Seed,s)
    #Return the matrix A of KSeeds scores
    A<-KSeedsScores(train,trainpharmat,num_clusters,Seed,s,clusters)





    #KSeeds predictions for drugs in the test set
    predizioni<-PredictionKSeeds(test,Seed,num_clusters,A)

    #Function for obtaining AUC
    vectorAUC<-AUC(predizioni,testpharmat,vectorAUC,"KSeeds")

    #Function for obtaining AUPR
    vectorAUPR<-AUPR(predizioni,testpharmat,vectorAUPR,"KSeeds")
  }
  vectorMedieAUC<-c(vectorMedieAUC,mean(vectorAUC))
  vectorMedieAUPR <-c(vectorMedieAUC,mean(vectorAUPR))
  j=j+1

}

AUCFinal<-mean(vectorMedieAUC)
AUPRFinal<-mean(vectorMedieAUPR)
return(list(AUCFinal,AUPRFinal))

}

#******************************************************************
#FinalDrugClustKSeedsEnrichment


#' DrugClustKSeedsEnrichment
#'
#' Function Implementing DrugClust with KSeeds and Enrichment
#'
#'
#' @name DrugClustKSeedsEnrichment
#' @param num_clusters number of clusters
#' @param features matrix of features
#' @param pharmat matrix of side effects
#' @return number of pathways for various clusters
#' @examples
#' #features is the features matrix
#' #resultSeeds<-DrugClustKSeedsEnrichment(4,features)
#' @export


DrugClustKSeedsEnrichment<-function(num_clusters,features,pharmat){
  source('enrichmentGO.R')
  #KSeeds clustering method application
  s<-RandomSeedGenerator(num_clusters) #function that generates randomly the numbers that will be the seeds of the cluster
  Seed<-SeedSelection(features,num_clusters,s)
  #Return the list of clusters
  clusters<-KSeedsClusters (features,num_clusters,Seed,s)
  all_pathways <-Enrichment_Proteins(features,num_clusters,clusters)

  return (all_pathways)
}

#*******************************************

