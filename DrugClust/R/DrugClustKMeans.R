#FinalDrugClust with KMeans algorithm
#' DrugClustKMeans
#'
#' Function Implementing DrugClust with KMeans algorithm
#'
#'
#' @name DrugClustKMeans
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
#' #features is the features matrix (see InitFeatures function)
#' # side effects is the matrix containing side effects (see InitSideEffects function)
#' #result<-DrugClustKMeans(num_folds,num_clusters,num_iterations,features,side_effects)
#' @export

DrugClustKMeans<-function(num_folds,num_clusters,num_iterations,features,side_effects){
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

  #KMeans clustering
  cl<-KMeans(train,num_clusters)

  #KMeans clustering method application

  A<-KMeansModel(train,trainpharmat,num_clusters,cl)

  #KSeeds predictions for drugs in the test set
  predizioni<-PredictionKMeans(A,cl,test)

  #Function for obtaining AUC
  vectorAUC<-AUC(predizioni,testpharmat,vectorAUC,"KMeans")

  #Function for obtaining AUPR
  vectorAUPR<-AUPR(predizioni,testpharmat,vectorAUPR,"KMeans")

  for(i in 1:(num_folds-1)){
    train = features[folds != i,]
    trainpharmat = side_effects[folds != i,]
    test = features[folds == i,]
    testpharmat = side_effects[folds == i,]

    #KMeans clustering and model
    cl<-KMeans(train,num_clusters)
    A<-KMeansModel(train,trainpharmat,num_clusters,cl)

    #KSeeds predictions for drugs in the test set
    predizioni<-PredictionKMeans(A,cl,test)

    #Function for obtaining AUC
    vectorAUC<-AUC(predizioni,testpharmat,vectorAUC,"KMeans")

    #Function for obtaining AUPR
    vectorAUPR<-AUPR(predizioni,testpharmat,vectorAUPR,"KMeans")
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
#DrugClustKMeansEnrichment


#' DrugClustKMeansEnrichment
#'
#' Function Implementing DrugClust with KMeans and Enrichment
#'
#'
#' @name DrugClustKMeansEnrichment
#' @param num_clusters number of clusters desired
#' @param features matrix features
#' @param pharmat matrix of side effects
#' @return number of pathways for various clusters
#' @examples
#' #features is the features matrix
#' #resultSeeds<-DrugClustKMeansEnrichment(4,features)
#' @export


DrugClustKMeansEnrichment<-function(num_clusters,features,pharmat){
  source('enrichmentGO.R')

  #Return the list of clusters
  clusters2<-KMeans(features,num_clusters)
  #clusters<-clusters2$cluster
  all_pathways <-Enrichment_Proteins(features,num_clusters,clusters2$cluster)

  return (all_pathways)
}



