#' InitFeatures
#'
#' Initialize the features matrix. The data needs to be binary matrices where each row is a
#' drug, and columns represents drugs features. If the element in position ij is 1 it means that the ith drug interacts with the jth element
#' (for example a protein). The same for the matrix where side effects are stored.
#'
#'
#' @name InitFeatures
#' @param namefeatures name of the file where the features are stored.The file needs to be in the same folder where you have the code.
#' @return The matrix containing drugs features
#' @examples
#' #Generate a sample features binary matrix
#' #for example you will find the file bioma2.txt which is a sample file for feature matrix
#' features<-InitFeatures("bio2mat.txt")
#' @export

InitFeatures<-function(namefeatures){
  p0<-paste0("./", namefeatures)
  features <- as.matrix(utils::read.delim(p0, sep="\t"))
  return(features)
}


#' InitSideEffect
#'
#' Initialize the matrix of features and Side Effects
#'
#'
#' @name InitSideEffect
#' @param nameSideEffects  name of the file where the side effects are stored. The format has to be a binary matrix, where the rows are the drugs and columns are the various side effects (1/0 meaning presence or absence of a certain side effect).
#' @return The matrix containing drugs side effects
#' @examples
#' #Generate a sample features binary matrix
#' #for example you will fin the file pharmat.txt which is a sample file of side_effects matrix
#' side_effects<-InitSideEffects("pharmat.txt")
#' @export

InitSideEffects<-function(nameSideEffects){
  p1<-paste0("./", nameSideEffects)
  pharmat  <- as.matrix(utils::read.delim(p1, sep="\t"))
  return(pharmat)
}


#' CreateFolds
#'
#' Create the folds given the features matrix
#'
#'
#' @name CreateFolds
#' @param features is the features matrix that has to be divided in folds for performing cross validation
#' @param num_folds number of folds desired
#' @return folds: the elements divided in folds
#' @examples
#' features<-InitFeatures("bio2mat.txt")
#' folds<-CreateFolds(features,4)
#' @export

CreateFolds<-function(features,num_folds){
  folds = sample(1:nrow(features)%%num_folds)
  return(folds)
}


#' RandomSeedGenerator
#'
#' Initizalize seeds for the KSeeds clustering algorithm
#'
#' @name RandomSeedGenerator
#' @param num_clusters number of clusters desired
#' @param numbrowfeatures number of rows of the features matrix
#' @return s list of seeds
#' @examples
#' features<-InitFeatures("bio2mat.txt")
#' s<-RandomSeedGenerator(4,nrow(features))
#' @export

RandomSeedGenerator<-function(num_clusters,numbrowfeatures){
s<-sample(1:numbrowfeatures,num_clusters,replace=FALSE)
return(s)
}


#' SeedSelection
#'
#' Given the seeds, it creates the submatrix of the features where the rows are just the seeds drugs
#'
#' @name SeedSelection
#' @param features train matrix of features (in the case of k-folding is the matrix of features)
#' @param num_clusters number of clusters desired
#' @param s the list of seeds
#' @return Seed subset of the feature matrix, where rows are the Seed drugs, and columns the relative features
#' @examples
#' features<-InitFeatures("bio2mat.txt")
#' s<-RandomSeedGenerator(num_clusters,nrow(features))
#' Seed<-SeedSelection(features,num_clusters,s)
#' @export

SeedSelection<-function(features,num_clusters,s){
  Seed<-matrix(nrow=num_clusters,ncol=ncol(features))
  Seed<-features[s,]
  return(Seed)
}



#'KSeedsClusters
#'
#' Function Implementing KSeeds. K-Seeds, firstly randomly chooses a number of
#' drugs (renamed Seeds) equal to the number of clusters desired. Then,
#' the other drugs are assigned to a cluster with respect to Hamming Distance between the drug and the seed of a certain cluster. Cluster
#' seeds are not recomputed at each iteration. This allows a speed up in
#' terms of computational complexity and the algorithm terminates when
#' all the drugs have been assigned.
#'
#'
#' @name KSeedsClusters
#' @param train train matrix of features
#' @param num_clusters number of clusters desired
#' @param Seed subset of drugs features matrix, with just the Seeds as rows
#' @param s the seeds of the clusters
#' @return clusters list indicating the cluster to which each drug belongs to
#' @examples
#' features<-InitFeatures("bio2mat.txt")
#' s<-RandomSeedGenerator(num_clusters,nrow(features))
#' Seed<-SeedSelection(features,num_clusters,s)
#' clusters<-KSeedsClusters (features,num_clusters,Seed,s)
#' @export


KSeedsClusters <- function(train,num_clusters,Seed,s){

  clusters<-numeric()
    for(j in 1:nrow(train)){
      min=e1071::hamming.distance(train[j,],Seed[1,])
      cluster=1
        for(i in 2:num_clusters){
          h<-e1071::hamming.distance(train[j,],Seed[i,])
          matrix1<-rbind(train[j,],Seed[i,])
            if(h<min){
              min=h
              cluster=i
            }
        }
      clusters<-c(clusters,cluster)
    }
  return (clusters)
}


#'KSeedsScores
#'
#' Function for obtaining the Bayesian prediction scores using KSeeds clustering
#'
#'
#' @name KSeedsScores
#' @param train train matrix of features
#' @param trainpharmat train matrix of side effects
#' @param num_clusters number of clusters desired
#' @param Seed subset of the features matrix containing only the Seeds drugs
#' @param s the seeds of the clusters
#' @param clusters the list of clusters where the various drugs are
#' @return A matrix containing prediction scores for each cluster
#' @examples
#' features<-InitFeatures("bio2mat.txt")
#' side_effects<-InitiSideEffects("pharmat.txt")
#' Seed<-SeedSelection(features,num_clusters,s)
#' clusters<-KSeedsClusters (features,num_clusters,Seed,s)
#' A<-KSeedsScores(features,side_effects,num_clusters,Seed,s,clusters)
#' @export


KSeedsScores <- function(train,trainpharmat,num_clusters,Seed,s,clusters){
  A <- matrix(nrow=num_clusters,ncol=ncol(trainpharmat))
  SESeeds<-trainpharmat[s,]
  somma<-colSums(trainpharmat)
  vettorepro<-somma/nrow(trainpharmat)
  numeroSostituzioni=0
    for(n in 1:num_clusters){
      l<-which(clusters==n)
      vettoreProva<-numeric()
        if(length(l)>1){
          drugs_cluster<-trainpharmat[l,]
          SommaSideEffect<-colSums(drugs_cluster)
          SommaTotSE<-colSums(trainpharmat)
          NumeroDrugsCluster<-length(l)
            for(y in 1:ncol(trainpharmat)){
              P1=as.numeric(SommaSideEffect[y])/(SommaTotSE[y])
              P2=as.numeric(SommaTotSE[y])/nrow(train)
              P3=NumeroDrugsCluster/nrow(train)
              A[n,y]<-P1*P2/P3
              vettoreProva<-c(vettoreProva,(P1*P2)/P3)
            }
        }
        else{
          A[n,]<-as.numeric(SESeeds[n,])
        }
    }
  return(A)
}

#' PredictionKSeeds
#'
#' Function implementing predictions for uncharacterized drugs
#'
#' @name PredictionKSeeds
#' @param test test drugs features matrix
#' @param Seed matrix of seeds initialize in the KSeed algorithm
#' @param num_clusters number of clusters desired
#' @param A matrix of Naive Bayes predictions scores, result of KSeedsScores function
#' @param numcolsideffects number of sideeffects
#' @return predizioni matrix containing predictions for the various uncharacterized drugs
#' @examples
#' features<-InitFeatures("bio2mat.txt")
#' side_effects<-InitSideEffects("pharmat.txt")
#' i=0
#' train = features[folds != i,]
#' trainpharmat = side_effects[folds != i,]
#' test = features[folds == i,]
#' testpharmat = side_effects[folds == i,]
#' s<-RandomSeedGenerator(num_clusters,nrow(train))
#' Seed<-SeedSelection(train,num_clusters,s)
#' clusters<-KSeedsClusters (train,num_clusters,Seed,s)
#' A<-KSeedsScores(train,trainpharmat,num_clusters,Seed,s,clusters)
#' predizioni<-PredictionKSeeds(test,Seed,num_clusters,A,ncol(side_effects))
#' @export

PredictionKSeeds <-function(test,Seed,num_clusters,A,numcolsideffects){
  predizioni=matrix(nrow=nrow(test),ncol=numcolsideffects)
  clusters_newdrugs<-numeric()
  h<-numeric()
    for(j in 1:nrow(test)){
      somma=test[j,]+Seed[1,]
      l<-length(which(somma==1))
      h<-c(h,l)
      jaccard=l
      assigned_cluster=1
        for(i in 2:num_clusters){
          somma=test[j,]+Seed[i,]
          l<-length(which(somma==1))
            if(l<jaccard){
              jaccard=l
              assigned_cluster=i
            }
          h<-c(h,l)
        }
      #writeLines(paste(assigned_cluster))
      clusters_newdrugs<-c(clusters_newdrugs,assigned_cluster) #lista dei cluster delle nuove drugs
    }
  for(i in 1:nrow(test)){
    num_cluster<-clusters_newdrugs[i]
    vettoreA<-A[num_cluster,]
    predizioni[i,]<-vettoreA
  }
  return(predizioni)
}


#AUC
#' AUC
#'
#' Function Implementing metrics calculation of AUC
#'
#' @name AUC
#' @param predizioni matrix of predictions
#' @param testpharmat matrix of test for the side effects
#' @param vectorAUC empty vector where the AUC values will be saved
#' @param name string stating the name of the clustering Algorithm used, KSeeds, Kmeans or PAM
#' @return vectorAUC vector containing the various AUC values for the various folds. Moreover the function draw the graph of AUC
#' @examples #' #Function for obtaining AUC
#' features<-InitFeatures("bio2mat.txt")
#' side_effects<-InitSideEffects("pharmat.txt")
#' i=0
#' train = features[folds != i,]
#' trainpharmat = side_effects[folds != i,]
#' test = features[folds == i,]
#' testpharmat = side_effects[folds == i,]
#' s<-RandomSeedGenerator(num_clusters,nrow(train))
#' Seed<-SeedSelection(train,num_clusters,s)
#' clusters<-KSeedsClusters (train,num_clusters,Seed,s)
#' A<-KSeedsScores(train,trainpharmat,num_clusters,Seed,s,clusters)
#' predizioni<-PredictionKSeeds(test,Seed,num_clusters,A,ncol(side_effects))
#' vectorAUC<-numeric()
#' vectorAUC<-AUC(predizioni,testpharmat,vectorAUC,"KSeeds")
#' @importFrom graphics plot abline legend
#' @export

AUC<-function(predizioni,testpharmat,vectorAUC,name){
  pred<-ROCR::prediction(as.numeric(predizioni),as.numeric(testpharmat))
  perf<-ROCR::performance(pred,"tpr","fpr")
  AUC<-ROCR::performance(pred,"auc")
  plot( perf, col = "blue",main="AUC Curves")
  abline(coef=c(0,1),col="red")
  if(name=="KSeeds"){
  names<-c("KSeeds","Reference")
  }
  if(name=="KMeans"){
    names<-c("KMeans","Reference")
  }
  if(name=="PAM"){
    names<-c("PAM","Reference")
  }
  legend('topleft', legend=names ,col=c("blue","red"),lwd=1, lty=c(1,1))
  vectorAUC<-c(vectorAUC, as.numeric(AUC@y.values))
  return(vectorAUC)
}


#' AUPR
#'
#' Function Implementing metrics calculation of AUPR
#'
#' @name AUPR
#' @param predizioni matrix of predictions
#' @param testpharmat matrix of test for the Side Effects
#' @param vectorAUPR empty vector to store AUPR
#' @param name name of the clustering algorithm used (KSeeds, KMeans,PAM)
#' @return vectorAUPR vector containing AUPR values for the various folds, the function also draws AUPR graphs
#' @importFrom graphics plot abline legend
#' @examples
#' #Function for obtaining AUC
#' features<-InitFeatures("bio2mat.txt")
#' side_effects<-InitSideEffects("pharmat.txt")
#' i=0
#' train = features[folds != i,]
#' trainpharmat = side_effects[folds != i,]
#' test = features[folds == i,]
#' testpharmat = side_effects[folds == i,]
#' s<-RandomSeedGenerator(num_clusters,nrow(train))
#' Seed<-SeedSelection(train,num_clusters,s)
#' clusters<-KSeedsClusters (train,num_clusters,Seed,s)
#' A<-KSeedsScores(train,trainpharmat,num_clusters,Seed,s,clusters)
#' predizioni<-PredictionKSeeds(test,Seed,num_clusters,A,ncol(side_effects))
#' vectorAUPR<-numeric()
#' vectorAUPR<-AUPR(predizioni,testpharmat,vectorAUPR,"KSeeds")
#' @export


AUPR <- function(predizioni,testpharmat,vectorAUPR,name){
  pred<-ROCR::prediction(as.numeric(predizioni),as.numeric(testpharmat))
  rec<-ROCR::performance(pred, "rec")
  prec<-ROCR::performance(pred, "prec")
  x<-rec@y.values
  y<-prec@y.values
  x<-unlist(x)
  y<-unlist(y)
  x[ is.nan(x) ] <- 0
  y[ is.nan(y) ] <- 0
  AUPR<-MESS::auc(x,y, type = 'linear')
  perf1 <-ROCR::performance(pred, "prec", "rec")
    if(name=="KSeeds"){
    plot(perf1,main="AUPR KSeeds")
    }
    if(name=="KMeans"){
      plot(perf1,main="AUPR KMeans")
    }
    if(name=="PAM"){
      plot(perf1,main="AUPR PAM")
    }
  vectorAUPR<-c(vectorAUPR,AUPR)
  return (vectorAUPR)
}




#' KMeans
#'
#' KMeans clustering algorithm
#'
#' @name KMeansClusteringAlgorithm
#' @param train matrix of train features
#' @param num_clusters number of clusters desired
#' @return cl list containing the clusters ownerships
#' @examples
#' #use the initFeatures to upload train feature matrix
#' num_clusters=4
#' features<-InitFeatures("bio2mat.txt")
#' i=0
#' train = features[folds != i,]
#' test = features[folds == i,]
#' cl<-KMeans(train,num_clusters)
#' @export

KMeans<-function(train,num_clusters){
cl<-cclust::cclust(train,num_clusters,verbose=TRUE,method="kmeans")
return(cl)
}


#' KMeansModel
#'
#' Function finding the Bayesian Model given the KMeans clustering algorithm
#'
#' @name KMeansModel
#' @param train matrix of train features
#' @param trainpharmat matrix of training of side_effects
#' @param num_clusters number of clusters desired
#' @param cl results of the KMeans model clustering function
#' @return A Bayesian matrix of model for predictions, given the KMeans clustering
#' @examples
#' #First call the KMeans function and obtain cl (list of clusters
#' cl<-KMeans(features,num_clusters)
#' features<-InitFeatures("bio2mat.txt")
#' pharmat<-InitSideEffects("pharmat.txt")
#' features<-InitFeatures("bio2mat.txt")
#' side_effects<-InitSideEffects("pharmat.txt")
#' i=0
#' train = features[folds != i,]
#' trainpharmat = side_effects[folds != i,]
#' test = features[folds == i,]
#' testpharmat = side_effects[folds == i,]
#' A<-KMeansModel(train,trainpharmat,4,cl)
#' @export

KMeansModel<-function(train,trainpharmat,num_clusters,cl){
    A <- matrix(nrow=num_clusters,ncol=1339) #matrice ADIJ
    v<-numeric(1339)
    for(n in 1:num_clusters){
      k<-which(cl$cluster==n)
      if(length(k)>1){
        drugs_cluster<-trainpharmat[k,]
        SommaSideEffect<-colSums(drugs_cluster)
        SommaTotSE<-colSums(trainpharmat)
        NumeroDrugsCluster=length(k)
        for(y in 1:1339){
          P1=as.numeric(SommaSideEffect[y])/(SommaTotSE[y])
          P2=as.numeric(SommaTotSE[y])/nrow(train)
          P3=NumeroDrugsCluster/nrow(train)
          A[n,y]<-P1*P2/P3
        }
      }
      else{
        A[n,]<-v
      }
    }
return(A)
}



#' PredictionKMeans
#'
#' Function finding the predictions for the uncharacterized drugs given the KMeans clustering algorithm
#'
#'
#' @name PredictionKMeans
#' @param A Bayesian model given by the application of KMeansModel algorithm
#' @param cl structure of clusters given by the KMeans function
#' @param test test matrix of drugs
#' @return predizioni matrix with a number of rows equal to the number of clusters and a number of columns equal to the features
#' @importFrom stats predict
#' @examples
#' # A will be the result of the previous call of KMeans model funcion
#' #cl will be the result of KMeans function
#' features<-InitFeatures("bio2mat.txt")
#' side_effects<-InitSideEffects("pharmat.txt")
#' i=0
#' train = features[folds != i,]
#' trainpharmat = side_effects[folds != i,]
#' test = features[folds == i,]
#' testpharmat = side_effects[folds == i,]
#' cl<-KMeans(train,4)
#' A<-KMeansModel(train,trainpharmat,4,cl)
#' predizioni<-PredictionKMeans(A,cl,test)
#' @export

PredictionKMeans<-function(A,cl,test){
  predizioni=matrix(nrow=nrow(test),ncol=1339)
  cl$centers[is.nan(cl$centers)]<-0
  ycl<-predict(cl,test,type="both")
      for(i in 1:nrow(test)){
        num_cluster<-ycl$cluster[i]
        vettoreA<-A[num_cluster,]
        predizioni[i,]<-vettoreA
      }
  return(predizioni)
}



#' PAM
#'
#' PAM clustering algorithm
#'
#'
#' @name PAM
#' @param train matrix of train features
#' @param num_clusters number of clusters desired
#' @return pamx structure with various values resulting from PAM clustering algorithm
#' @examples
#' features<-InitFeatures("bio2mat.txt")
#' i=0
#' train = features[folds != i,]
#' pamx<-PAM(train,4)
#' @export


PAM<-function(train,num_clusters){
  pamx <- cluster::pam(train, num_clusters)
  return(pamx)
}


#' PAM_Model
#'
#' PAM clustering algorithm Model
#'
#' @name PAM_Model
#' @param pamx result of pam clustering algorithm
#' @param num_clusters number of clusters desired
#' @param trainpharmat matrix of training for side effects
#' @param train matrix of train features
#' @return A matrix of model for prediction of uncharacterised drugs, given PAM clustering
#' @examples
#' #pamx is the result of the PAM function
#' #trainpharmat is the side effect train matrix
#' #train is the feature train matrix
#' features<-InitFeatures("bio2mat.txt")
#' side_effects<-InitSideEffects("pharmat.txt")
#' i=0
#' train = features[folds != i,]
#' trainpharmat = side_effects[folds != i,]
#' pamx<-PAM(train,4)
#' A<-PAM_Model(pamx,4,trainpharmat,train)
#' @export


PAM_Model<-function(pamx,num_clusters,trainpharmat,train){
    A <- matrix(nrow=num_clusters,ncol=ncol(trainpharmat))
    I <- matrix(nrow=num_clusters,ncol=ncol(trainpharmat))
    v<-numeric(ncol(trainpharmat))
      for(n in 1:num_clusters){
        l<-which(pamx$clustering==n)
          if(length(l)>1){
            drugs_cluster<-trainpharmat[l,]
            SommaSideEffect<-colSums(drugs_cluster)
            I[n,]<-SommaSideEffect
            SommaTotSE<-colSums(trainpharmat)
            NumeroDrugsCluster=length(l)
              for(y in 1:ncol(trainpharmat)){
                P1=as.numeric(SommaSideEffect[y])/(SommaTotSE[y])
                P2=as.numeric(SommaTotSE[y])/nrow(train)
                P3=NumeroDrugsCluster/nrow(train)
                A[n,y]<-P1*P2/P3
              }
          }
          else{
            A[n,]<-v
          }
      }
return(A)
}

#' Prediction_PAM
#'
#' PAM prediction models
#'
#' @name Prediction_PAM
#' @param pamx result of pam clustering algorithm
#' @param A prediction scores matrix
#' @param test test features matrix
#' @param numb_sideEffects number of side effects
#' @return predizioni matrix of predictions given PAM clustering
#' @examples
#' #A is the result of PAM_Model function
#' #pamx comes from the PAM function
#' #test is the feature test matrix
#' features<-InitFeatures("bio2mat.txt")
#' side_effects<-InitSideEffects("pharmat.txt")
#' i=0
#' train = features[folds != i,]
#' test = features[folds == i,]
#' testpharmat = side_effects[folds == i,]
#' trainpharmat = side_effects[folds != i,]
#' pamx<-PAM(train,4)
#' A<-PAM_Model(pamx,4,trainpharmat,train)
#' predizioni<-PredictionPAM(A,pamx,test)
#' @export

PredictionPAM<-function(A,pamx,test,numb_sideEffects){
  predizioni=matrix(nrow=nrow(test),ncol=numb_sideEffects)
  clusters_newdrugs<-numeric()
  di<-numeric()
  h<-numeric()
  for(j in 1:nrow(test)){
    d<-stats::dist(rbind(test[j,], as.numeric(pamx$medoids[1,])),method="euclidean")
    d1<-stats::dist(rbind(test[j,],as.numeric(pamx$medoids[2,])),method="euclidean")
    if(d>d1){
      clusters_newdrugs<-c(clusters_newdrugs,2)
    }
    else{
      clusters_newdrugs<-c(clusters_newdrugs,1)
    }
  }
  for(i in 1:nrow(test)){
    num_cluster<-clusters_newdrugs[i]
    vettoreA<-A[num_cluster,]
    predizioni[i,]<-vettoreA
  }
  return(predizioni)
}


#' Enrichment_Proteins
#'
#' Function Performing Proteins Enrichment using Gene Ontology
#'
#'
#' @name Enrichment_Proteins
#' @param features matrix of features
#' @param num_clusters number of clusters
#' @param clusters clusters returned from the clustering algorithms
#' @return vector_numb_pathway return a vector telling in how many pathways the various clusters are involved
#' @examples
#' #feature is the feature matrix
#' features<-InitFeatures("bio2mat.txt")
#' side_effects<-InitSideEffects("pharmat.txt")
#' i=0
#' train = features[folds != i,]
#' test = features[folds == i,]
#' testpharmat = side_effects[folds == i,]
#' trainpharmat = side_effects[folds != i,]
#' pamx<-PAM(train,4)
#' A<-PAM_Model(pamx,4,trainpharmat,train)
#' predizioni<-PredictionPAM(A,pamx,test)
#' # pamx$clustering gives the list assigning each element to a certain cluster
#' all_pathways<-Enrichment_Proteins(features,4,pamx$clustering)
#' @export

# This function returns a vector with the number of pathways in which the various
#clusters are involved
Enrichment_Proteins<-function(features,num_clusters,clusters){
  CoresProtein<-matrix(0,nrow=num_clusters,ncol=ncol(features))
  vector_numb_elements<-numeric()
  vector_numb_pathways<-numeric()
    for(p in 1:num_clusters){
      a<-which(clusters==p)
      numb_elements<-length(a)
      vector_numb_elements<-c(vector_numb_elements,numb_elements)
    }
    for(l in 1:length(clusters)){
      cluster<-clusters[l]
      CoresProtein[cluster,]<-CoresProtein[cluster,]+features[l,]
    }
    for(k in 1:nrow(CoresProtein)){
      CoresProtein[k,]<-CoresProtein[k,]/vector_numb_elements[k]
    }
  EntrezList <- as.vector(utils::read.delim("./EntrezList.txt", sep="\n"))
  EntrezList<-EntrezList[[1]]
    for(t in 1:nrow(CoresProtein)){
      ClusterI_simo<-(CoresProtein[t,]>0.01)
      index_elements<-which(ClusterI_simo==TRUE)
      InputEnrichment<-EntrezList[index_elements]
      EnrichmentResults<-enrichmentGO(InputEnrichment,0.05,"BP",1,1)
      vector_numb_pathways<-c(vector_numb_pathways,nrow(EnrichmentResults))
    }
  return (vector_numb_pathways)
}



