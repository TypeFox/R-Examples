RWBP.default <-
function(x,...,nn_k=10,min.clusters=8,clusters.iterations=6,clusters.stepSize=2,alfa=0.5,dumping.factor=0.9){
#library(RANN)
#library(igraph)
#library(lsa)
#library(SnowballC)
RW <- NULL
# remove records with missing data 
 newdata <- na.omit(as.data.frame(x)) 
RW$data <- newdata 
#spatial attributes
RW$X <- newdata [,c(1,2)]
#non-spatial attributes
RW$Y <- newdata [,-c(1,2)]
class(RW) <- "RWBP"
#index records according to their original order
RW$ID <- row(newdata)[,1]
#records amount before and after removing records with empty fields
RW$n <- nrow(newdata)
RW$n.orig <- nrow(as.data.frame(x))
#neighbourhood size parameter
RW$nn_k <- nn_k
#the initial amount of clusters (as used in the first clustering process)
RW$k <- min.clusters
#increase the following clustering process's clusters amount by this size
RW$clusters.stepSize<- clusters.stepSize
#amount of clustering processes 
RW$h <- clusters.iterations
#helps to compute more accurate edge value (distance between object and cluster)
RW$alfa<- alfa
#dumping factor (the probability to return to the original node during each step along a random walk)
RW$c <- dumping.factor 

#Find nn_k nearest neighbours for each spatial object
RW$nearest.indexes <- matrix(c(0),RW$n,RW$nn_k)
for (i in 1:RW$n){
nearest<- nn2(RW$X,RW$X[i,],k=(RW$nn_k + 1),treetype="kd",searchtype="standard")
RW$nearest.indexes[i,] <- nearest[[1]][nearest[[1]]!=i]
}

#Calculate clusters according to numeric attributes.
#the clusters amount will be set to the amount of unique attributes values if it is less than the intended clusters amount
if (!is.vector(RW$Y)) n.unique<- nrow(unique(RW$Y))
else n.unique<- length(unique(RW$Y))
cur_k <- min(RW$k,n.unique)

#clusters is a matrix that holds the clusters achieved in the current iteration(out of h iterations)
#clusteredData holds the spatial attributes in the dataset(X) and the cluster it belongs to according to clustering its numeric attributes(Y)
clusters <- kmeans(RW$Y,cur_k)
clusteredData<- data.frame(RW$ID,RW$Y,clusters$centers[clusters$cluster,],clusters$cluster)
clusteredData$clusters.cluster<- ifelse(clusteredData$clusters.cluster >= 0, -1* clusteredData$clusters.cluster,clusteredData$clusters.cluster)
#perform next clustering processes
if (RW$h > 1) {
 for (i in 2:RW$h) {
  cur_k<- min((cur_k+RW$clusters.stepSize),n.unique)
  clusters<-kmeans(RW$Y,cur_k)
  clusteredData<- rbind(clusteredData,data.frame(RW$ID,RW$Y,clusters$centers[clusters$cluster,],clusters$cluster))
  clusteredData$clusters.cluster<-  ifelse(clusteredData$clusters.cluster >= 0,-1000*(i-1) - clusteredData$clusters.cluster,clusteredData$clusters.cluster)
 }
}
#calculate distance between non-spatial attributes and cluster's centroid
obj.attr.end <- 2
attr.amount <- 1
if (!is.vector(RW$Y)){
obj.attr.end <- 1 + ncol(RW$Y)
attr.amount <- ncol(RW$Y)
}
cluster.attr.start <- obj.attr.end + 1
attributesDist<- 1/exp( abs(clusteredData[,2:obj.attr.end]  - clusteredData[,cluster.attr.start : cluster.attr.start + attr.amount-1])^RW$alfa)
#an extension of the original algorithm to support more than one attribute
if (attr.amount > 1)
clusteredData$avgDist <-  apply(attributesDist,1,mean)
else clusteredData$avgDist <- attributesDist

#building a bipartite graph according to the connections between spatial objects and clusters
clusteredData <- clusteredData[,c(1,ncol(clusteredData)-1,2:(1+attr.amount),ncol(clusteredData))]
RW$clusteredData <- clusteredData
RW$igraph <- graph.data.frame(clusteredData,directed=FALSE)
V(RW$igraph)$type <- V(RW$igraph)$name %in% RW$ID
E(RW$igraph)$weight <- E(RW$igraph)$avgDist #Construct the relation matrix of the bipartite graph

#building adjacency matrix 
M <- as.matrix(RW$igraph, "adjacency") #Construct adjacent matrix
M[is.na(M)] <- 0
M<- apply(M[,],c(1,2),as.numeric)
 #Column normalization (there was an error in formula 5, fixed here)
W<- M/colSums(M)[col(M)] 
W[is.na(W)] <- 0

#Compute similarity vector for each object
#S is a matrix where each row is a similarity vector of a point
I <- matrix(rep(1), nrow(W),ncol(W))
e <- matrix(0,nrow(W),ncol(W))
e[col(e)==row(e)]<- 1
S <- (1-RW$c)*((I-RW$c*W)^(-1))%*%e

#compute the relevance scores between specified object and its neighbours
sim <-matrix(0,RW$n,RW$n)
for (i in 1:RW$n) {
  for (j in 1:RW$nn_k) {
   nb <-RW$nearest.indexes[i,j]
   sim[i,nb] <- cosine(S[i,],S[nb,])
  }
 }
#Compute the Outlierness for each spatial object
RW$OutScore <-matrix(rep(0),RW$n,2)
for (i in 1:RW$n) {
 a <- sim[i,][sim[i,]>0]
 RW$OutScore[i,2] <- prod(a)^(1/RW$nn_k) #calculate geometric mean
 RW$OutScore[i,1] <- i
}
#sort outlier scores ascending according to the score value (so top outliers could be easily extracted)
RW$OutScore<- RW$OutScore[order(RW$OutScore[,2],decreasing=F), ] 
colnames(RW$OutScore) <- c("row_num","outlierScore")
RW$objects.similarity <- sim
RW
}
