#Clustering example on the "Crabs" dataset :
#Here two models are fitted and different number of clusters 
#are tested. The model and the number of clusters that 
#maximize the BIC criterion are selected.

data(Crabs)
A <- Crabs[,-1]
cls <- Crabs[,1]
prms <- hddc(A,4)
res <- predict(prms,A,cls)

#selection of the number of clusters
prms <- hddc(A,graph=TRUE)
res<-predict(prms,A,cls)

#Now a PCA on the "Crabs" dataset is done on the two first 
#principal axis.
#It is an illustration of the clustering process. 
#The means of the clusters are represented by the points 
#and the directions of each group is represented by a line. 
#All the steps of the algorithm are shown.
#The algorithm, the initialization and the model can be chosen.
demo_hddc()




