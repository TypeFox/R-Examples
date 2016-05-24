#partition a relational matrix to look for a set number of groups, minimizing the inner dissimilarity
`relational.clustering` <- function(dist, clusters = 2) {
  c<-clusters
  D<-as.matrix(dist)
  n<-dim(D)[1]
  #this function is only set up to deal with situations where there could be a theoretical minimum of 2 possible samples per group (ie c<=n/2)
  if (c > n/2) return('Please use a smaller number of groups')
  #The initial partitioning is set so that if there is a group with n==1 members, it will run again
  #however, this will only go on for 1000 runs, so it gives it a fair chance but after 1000 tries it stops
  runs <- 0
  while (runs < 1000) {
    if (runs == 999) return('Sorry, but the initial partitioning has failed; this is often due to too many groups being requested')
    runs <- runs + 1
    
    ins<-sample(1:n,c)
    outs<-1:n
    outs<-outs[-ins]
    #set up a identity matrix for the sites, starting with the 'c' initially picked
    id.mat<-rep(0,n)
    #assign the three randos to the 'c' different groups; obviously, if there are a large numbr of groups, it would be a good idea to run this a number of times to get an initial start state with a sample from each group, probably on the order of 1000 times or something; I think this whole thing could be called a 'branch and bound' style search
    id.mat[ins]<-1:c
    #the -1 in the for loop is because a 1X(n-1) matrix is a vector, therfore has only 1 dimension and should becalled appropriately
    for (i in 1:(n-c-1)) {
      minD <- which(D[outs,ins]==min(D[outs,ins]),TRUE) #check me later
      if (dim(minD)[1]>1) minD<-minD[sample(1:dim(minD)[1],1),]
      pick<-outs[minD[1]]
      ins<-c(ins,pick)
      outs<-outs[-minD[1]]
      id.mat[pick]<-id.mat[ins[minD[2]]]
    }
    minD<-which(D[outs,ins]==min(D[outs,ins]))
    if (length(minD)>1) minD<-minD[sample(1:length(minD),1)]
    id.mat[outs]<-id.mat[ins[minD]]
    #at this point, if any of the groups after the initial placement has a size == 1, just quit and start again, 
    #b/c it's easier tan trying to deal with a group of 1
    if (all(table(id.mat)>1)) runs <- 1000
  }

  #once the initial divisions are made, we can evaluate which localities are poorly placed and start switching those around; 
  #probably do this by seeing which ones have the higher avg in group vs out group values.
  #this matrix holds the values for the ins and outs, just simply as groups 1 through 'c' for each of the columns
  #we're gonna recalculate and do this a bunch of times, once fro every time we switch a locality's group (in id.mat)
  id.mat <- localoptima(dist, id.mat)
  #as a way to evaluate if internal connections are decreasing vs external
  return(id.mat)
}


