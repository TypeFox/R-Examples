nppbib  <-  function (datafilename){

#######################################################################
# Load the data                                                       #
#######################################################################
### Check if the data file exists:
 if(!file.exists(datafilename))
  stop('File \'',datafilename,'\' not found')

### Open the file:
 inputfile  <-  file(datafilename,"r")
### Read the data:
###  Note that comments (starting with #) are allowed
 rawdata  <-  read.table(inputfile, header=FALSE)
### Close the file:
 close(inputfile)

### Extract the meta variables from the data:
 ### Nproducts is the number of columns:
  Nproducts <- ncol(rawdata)
 ### Njudges is the number of rows:
  Njudges <- nrow(rawdata)
 ### Nperjudgement is the number of non-NA values in each row:
  Nperjudgement <- length(rawdata[1,!is.na(rawdata[1,])])

### Convert scores to ranks (if already ranks this has no effect):
 rankdata <- matrix(NA,nrow=Njudges,ncol=Nproducts) # initialise
 for(i in 1:Njudges){
  r<-rank(rawdata[i,])  # convert to ranks (NAs are included, though)
  for(j in 1:Nproducts)
   if(!is.na(rawdata[i,j])) # preserve NAs in the ranked data
    rankdata[i,j]<-r[j]
 }


#######################################################################
# Compute the statistic                                               #
#######################################################################
### For each product, compute the number of tests of it and its rank sum:
 R      <- vector(length=Nproducts)       # initialise
 Ntests <- vector(length=Nproducts)       # initialise
 for(i in 1:Nproducts){
  Ntests[i] <- sum(!is.na(rankdata[,i]))      # number of tests involving product i
  R[i]      <- sum(rankdata[,i],na.rm=TRUE)   # sum of ranks assigned to product i
 }
### Compute the adjusted sums (aka product effects):
 adjustedSum <- vector(length=Nproducts)  # initialise
 for(i in 1:Nproducts)
  adjustedSum[i] <- sqrt(12/(Nperjudgement+1))*(R[i]-(Nperjudgement+1)*Ntests[i]/2)

### Construct the variance-covariance matrix:
 V <- matrix(0,nrow=Nproducts,ncol=Nproducts) # initialise
 for(i in 1:Njudges)
  for(j in 1:(Nproducts-1))
   if(!is.na(rankdata[i,j]))
    for(k in (j+1):Nproducts)
     if(!is.na(rankdata[i,k])){
      V[j,k] <- V[j,k]-1
      V[k,j] <- V[j,k]
     }
### Fill in the diagonals:
 for(i in 1:Nproducts)
  V[i,i] <- -sum(V[i,])
### Drop the last row and column of the matrix and then compute its inverse:
 Vinverse <- solve(V[1:(Nproducts-1),1:(Nproducts-1)])

### Compute the Skilling-Mack test statistic:
 SM <- adjustedSum[1:(Nproducts-1)]%*%Vinverse%*%t(t(adjustedSum[1:(Nproducts-1)]))
### Find its p-value:
 pValue <- 1-pchisq(SM,Nproducts-1)

### Return some variables of interest:
 list(
  Nproducts=Nproducts,
  Njudges=Njudges,
  Nperjudgement=Nperjudgement,
  rawdata=rawdata,
  rankdata=rankdata,
  varCovarMatrix=V,
  adjustedSum=adjustedSum,
  SMstatistic=SM,
  pValue=pValue
 )

} ### end of function
