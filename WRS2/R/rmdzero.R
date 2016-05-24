rmdzero<-function(x,est=onestep,grp=NA,nboot=500,...){
#
#   Do ANOVA on dependent groups
#   using #   depth of zero among  bootstrap values
#   based on difference scores.
#
#   The data are assumed to be stored in x in list mode
#   or in a matrix. In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, columns correspond to groups.
#
#   grp is used to specify some subset of the groups, if desired.
#   By default, all J groups are used.
#
#   The default number of bootstrap samples is nboot=500
#
if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
if(is.list(x)){
# put the data in an n by J matrix
mat<-matrix(0,length(x[[1]]),length(x))
for (j in 1:length(x))mat[,j]<-x[[j]]
}
if(is.matrix(x))mat<-x
if(!is.na(grp[1])){
mat<-mat[,grp]
}
mat<-elimna(mat) # Remove rows with missing values.
J<-ncol(mat)
jp<-0
Jall<-(J^2-J)/2
dif<-matrix(NA,nrow=nrow(mat),ncol=Jall)
ic<-0
for(j in 1:J){
for(k in 1:J){
if(j<k){
ic<-ic+1
dif[,ic]<-mat[,j]-mat[,k]
}}}
dif<-as.matrix(dif)
#if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#print("Taking bootstrap samples. Please wait.")
data <- matrix(sample(nrow(mat), size = nrow(mat) * nboot, replace = T), nrow = nboot)
bvec <- matrix(NA, ncol = ncol(dif), nrow = nboot)
        
for(j in 1:ncol(dif)) {
    #paste("Working on contrast ",j, "of ",ncol(dif))
    temp <- dif[, j]
    bvec[, j] <- apply(data, 1, rmanogsub, temp, est)
}  #bvec is an nboot by Jm matrix
center<-apply(dif,2,est,...)
bcen<-apply(bvec,2,mean)
cmat<-var(bvec-bcen+center)
zvec<-rep(0,Jall)
m1<-rbind(bvec,zvec)
bplus<-nboot+1
discen<-mahalanobis(m1,center,cmat)
sig.level<-sum(discen[bplus]<=discen[1:nboot])/nboot
list(discen = discen, p.value=sig.level,center=center)
}
