############################
# observed counterexamples #
############################

#########################################
#                                       # 
# This function computes from a dataset #
# for all item pairs the corresponding  #
# numbers of counterexamples.           #
#                                       # 
#########################################

ob_counter<-function(dataset){
m<-ncol(dataset)
n<-nrow(dataset)
b<-matrix(0,ncol = m, nrow = m)
for(i in 1:m){
for(j in 1:m){
if(i != j) b[i,j]<-sum(dataset[,i] == 0 & dataset[,j] == 1)
}
}
return(b)
}