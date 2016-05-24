canonCommonality <-
function(A,B,nofns=1){
################################################################# 
##DESCRIPTION
##Returns commonality data for requested functions

##REQUIRED ARGUMENTS
##A		Matrix containing variable set A
##B		Matrix containing variable set B
 
##OPTIONAL ARGUMENTS
##nofns	number of canonical functions to analyze - default to 1

CCdata<-vector("list",2)
CCdata[[1]]<-canonVariate(A,B,nofns)
CCdata[[2]]<-canonVariate(B,A,nofns)
return(CCdata)
}

