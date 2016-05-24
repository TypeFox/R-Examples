canonVariate <-
function(A,B,nofns){

##Perform canonical correlation
canon<-cca(A,B)

##Determine names of variables in set B
Blist<-dimnames(B)[[2]]

##Build output object
CCdata<-vector("list",nofns)

##For each function requested:
##  Get commonality data
for (i in 1:nofns)
{
   CVA<-canon$canvarx[,i]
   data<-cbind(B,CVA)
   CCdata[[i]]<-commonalityCoefficients(data,"CVA", Blist, "F")
}
return(CCdata)
}

