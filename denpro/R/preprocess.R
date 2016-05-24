preprocess<-function(dendat, type="copula")
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]
prodendat<-matrix(0,n,d)

if (type=="sphering"){

   cova<-cov(dendat)
   eig<-eigen(cova,symmetric=TRUE)
   sigsqm<-eig$vectors%*%diag(eig$values^(-1/2)) 
   prodendat<-t(t(sigsqm)%*%t(dendat-mean(dendat)))   # dendat%*%sigsqm 

}
else if (type=="sd"){
   for (ii in 1:d){
        prodendat[,ii]<-(dendat[,ii]-mean(dendat[,ii]))/sd(dendat[,ii])
   }

}
else if (type=="standardcopula"){
   for (ii in 1:d){
        or<-order(dendat[,ii])
        mones<-matrix(0,n,1)
        for (i in 1:n) mones[or[i]]<-i
        prodendat[,ii]<-mones/n
   }

}
else{
   for (ii in 1:d){
        or<-order(dendat[,ii])
        mones<-matrix(0,n,1)
        for (i in 1:n) mones[or[i]]<-i
        prodendat[,ii]<-qnorm(mones/n)
   }
}

return(prodendat)
}
