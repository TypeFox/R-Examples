prunelev<-function(bt,lambda=NULL,n=NULL){

len<-length(bt$mean)
bt$mean<-matrix(1,len,1)

if (!is.null(lambda)) bt$ssr<-exmavec(bt$volume,bt$nelem,n,lambda)

ini1<-preprocess(bt$ssr,bt$left,bt$right,bt$mean)
bt$S<-t(ini1$S)
bt$mean<-t(ini1$mean)
bt$left<-ini1$left
bt$right<-ini1$right


#ini<-initial(bt$ssr,bt$left,bt$right)
#bt$left<-ini$left
#bt$right<-ini$right

treeseq<-pruseqlev(bt)

return(treeseq)
}


