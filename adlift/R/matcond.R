"matcond" <-
function(x,f,Pred,neigh,int,clo,keep){

a<-fwtnp(x,f,LocalPred=Pred,neighbours=neigh,intercept=int,closest=clo,nkeep=keep,varonly=FALSE,do.W=TRUE)

W<-a$W

cno<-condno(W,type="F")
 
v<-svd(W)$d[1]/(svd(W)$d[nrow(W)])

return(list(cno=cno,v=v,a=a))

}
