ES.May<-function(mat,digs=3){
	Nk<-apply(mat,1,sum)
        i.names<-paste("i=",1:nrow(mat),sep="")
##calculate Pki matrix        
	Pk<-apply(mat,2,function(x){x/Nk})
	rownames(Pk)<-rownames(mat)
        colnames(Pk)<-i.names
##calculate N matrix            
	N.mat<-matrix(nrow=nrow(mat),ncol=nrow(mat))	
	     for(i in 1:nrow(mat)){
	N.mat[,i]<-rep((1/i),nrow(mat))*(Pk[,i]*Nk)}
	rownames(N.mat)<-rownames(mat)
        colnames(N.mat)<-i.names
##calculate fk vector
	fk.mat<-matrix(nrow=nrow(mat),ncol=nrow(mat))	
	     for(i in 1:nrow(mat)){
	fk.mat[,i]<-rep((1/i),nrow(mat))*Pk[,i]}
        fk.vec<-apply(fk.mat,1,sum)
	names(fk.vec)<-rownames(mat)
	rownames(fk.mat)<-rownames(mat)
        colnames(fk.mat)<-i.names
 
##Nf
	Nk.vec<-fk.vec*Nk
        names(Nk.vec)<-rownames(mat)


	res<-list()
        res$E.S_coefficents<-cbind(N=sum(N.mat),F=sum(N.mat)/sum(Nk))
	res$Nk<-Nk
	res$Pki.matrix<-round(Pk,digs)
        res$N.matrix<-round(N.mat,digs)
        res$fk.matrix<-round(fk.mat, digs)
        res$fk.vector<-round(fk.vec,digs)
        res$Nk.vector<-round(Nk.vec,digs)
        
        res
}