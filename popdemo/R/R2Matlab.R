R2Matlab <-
function(A,noquote=FALSE){
rows<-rep(0,nrow(A))
for(i in 1:nrow(A)){
    rows[i]<-paste(as.vector(A[i,]),collapse=" ")
}
matlabstr<-paste(noquote(rows),collapse=";")
matlabstr<-paste(noquote(c("[",matlabstr,"]")),collapse="")
if(noquote) matlabstr<-noquote(matlabstr)
return(matlabstr)
}

