contrastMatrix<-function(nrFocal,model){
np<-switch(model,"1PL"=1,"2PL"=2,"3PL"=3,"3PLc"=2)
C<-matrix(0,nrFocal*np,(nrFocal+1)*np)
for (i in 1:np){
seq.row<-((i-1)*nrFocal+1):(i*nrFocal)
C[seq.row,i]<-1
for (j in 1:nrFocal){
C[seq.row[j],i+j*np]<-(-1)
}
}
return(C)
}


