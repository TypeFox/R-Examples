View_objective<-function(fba_object){
ix<-which(fba_object$mat[,which(fba_object$obj==1)]!=0)
neg_ix<-which(fba_object$mat[ix,which(fba_object$obj==1)]<0)
pos_ix<-which(fba_object$mat[ix,which(fba_object$obj==1)]>0)

message("Objective Function Substrates")
print(paste(fba_object$metabolite_name[ix[neg_ix]],fba_object$mat[ix[neg_ix],which(fba_object$obj==1)],sep=" -> "))

message("Objective Function Products")
print(paste(fba_object$metabolite_name[ix[pos_ix]],fba_object$mat[ix[pos_ix],which(fba_object$obj==1)],sep=" -> "))
}
