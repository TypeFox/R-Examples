CEL <-
function(E,C,Names=NULL){

n<-length(C)

L<-sum(C)-E

R<-NULL
R[1]<-min(L/n,C[1])

for (i in 2:n){

	if (R[i-1] >= C[i-1]) {R[i] <- min(C[i], (L - sum(C[1:(i-1)])) / (n+1-i)  )}
	else { R[i] <- R[i-1]  }
	    }
 
 R<-C-R
  
 

 R<-as.matrix(R)
 rownames(R)<- Names
 colnames(R)<-"CEL"
 Output<-list(Results=R,Claims=C,Method="Constrained Equal Losses",Short="CEL",E=E,Names=Names)
 class(Output)<-"ClaimsRule"
 return(Output) 
 
}
