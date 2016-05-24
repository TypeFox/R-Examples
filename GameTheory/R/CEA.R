CEA <-
function(E,C,Names=NULL){


n<-length(C)

R<-NULL
R[1]<-min(E/n,C[1])

for (i in 2:n){

	if (R[i-1] >= C[i-1]) {R[i] <- min(C[i], (E - sum(C[1:(i-1)])) / (n+1-i)  )}
	else { R[i] <- R[i-1]  }
	    }
 

 R<-as.matrix(R)
 rownames(R)<-Names
 colnames(R)<-"CEA"
 Output<-list(Results=R,Claims=C,Method="Constrained Equal Awards",Short="CEA",E=E,Names=Names)
 class(Output)<-"ClaimsRule"
 return(Output)   
}
