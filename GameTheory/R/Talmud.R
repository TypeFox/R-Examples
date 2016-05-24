Talmud <-function(E,C,Names=NULL){

R<-NULL
	if (E <= 0.5*sum(C)){ R <- CEA(E,C/2)[[1]]}
	else {J <- CEL(E - 0.5*sum(C),C/2)
		  R <- C/2 + as.matrix(J[[1]])
		  }
R<-as.matrix(R)
colnames(R)<-"Talmud"
rownames(R)<-Names

Output<-list(Results=R,Claims=C,Method="Proportional Rule",Short="T",E=E,Names=Names)
class(Output)<-"ClaimsRule"
return(Output) 
	
}
