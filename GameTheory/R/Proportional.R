Proportional <-
function(E,C,Names=NULL){
	
lambda<-E/sum(C)
R<-lambda*C

Output<-list(Results=R,Claims=C,Method="Proportional Rule",Short="P",E=E,Names=Names)
class(Output)<-"ClaimsRule"
return(Output) 
		
}
