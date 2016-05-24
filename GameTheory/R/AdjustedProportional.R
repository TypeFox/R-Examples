AdjustedProportional <-
function(E,C,Names=NULL){

#Acabada
	
	# Nuevo claim
	N<-length(C)
	
	Nc<-NULL
	for (i in 1:N){
		R<-max(0,E - sum(C) + C[i])
		Nc<-rbind(Nc,R)
		
	}
	
	print(Nc)
	
	NE <- E - sum(Nc)
	CP <- C - Nc
	
	print(NE)
	print(CP)
	
	NC <- Proportional(NE,CP)$Results
	print(NC)
	print(CP)
	R<- NC + Nc
    
    
	
	Output<-list(Results=R,Claims=C,Method="Adjusted Proportional Rule",Short="AP",E=E,Names=Names)
	class(Output)<-"ClaimsRule"
	return(Output) 
	
}
