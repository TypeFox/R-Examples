AllRules <-
function(E,C,Names=NULL,pct=0,r=2){
	A<-Proportional(E,C)$Results
	B<-CEA(E,C)$Results
	H<-CEL(E,C)$Results
	D<-Talmud(E,C)$Results
	W<-RandomArrival(E,C)$Results
	Res<-cbind(C,A,B,H,D,W)
	G<-apply(Res,2,Gini)
	
	if (pct==1){Res<-cbind(C,A,B,H,D,W)/C}
	
	
	
	Res[,1]<-C
	Res<-round(Res,r)
	
	
	
	G<-round(G,r)
	Res<-rbind(Res,G)
	
	Res[length(C)+1,1]<-NA
	
	colnames(Res)<-c("Claim","Proportional","CEA","CEL","Talmud","RA")
	
	rownames(Res)<-c(Names,"Gini Index")
	class(Res)<-"ClaimsRules"
     
	return(Res)
}
