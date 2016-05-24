simple2SO<-function(x)
{
individualsNo<-dim(x)[1]
indivA<-array("",c(individualsNo,3))
for(i in 1:individualsNo)
{
	indivA[i,1]<-i
	indivA[i,2]<-paste("indiv",i,sep="_")
	indivA[i,3]<-i
}
indiv<-as.data.frame(indivA)
names(indiv)<-c("num","name","label")

variablesNo<-dim(x)[2]
variablesICNo<-dim(x)[2]
	detailsICA<-array("",c(variablesICNo,4))
variablesA<-array("",c(variablesNo,5))
detailsICNo<-0
for(i in 1:variablesNo)
{

	variablesA[i,1]<-i
	variablesA[i,2]<-paste("var",i,sep="_")
	variablesA[i,3]<-paste("var",i,sep="_")
		variablesA[i,4]<-"IC"
		detailsICNo<-detailsICNo+1
		detailsICA[detailsICNo,1]<-0
		detailsICA[detailsICNo,2]<-0
		detailsICA[detailsICNo,3]<-min(x[,i,1])
		detailsICA[detailsICNo,4]<-max(x[,i,2])
		variablesA[i,5]<-detailsICNo
}
variables<-as.data.frame(variablesA)
names(variables)<-c("num","name","label","type","details")
detailsIC<-as.data.frame(detailsICA)
names(detailsIC)<-c("na","nu","min","max")
indivICA<-array(0,c(individualsNo,variablesNo,2))
for(i in 1:individualsNo)
{
	for(j in 1:variablesNo)
	{
      indivICA[i,j,1]<-x[i,j,1]
      indivICA[i,j,2]<-x[i,j,2]
	}
}
resul<-list(individuals=indiv,variables=variables,detailsIC=detailsIC,detailsN=NULL,detailsListNom=NULL,detailsNM=NULL,detailsListNomModif=NULL,indivIC=indivICA,indivN=NULL,indivListNom=NULL,indivNM=NULL,indivListNomModif=NULL)

class(resul)<-"symbolic"
resul

}
