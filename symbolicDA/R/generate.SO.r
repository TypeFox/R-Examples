generate.SO<-function(numObjects,numClusters,numIntervalVariables,numMultivaluedVariables){
ile<-numObjects
liczebnosc_klasy<-numClusters
individualsNo<-liczebnosc_klasy*ile
indivA<-array("",c(individualsNo,3))
for(i in 1:individualsNo)
{
	indivA[i,1]<-i
	indivA[i,2]<-paste("indiv",i,sep="_")
	indivA[i,3]<-i
}
indiv<-as.data.frame(indivA)
names(indiv)<-c("num","name","label")

variablesNo<-numMultivaluedVariables+numIntervalVariables
variablesICNo<-numIntervalVariables
variablesNNo<-numMultivaluedVariables
variablesA<-array("",c(variablesNo,5))
if(variablesICNo!=0)
	detailsICA<-array("",c(variablesICNo,4))
else 
	detailsICA<-NULL
if(variablesNNo!=0)
	detailsNA<-array("",c(variablesNNo,3))
else
	detailsNA<-NULL
detailsICNo<-0
detailsNNo<-0
detailsListNomA<-NULL
for(i in 1:variablesNo)
{

	variablesA[i,1]<-i
	variablesA[i,2]<-paste("var",i,sep="_")
	variablesA[i,3]<-paste("var",i,sep="_")
	if (i<=variablesNNo)
	{
		variablesA[i,4]<-"IC"
		detailsICNo<-detailsICNo+1
		detailsICA[detailsICNo,1]<-0
		detailsICA[detailsICNo,2]<-0
		detailsICA[detailsICNo,3]<-1
		detailsICA[detailsICNo,4]<-10
		variablesA[i,5]<-detailsICNo
	}
	else
	{
		variablesA[i,4]<-"MN"
		detailsNNo<-detailsNNo+1
		print(detailsNA)
		print(detailsNNo)
		detailsNA[detailsNNo,1]<-0
		detailsNA[detailsNNo,2]<-0
		detailsNA[detailsNNo,3]<-0
		variablesA[i,5]<-detailsNNo
		for(j in 1:(as.integer(detailsNA[detailsNNo,3])))
		{
			cNom<-c(detailsNNo,"a1","a1","a1")
		if (is.null(detailsListNomA))
			detailsListNomA <- cNom
		else
			detailsListNomA <- rbind(detailsListNomA,cNom)
		}
	}
}
variables<-as.data.frame(variablesA)
names(variables)<-c("num","name","label","type","details")
if (!is.null(detailsICA))
{
	detailsIC<-as.data.frame(detailsICA)
	names(detailsIC)<-c("na","nu","min","max")
}
else
{
	detailsIC<-NULL
}
if (!is.null(detailsNA))
{
	detailsN<-as.data.frame(detailsNA)
	names(detailsN)<-c("na","nu","modals")
}
else
	detailsN<-NULL
if (!is.null(detailsListNomA))
{
	row.names(detailsListNomA)<-NULL
	detailsListNom<-as.data.frame(detailsListNomA)
	names(detailsListNom)<-c("details_no","num","name","label")
}
else
	detailsListNom<-NULL
detailsListNomModifA<-NULL
if (!is.null(detailsListNomModifA))
{
	row.names(detailsListNomModifA)<-NULL
	detailsListNomModif<-as.data.frame(detailsListNomModifA)
	names(detailsListNomModif)<-c("details_no","num","name","label")
}
else
{
  print("detailsListNomModif")
	detailsListNomModif<-NULL
}
indivICA<-array(0,c(individualsNo,variablesNo,2))
indivNNo<-0
indivNA<-NULL
for(i in 1:individualsNo)
{
	for(j in 1:variablesNo)
	{
		if (j==1 && i<=4)
		{
			od<-NULL
			do<-NULL
			for (z in 1:ile)
			{
				od<-cbind(od,rnorm(4,z*4,1))
				do<-cbind(do,rnorm(4,z*4+1,1))
			}
		}
		if (variables[j,"type"]=="IC")
		{
			
			indivICA[i,j,1]<-as.vector(od)[i]
			indivICA[i,j,2]<-as.vector(do)[i]
		}
		if (variables[j,"type"]=="MN")
		{
			kk<-.los(4)
			for(k in 1:kk)
			{
				cNom<-c(i,j,paste("a",.los(3)+4*as.integer(i/4),sep=""))
				if (is.null(indivNA))
					indivNA<-cNom
				else
					indivNA<-rbind(indivNA,cNom)	
			}
		}


	}
}
row.names(indivNA)<-NULL
if (!is.null(indivNA))
{
	indivN<-as.data.frame(indivNA)
	names(indivN)<-c("indiv","variable","value")
}
else
	indivN<-NULL
resul<-list(individuals=indiv,variables=variables,detailsIC=detailsIC,detailsN=detailsN,detailsListNom=detailsListNom,detailsNM=NULL,detailsListNomModif=detailsListNomModif,indivIC=indivICA,indivN=indivN,indivListNom=NULL,indivNM=NULL,indivListNomModif=NULL)

class(resul)<-"symbolic"
resul
}


.los<-function(i)
{
	t<-rnorm(c(0,i),i/2,i/2)
	l<-as.integer(t[1])
	if (t[1]<=1) l<-1

	if (t[1]>=i) l<-i
	l	
}


.is.symbolic<-function(table.Symbolic){
	if (class(table.Symbolic)=="symbolic") is.Symbolic=TRUE else is.Symbolic=FALSE
	is.Symbolic
}


#symbolic.class<-function(individuals,variables,detailsIC,detailsN,detailsListNom,detailsNM,detailsListNomModif,indivIC,indivN,indivListNom,indivNM,indivListNomModif){

.summary.symbolic=function(object,...){
    paste("Number of individuals",nrow(object$individuals),"\n",
    "Number of individuals",nrow(object$variables))
}


#symbolicObjects <-generateSO(30,3,3,3)
#print(summary(symbolicObjects))

