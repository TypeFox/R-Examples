parse.SO<-function(file){
dom<-xmlTreeParse(paste(file,".xml",sep=""))
root<-xmlRoot(dom)
containsXML<-root[["contains"]]
##print(containsXML)
if (xmlAttrs(containsXML)[["HEADER"]]=="YES")
{
headerXML<-root[["header"]]
}
individualsXML<-root[["individus"]]
individualsNo<-as.integer(xmlAttrs(headerXML)[["indiv_nb"]])
indivA<-array("",c(individualsNo,3))
##print(individualsNo)
for(i in 1:individualsNo)
{
	indivA[i,1]<-as.integer(xmlValue(individualsXML[[i]][[1]]))
	indivA[i,2]<-xmlValue(individualsXML[[i]][[2]])
	indivA[i,3]<-xmlValue(individualsXML[[i]][[3]])
}
indiv<-as.data.frame(indivA)
names(indiv)<-c("num","name","label")

variablesXML<-root[["variables"]]
variablesNo<-as.integer(xmlAttrs(headerXML)[["var_nb"]])
variablesA<-array("",c(variablesNo,5))

variablesICNo<-0
variablesCNo<-0
variablesNNo<-0
variablesNMNo<-0
for(i in 1:variablesNo)
{
	##print(xmlValue(variablesXML[[i]][[1]][[1]]))
	##print(variablesXML[[i]])
	variablesA[i,1]<-as.integer(xmlValue(variablesXML[[i]][[1]][[1]]))
	variablesA[i,2]<-xmlValue(variablesXML[[i]][[1]][[2]])
	variablesA[i,3]<-xmlValue(variablesXML[[i]][[1]][[3]])
	if (!is.null(variablesXML[[i]][["mult_nominal"]])  && is.null(variablesXML[[i]][["mult_nominal_Modif"]]))
	{
		variablesNNo<-variablesNNo+1
	}
	if (!is.null(variablesXML[[i]][["inter-cont"]]))
	{
		variablesICNo<-variablesICNo+1
	}
	if (!is.null(variablesXML[[i]][["continue"]]))
	{
		variablesCNo<-variablesCNo+1
	}
	if (!is.null(variablesXML[[i]][["nominal"]]))
	{
		variablesNNo<-variablesNNo+1
	}
	if (!is.null(variablesXML[[i]][["mult_nominal_Modif"]]))
	{
		variablesNMNo<-variablesNMNo+1
	}
}





variablesA<-array("",c(variablesNo,5))
#print(paste("variablesNNo",variablesNNo))
#print(paste("variablesNMNo",variablesNMNo))
if(variablesICNo!=0)
	detailsICA<-array("",c(variablesICNo,4))
else 
	detailsICA<-NULL
if(variablesCNo!=0)
	detailsCA<-array("",c(variablesCNo,4))
else 
	detailsCA<-NULL
if(variablesNNo!=0)
	detailsNA<-array("",c(variablesNNo,3))
else
	detailsNA<-NULL
if(variablesNMNo!=0){
	detailsNMA<-array("",c(variablesNMNo,3))
	dim(detailsNMA)<-c(variablesNMNo,3)
	}
else{
	detailsNMA<-NULL
}

#print(paste("detailsNMA",detailsNMA))
detailsICNo<-0
detailsCNo<-0
detailsNNo<-0
detailsNMNo<-0
#print(paste("detailsNMNo",detailsNMNo))

detailsListNomA<-NULL
detailsListNomModifA<-NULL
for(i in 1:variablesNo)
{
	##print(xmlValue(variablesXML[[i]][[1]][[1]]))
	##print(variablesXML[[i]])
	variablesA[i,1]<-as.integer(xmlValue(variablesXML[[i]][[1]][[1]]))
	variablesA[i,2]<-xmlValue(variablesXML[[i]][[1]][[2]])
	variablesA[i,3]<-xmlValue(variablesXML[[i]][[1]][[3]])
	if (!is.null(variablesXML[[i]][["mult_nominal"]])  && is.null(variablesXML[[i]][["mult_nominal_Modif"]]))
	{
		variablesA[i,4]<-"MN"
		detailsNNo<-detailsNNo+1
		t<-variablesXML[[i]][["mult_nominal"]]
		detailsNA[detailsNNo,1]<-as.integer(xmlAttrs(t[["nominal-desc"]])["nbna"])
		detailsNA[detailsNNo,2]<-as.integer(xmlAttrs(t[["nominal-desc"]])["nbnu"])
		detailsNA[detailsNNo,3]<-as.integer(xmlAttrs(t[["nominal-desc"]])["nbmoda"])
		variablesA[i,5]<-detailsNNo
		for(j in 1:(as.integer(detailsNA[detailsNNo,3])))
		{
			t1<-t[["nominal-desc"]]
			cNom<-c(detailsNNo,t1[j]$"list-nom"[1]$num[1]$text$value,t1[j]$"list-nom"[2]$name[1]$text$value,t1[j]$"list-nom"[3]$label[1]$text$value)
		if (is.null(detailsListNomA))
			detailsListNomA <- cNom
		else
			detailsListNomA <- rbind(detailsListNomA,cNom)
		}
	}
	if (!is.null(variablesXML[[i]][["inter-cont"]]))
	{
		variablesA[i,4]<-"IC"
		detailsICNo<-detailsICNo+1
		t<-variablesXML[[i]][["inter-cont"]]
		detailsICA[detailsICNo,1]<-as.integer(xmlAttrs(t[["continue-desc"]])["nbna"])
		detailsICA[detailsICNo,2]<-as.integer(xmlAttrs(t[["continue-desc"]])["nbnu"])
		detailsICA[detailsICNo,3]<-as.integer(xmlAttrs(t[["continue-desc"]])["min"])
		detailsICA[detailsICNo,4]<-as.integer(xmlAttrs(t[["continue-desc"]])["max"])
		variablesA[i,5]<-detailsICNo
	}
	if (!is.null(variablesXML[[i]][["continue"]]))
	{
		variablesA[i,4]<-"C"
		detailsCNo<-detailsCNo+1
		t<-variablesXML[[i]][["continue"]]
		detailsCA[detailsCNo,1]<-as.integer(xmlAttrs(t[["continue-desc"]])["nbna"])
		detailsCA[detailsCNo,2]<-as.integer(xmlAttrs(t[["continue-desc"]])["nbnu"])
		detailsCA[detailsCNo,3]<-as.integer(xmlAttrs(t[["continue-desc"]])["min"])
		detailsCA[detailsCNo,4]<-as.integer(xmlAttrs(t[["continue-desc"]])["max"])
		variablesA[i,5]<-detailsCNo
	}
	if (!is.null(variablesXML[[i]][["nominal"]]))
	{
		variablesA[i,4]<-"N"
		detailsNNo<-detailsNNo+1
		t<-variablesXML[[i]][["nominal"]]
		detailsNA[detailsNNo,1]<-as.integer(xmlAttrs(t[["nominal-desc"]])["nbna"])
		detailsNA[detailsNNo,2]<-as.integer(xmlAttrs(t[["nominal-desc"]])["nbnu"])
		detailsNA[detailsNNo,3]<-as.integer(xmlAttrs(t[["nominal-desc"]])["nbmoda"])
		variablesA[i,5]<-detailsNNo
		for(j in 1:(as.integer(detailsNA[detailsNNo,3])))
		{
			t1<-t[["nominal-desc"]]
			cNom<-c(detailsNNo,t1[j]$"list-nom"[1]$num[1]$text$value,t1[j]$"list-nom"[2]$name[1]$text$value,t1[j]$"list-nom"[3]$label[1]$text$value)
		if (is.null(detailsListNomA))
			detailsListNomA <- cNom
		else
			detailsListNomA <- rbind(detailsListNomA,cNom)
		}
	}
	if (!is.null(variablesXML[[i]][["mult_nominal_Modif"]]))
	{
		#print("mult_nominal_Modif")
    #print(paste("detailsNMNo",detailsNMNo))
    #print(dim(detailsNMA))
    #print(detailsNMA)
		
		variablesA[i,4]<-"NM"
		detailsNMNo<-detailsNMNo+1
    #print(paste("detailsNMNo",detailsNMNo))
		t<-variablesXML[[i]][["mult_nominal_Modif"]]
		#print(as.integer(xmlAttrs(t[["nominal-desc"]])["nbna"]))
		##print(as.integer(xmlAttrs(t[["nominal-desc"]])["nbna"]))
		detailsNMA[detailsNMNo,1]<-as.integer(xmlAttrs(t[["nominal-desc"]])["nbna"])
		detailsNMA[detailsNMNo,2]<-as.integer(xmlAttrs(t[["nominal-desc"]])["nbnu"])
		detailsNMA[detailsNMNo,3]<-as.integer(xmlAttrs(t[["nominal-desc"]])["nbmoda"])
		variablesA[i,5]<-detailsNMNo
		for(j in 1:(as.integer(detailsNMA[detailsNMNo,3])))
		{
			t1<-t[["nominal-desc"]]
			cNom<-c(detailsNMNo,t1[j]$"list-nom"[1]$num[1]$text$value,t1[j]$"list-nom"[2]$name[1]$text$value,t1[j]$"list-nom"[3]$label[1]$text$value)
			##print(cNom)
		if (is.null(detailsListNomModifA))
			detailsListNomModifA<- cNom
		else
			detailsListNomModifA <- rbind(detailsListNomModifA,cNom)
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
	detailsIC<-NULL
if (!is.null(detailsCA))
{
	detailsC<-as.data.frame(detailsCA)
	names(detailsC)<-c("na","nu","min","max")
}
else
	detailsC<-NULL
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
if (!is.null(detailsNMA))
{
	detailsNM<-as.data.frame(detailsNMA)
	names(detailsNM)<-c("na","nu","modals")
}
else
	detailsNM<-NULL

if (!is.null(detailsListNomModifA))
{
	row.names(detailsListNomModifA)<-NULL
	detailsListNomModif<-as.data.frame(detailsListNomModifA)
	names(detailsListNomModif)<-c("details_no","num","name","label")
}
else
	detailsListNomModif<-NULL

indivICA<-array(0,c(individualsNo,variablesNo,2))
indivCA<-array(0,c(individualsNo,variablesNo,1))
indivNNo<-0
indivNA<-NULL
indivNMNo<-0
indivNMA<-NULL
indivMatXML<-root[["indiv_mat"]]
for(i in 1:individualsNo)
{
	for(j in 1:variablesNo)
	{
		if (variables[j,"type"]=="IC")
		{
			t<-indivMatXML[i]$ligmat[j]$valmat[1]
			indivICA[i,j,1]<-as.integer(t$val_interv[1]$pmin[1]$text$value)
			indivICA[i,j,2]<-as.integer(t$val_interv[2]$pmax[1]$text$value)			
		}
		if (variables[j,"type"]=="C")
		{
			t<-indivMatXML[i]$ligmat[j]$valmat[1]
      if(!is.null(t$val_conti)){
			indivCA[i,j,1]<-as.integer(t$val_conti[1]$text$value)
			}
		}
		if (variables[j,"type"]=="MN")
		{
			for(k in 1:xmlSize(indivMatXML[i]$ligmat[j]$valmat))
			{
				t<-indivMatXML[i]$ligmat[j]$valmat[k]
				#print(paste("val_modal",t$val_modal[1]$text$value))
				cNom<-c(i,j,t$val_modal[1]$text$value)
				if (is.null(indivNA))
					indivNA<-cNom
				else
					indivNA<-rbind(indivNA,cNom)	
			}
		}

		if (variables[j,"type"]=="NM")
		{
			for(k in 1:xmlSize(indivMatXML[i]$ligmat[j]$valmat))
			{
				t<-indivMatXML[i]$ligmat[j]$valmat[k]
				##print(as.double(t$val_list_modal[2]$frequency[1]$text$value))
         #print(paste("no_moda",as.integer(t$val_list_modal[1]$no_moda[1]$text$value)))
         if(!is.null(t$val_list_modal[1]$no_moda[1]$text$value)){
				cNom<-c(i,j,as.integer(t$val_list_modal[1]$no_moda[1]$text$value),as.double(t$val_list_modal[2]$frequency[1]$text$value))
				if (is.null(indivNMA))
					indivNMA<-cNom
				else
					indivNMA<-rbind(indivNMA,cNom)	
					}
			}
		}

		if (variables[j,"type"]=="N")
		{
			t<-indivMatXML[i]$ligmat[j]$valmat[1]
      if(!is.null(t$val_nomina[1]$text$value)){
			cNom<-c(i,j,t$val_nomina[1]$text$value)
			#print(paste("val_nomina",t$val_nomina[1]$text$value))
			if (is.null(indivNA) )
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
if (!is.null(indivNMA))
{
	row.names(indivNMA)<-NULL
	indivNM<-as.data.frame(indivNMA)
	names(indivNM)<-c("indiv","variable","value","frequency")
}
else
	indivNM<-NULL
resul<-list(individuals=indiv,variables=variables,detailsIC=detailsIC,detailsC=detailsC,detailsN=detailsN,detailsListNom=detailsListNom,detailsNM=detailsNM,detailsListNomModif=detailsListNomModif,indivIC=indivICA,indivC=indivCA,indivN=indivN,indivNM=indivNM)
class(resul)<-"symbolic"
resul
}

#sdt<-parse.SO("wine")