save.SO<-function(sdt,file){
doc=newXMLDoc()
#top=newXMLNode("assofile",doc=doc)
top=newXMLNode("assofile",doc=doc,attrs=c("xmlns:xsi"="http://www.w3.org/2000/10/XMLSchema-instance","xsi:noNamespaceSchemaLocation"="asso.xsd"),
namespaceDefinitions = c("xsi" = "http://www.w3.org/2000/10/XMLSchema-instance"))
newXMLNode("contains",attrs=c(DIST_MATRIX="NO",FILES="YES",HEADER="YES",HIERARCHIE="NO",INDIVIDUALS="YES",MATCH_MATRIX="NO",RECTANGLE_MATRIX="NO",RULES="NO",VARIABLES="YES"),parent=top)
filed=newXMLNode("filed",parent=top)
                                                                                                                                                                                                                                             newXMLNode("Procedure","SDAEdit",parent=filed)
newXMLNode("version","0.99",parent=filed)
newXMLNode("date_creat",date(),parent=filed)
#indivNo<-0
indivNo<-nrow(sdt$individuals)
varNo<-nrow(sdt$variables)
ICNo<-sum(sdt$variables[,"type"]=="IC")
NNo<-sum(sdt$variables[,"type"]=="N")
NMNo<-sum(sdt$variables[,"type"]=="NM")
MNNo<-sum(sdt$variables[,"type"]=="MN")
newXMLNode("header",parent=top,attrs=c("indiv_nb"=indivNo,
nb_hierarchies="0",nb_na="0",nb_nu="0",nb_var_cont="0",nb_var_cont_symb=ICNo,nb_var_nom=NNo,nb_var_nom_mod=NMNo,nb_var_nom_symb=MNNo,nb_var_text="0",rules_nb="0",sub_title="auto",title="auto",var_nb=varNo))
individus<-newXMLNode("individus",parent=top)
for(i in 1:nrow(sdt$individuals)){
indiv<-newXMLNode("stindiv",parent=individus)
newXMLNode("num",sdt$individuals[i,1],parent=indiv)
newXMLNode("name",sdt$individuals[i,2],parent=indiv)
newXMLNode("label",sdt$individuals[i,3],parent=indiv)
}
variables<-newXMLNode("variables",parent=top)
for(i in 1:nrow(sdt$variables)){
var<-newXMLNode("stvar",parent=variables)
ident<-newXMLNode("ident",parent=var)
newXMLNode("num",sdt$variables[i,1],parent=ident)
newXMLNode("name",sdt$variables[i,2],parent=ident)
newXMLNode("label",sdt$variables[i,3],parent=ident)
if(sdt$variables[i,"type"]=="IC"){
varDesc<-newXMLNode("inter-cont",parent=var)
newXMLNode("continue-desc", attrs=c(max=max(sdt$indivIC[,i,2]),min=min(sdt$indivIC[,i,1]),nbna="0",nbnu="0"),parent=varDesc)
}
if(sdt$variables[i,"type"]=="N"){
nominal<-newXMLNode("nominal",parent=var)
nominalDesc<-newXMLNode("nominal-desc",parent=nominal,attrs=c(nbmoda=sum(as.matrix(sdt$detailsListNom)[,"details_no"]==as.matrix(sdt$variables)[i,"details"]), nbna="0" ,nbnu="0"))
for(j in 1:nrow(sdt$detailsListNom)){
  if(as.matrix(sdt$detailsListNom)[j,"details_no"]==as.matrix(sdt$variables)[i,"details"]){
    listNom=newXMLNode("list-nom",parent=multNominal)
    newXMLNode("num",sdt$detailsListNom[j,2],parent=listNom)
    newXMLNode("name",sdt$detailsListNom[j,3],parent=listNom)
    newXMLNode("label",sdt$detailsListNom[j,4],parent=listNom)
  }
}

}
if(sdt$variables[i,"type"]=="MN"){
multNominal<-newXMLNode("mult_nominal",parent=var)
nominalDesc<-newXMLNode("nominal-desc",parent=multNominal,attrs=c(nbmoda=sum(as.matrix(sdt$detailsListNom)[,"details_no"]==as.matrix(sdt$variables)[i,"details"]),nbna="0", nbnu="0"))
for(j in 1:nrow(sdt$detailsListNom)){
  if(as.matrix(sdt$detailsListNom)[j,"details_no"]==as.matrix(sdt$variables)[i,"details"]){
    listNom=newXMLNode("list-nom",parent=nominalDesc)
    newXMLNode("num",sdt$detailsListNom[j,2],parent=listNom)
    newXMLNode("name",sdt$detailsListNom[j,3],parent=listNom)
    newXMLNode("label",sdt$detailsListNom[j,4],parent=listNom)
  }
}

}
if(sdt$variables[i,"type"]=="NM"){
multNominal<-newXMLNode("mult_nominal_Modif",parent=var)
newXMLNode("<type_modif","proba",parent=var)
nominalDesc<-newXMLNode("nominal-desc",parent=multNominal,attrs=c(nbmoda=sum(as.matrix(sdt$detailsListNomModif)[,"details_no"]==as.matrix(sdt$variables)[i,"details"]) ,nbna="0" ,nbnu="0"))
for(j in 1:nrow(sdt$detailsListNom)){
  if(as.matrix(sdt$detailsListNomModif)[j,"details_no"]==as.matrix(sdt$variables)[i,"details"]){
    listNom=newXMLNode("list-nom",parent=multNominal)
    newXMLNode("num",sdt$detailsListNomModif[j,2],parent=listNom)
    newXMLNode("name",sdt$detailsListNomModif[j,3],parent=listNom)
    newXMLNode("label",sdt$detailsListNomModif[j,4],parent=listNom)
  }
}

}
}
if(!is.null(sdt$indivN))indivN<-as.matrix(sdt$indivN)
if(!is.null(sdt$indivNM))indivNM<-as.matrix(sdt$indivNM)
#print("po konwersji")
indivMat<-newXMLNode("indiv_mat",parent=top)
for(i in 1:nrow(sdt$individuals)){
ligmat<-newXMLNode("ligmat",parent=indivMat)
for(j in 1:nrow(sdt$variables)){
valMat<-newXMLNode("valmat",parent=ligmat)
#print(paste(i,j))
if(sdt$variables[j,"type"]=="IC"){
  valInterv<-newXMLNode("val_interv",parent=valMat)
  newXMLNode("pmin",sdt$indivIC[i,j,1],parent=valInterv)
  newXMLNode("pmax",sdt$indivIC[i,j,2],parent=valInterv)
}
if(sdt$variables[j,"type"]=="N"){
  for(k in 1:nrow(indivN)){
    if((indivN[k,"indiv"]==i) && (indivN[k,"variable"]==j)){
      newXMLNode("val_nominal",indivN[k,"value"],parent=valMat)
    }
  }
}
if(sdt$variables[j,"type"]=="MN"){
  for(k in 1:nrow(indivN)){
    if((indivN[k,"indiv"]==i) && (indivN[k,"variable"]==j)){
      newXMLNode("val_modal",indivN[k,"value"],parent=valMat)
    }
  }
}
if(sdt$variables[j,"type"]=="NM"){
  for(k in 1:nrow(indivNM)){
    if((indivNM[k,"indiv"]==i) && (indivNM[k,"variable"]==j)){
      vlm<-newXMLNode("val_list_modal",parent=valMat)
      newXMLNode("no_moda",indivN[k,"value"],parent=vlm)
      newXMLNode("frequency",indivN[k,"frequency"],parent=vlm)
    }
  }
}
}
}

saveXML(doc,file=file,prefix='<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n<?xml-stylesheet type="text/xsl" href="asso2.1.xsl"?>\n')
}

