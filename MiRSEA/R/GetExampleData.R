

GetExampleData<-function(exampleData){

if(!exists("envData")) {
envData<-initializeMiRSEA()
}
if (exampleData=="dataset")
{
dataset<- get("example.GCT",envir=envData)

return(dataset)
}



if(exampleData=="class.labels")
{

class.labels<- get("example.CLS",envir=envData)

return(class.labels)
}

if (exampleData=="p_value")
{
p_value<- get("p",envir=envData)

return(p_value)
}

if (exampleData=="p2miR")
{
p2miR<- get("p2miR",envir=envData)

return(p2miR)
}
if(exampleData=="miRList")
{
 miRList<-get("miRList",envir=envData)
 MiRList2<-list()
 
line1<-noquote(unlist(strsplit(miRList[1], "\t")))
line1<-line1[line1!=""]
N<-as.numeric(line1)
line2<-noquote(unlist(strsplit(miRList[2], "\t")))
Obs.RES<-as.numeric(line2[line2!=""])
line3<-noquote(unlist(strsplit(miRList[3], "\t")))
obs.s2n<-as.numeric(line3[line3!=""])
line4<-noquote(unlist(strsplit(miRList[4], "\t")))
Obs.ES<-as.numeric(line4[line4!=""])
line5<-noquote(unlist(strsplit(miRList[5], "\t")))
size.M<-as.numeric(line5[line5!=""])
line6<-noquote(unlist(strsplit(miRList[6], "\t")))
Obs.arg.ES<-as.numeric(line6[line6!=""])
line7<-noquote(unlist(strsplit(miRList[7], "\t")))
Obs.indicator<-as.numeric(line7[line7!=""])
line8<-noquote(unlist(strsplit(miRList[8], "\t")))
phen1<-as.character(line8[line8!=""])
line9<-noquote(unlist(strsplit(miRList[9], "\t")))
phen2<-as.character(line9[line9!=""])
line10<-noquote(unlist(strsplit(miRList[10], "\t")))
obs.index<-as.numeric(line10[line10!=""])
line11<-noquote(unlist(strsplit(miRList[11], "\t")))
obs.miR.labels<-as.character(line11[line11!=""])
line12<-noquote(unlist(strsplit(miRList[12], "\t")))
MsNAME<-as.character(line12[line12!=""])
MiRList2<-list(N,t(Obs.RES),t(obs.s2n),Obs.ES,size.M,Obs.arg.ES,t(Obs.indicator),phen1,phen2,t(obs.index),t(obs.miR.labels),MsNAME)

return(MiRList2)
}
}