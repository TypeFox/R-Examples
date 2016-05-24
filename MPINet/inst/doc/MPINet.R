### R code from vignette source 'MPINet.Rnw'

###################################################
### code chunk number 1: MPINet.Rnw:49-50
###################################################
library(MPINet)


###################################################
### code chunk number 2: MPINet.Rnw:70-95
###################################################
#example 1
##########get example data
risk<-GetExampleData(dataset="prostate")

###########calculate the CGNB score
pss<-getPSS(risk ,plot=F)
CGNBscore<-pss[,"CGNB"]
names(CGNBscore)<-rownames(pss)
##########print the CGNB score of some metabolites to screen
head(CGNBscore)

#example 2
#get example data from file
risk<-read.table(paste(system.file(package="MPINet"),"/localdata/prostate.txt",sep=""),
header=F,sep="\t","\"")

####convert the data to a character vector
risk<-as.character(risk[[1]])

###########calculate the CGNB score
pss<-getPSS(risk ,plot=F)
CGNBscore<-pss[,"CGNB"]
names(CGNBscore)<-rownames(pss)
##########print the CGNB score of some metabolites to screen
head(CGNBscore)


###################################################
### code chunk number 3: MPINet.Rnw:117-157
###################################################
#example 1
#### get the metastatic prostate cancer interesting metabolite data set
risk<-GetExampleData(dataset="prostate")
#### integrate the non-equivalence of metabolites and the character of 
#### differential metabolites by the monotonic spline model 
pss<-getPSS(risk,plot=F)        

#identify dysregulated pathways
anncpdpre<-identifypathway(risk,pss,pathType="KEGG",method="MPINet",annlim=1,bglim=6)
#convert ann to data.frame
result<-printGraph(anncpdpre,pathType="KEGG",method="MPINet")
head(result)


#example 2
########get example data from file
risk<-read.table(paste(system.file(package="MPINet"),"/localdata/prostate.txt",sep=""),
header=F,sep="\t","\"") 

####convert the data to a character vector
risk<-as.character(risk[[1]])
pss<-getPSS(risk,plot=F)  

#identify dysregulated pathways
anncpdpre<-identifypathway(risk,pss,pathType="KEGG",method="MPINet",annlim=1,bglim=6)
#convert ann to data.frame
result<-printGraph(anncpdpre,pathType="KEGG",method="MPINet")
head(result)

#example 3
#### get the metastatic prostate cancer interesting metabolite data set
risk<-GetExampleData(dataset="prostate") 
pss<-getPSS(risk,plot=F)  

#identify dysregulated Reactome and KEGG pathways
anncpdpre<-identifypathway(risk,pss,pathType=c("KEGG","Reactome"),
                method="MPINet",annlim=1,bglim=6)
#convert ann to data.frame
result<-printGraph(anncpdpre,pathType=c("KEGG","Reactome"),method="MPINet")
head(result)


###################################################
### code chunk number 4: MPINet.Rnw:178-208
###################################################
#example 1
#### get the metastatic prostate cancer interesting metabolite data set
risk<-GetExampleData(dataset="prostate")
#### integrate the global non-equivalence of metabolites and the character of 
####differential metabolites by the monotonic spline model 
pss<-getPSS(risk,plot=F)        

#identify dysregulated pathways
anncpdpre<-identifypathway(risk,pss,pathType="KEGG",method="MPINet",annlim=1,bglim=6)
#convert ann to data.frame
result<-printGraph(anncpdpre,pathType="KEGG",method="MPINet")
#print part of the results to screen
head(result)
result1<-printGraph(anncpdpre,pathType="KEGG",method="MPINet",detail=TRUE)


#example 2
#### get the metastatic prostate cancer interesting metabolite data set
risk<-GetExampleData(dataset="prostate") 

pss<-getPSS(risk,plot=F)

#identify dysregulated pathways
anncpdpre<-identifypathway(risk,pss,pathType="Reactome",method="MPINet",annlim=1,bglim=6)
#convert ann to data.frame
result<-printGraph(anncpdpre,pathType="Reactome",method="MPINet")
#print part of the results to screen
head(result)

result1<-printGraph(anncpdpre,pathType="Reactome",method="MPINet",detail=TRUE)


###################################################
### code chunk number 5: sessionInfo
###################################################
sessionInfo()


