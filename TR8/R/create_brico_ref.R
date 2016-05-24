## run this script if you want to recreate the lookup table which contains
## species names withouth authors for briconet dataset


## url<-"http://bricol.net/downloads/data/PLANTSdatabase/09-02-02PLANTSdata.csv"
## brico<-read.table(url,fileEncoding="iso-8859-1")

## ## remove commas in species names since they conflict with TNRS

## brico$Scientific.Name<-as.character(brico$Scientific.Name)
## brico$Scientific.Name<-gsub(","," ",brico$Scientific.Name)

## go<-tnrs(brico$Scientific.Name ,source = "iPlant_TNRS")
## go$score<-as.numeric(go$score)
## dai<-go[,c("submittedname","acceptedname","score","matchedname","authority")]
## dai<-dai[dai$score>=.6,]
## names(dai)<-revalue(names(dai),c("submittedname"="Scientific.Name"))
## ref_brico<-merge(dai,brico)

## ref_brico<-ref_brico[,c("Scientific.Name","acceptedname","score","matchedname","authority","Symbol","Synonym.Symbol")]

## save(ref_brico,file="~/Work/src/tr8_github/data/ref_brico.Rda")
