getFields<-function(seerHome="~/data/SEER"){
	sas=readLines(dir(pattern="*.sas",path=file.path(seerHome,"incidence"),full.names=TRUE)[1],warn=F)
	sas=sas[-c(1:5)]
	(sas=sas[nchar(sas)>15]) # remove trailing rows
	sas=gsub("[\\@\\$\\.]","",sas)  # an OR of three characters, @, $ or .
	sas=sub("char","",sas,ignore.case = TRUE)      # note one Char (not char) typo so ignore case
	sas=sub("/\\*","",sas); sas=sub("\\*/","",sas) # remove SAS comments
	sas=sub('^ +','',sas); 	sas=sub(' +$','',sas)  # remove front and back space
	sas=sub('/$','',sas)     # remove typo of SAS comment not closed at the end
	sas=gsub(' +', ' ',sas)  # take internal spaces of more than one to one
	sas=strsplit(sas," ")
	start=as.numeric(sapply(sas,function(x) x[[1]]))
	names=sapply(sas,function(x) x[[2]])
	width=as.numeric(sapply(sas,function(x) x[[3]]))
	desc=sapply(sas,function(x) paste(x[4:length(x)],collapse=" "))
	sas=SAS=data.frame(start,width,names,desc,stringsAsFactors=F)
	sas$names=tolower(gsub("_","",sas$names)); sas$names=gsub("v$","",sas$names)  # clean names a little
	# these are ones I disliked enough to change. Mostly I'm taking the SAS file names of fields.
 	# sas[which(sas$names=="datemo"),"names"]="modx" 
 	sas[which(sas$names=="mdxrecmp"),"names"]="modx" 
	sas$names[which(sas$names=="yeardx")]="yrdx" 
	# sas$names[which(sas$names=="dateyr")]="yrdx" 
	sas$names[which(sas$names=="icdoto9")]="ICD9" 
	sas$names[which(sas$names=="icdot10")]="ICD10" 
	# sas$names[which(sas$names=="icd5dig")]="COD" 
	sas$names[which(sas$names=="codpub")]="COD" 
	sas$names[which(sas$names=="age1rec")]="agerec" 
	sas$names[which(sas$names=="race1")]="race" 
	sas$names[which(sas$names=="pubcsnum")]="casenum" 
	sas$names[which(sas$names=="srvtimemon")]="surv" 
	# to spell out Collaborative Stage = CS, change a desc
  sas$desc[which(sas$names=="cssize")]="Collaborative Stage (CS) Tumor Size"

###  uncomment and run this when field names and/or positions change
#   setwd("inst/doc")
# 	tmp=cbind(SAS[,1:3],sas[,3:4]); colnames(tmp)[3:4]<-c("SAS","SEERaBomb");tmp
#   library(hwriter); hwrite(tmp, 'fieldNames.html', row.bgcolor='#ffdc98')#=file in doc directory
	sas
}

#pathPrep()
#fid=file("/Users/radivot/case/active/seer/SEERaBomb/R/getFields.R")
#ss=readLines(fid)
#iconv(ss, "latin1", "ASCII", "byte") 
