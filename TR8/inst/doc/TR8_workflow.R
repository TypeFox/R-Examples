### R code from vignette source 'TR8_workflow.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: dryad (eval = FALSE)
###################################################
## ## the readxl package is needed
## ## library(readxl)
## ## store the  url of the dryad package
## url<-"http://datadryad.org/bitstream/handle/
##     10255/dryad.65646/MEE-13-11-651R2_data.xlsx?sequence=1"
## ## choose the extension for the temp file where 
## ## data will be stored
## tmp = tempfile(fileext = ".xlsx")
## ## download the data
## download.file(url = url, destfile = tmp)
## 
## ## we first read the "metadata" sheet from the xlsx file
## ## (the row containing the species names start from 
## ## row 13
## metadata<-read_excel(path=tmp,sheet="metadata",skip=12,col_names=F)
## ## lets rename the column of this dataset 
## names(metadata)<-c("Col1","Col2")
## 
## ## then read the vegetation data
## veg_data <-readWorksheetFromFile(file = tmp, sheet = "data.txt")
## ## only the columns from 11 to 123 contains the species data
## veg_data<-veg_data[,11:123]
##        ## round veg_data numbers to the second digit
## veg_data<-round(veg_data,digits = 2)
## ## read the dataset with the environmental variables
## env_data<-read_excel(path = tmp, sheet = "data.txt")
## ## and select only the column from 1 to 4 which contain
## ## the data of interest
## env_data<-env_data[,1:4]


###################################################
### code chunk number 2: taxize (eval = FALSE)
###################################################
## library(taxize)
## check_names<-tnrs(metadata$Col2,source="iPlant_TNRS")


###################################################
### code chunk number 3: discarded (eval = FALSE)
###################################################
## setdiff(metadata$Col2,check_names$submittedname)


###################################################
### code chunk number 4: taxize2 (eval = FALSE)
###################################################
## issues<-with(check_names,check_names[score!="1",])
## issues[,c("submittedname","score","acceptedname","authority")]


###################################################
### code chunk number 5: substitution (eval = FALSE)
###################################################
## library(plyr)
## ## we use the revalue function in the plyr package
## ## to fix all the above mentioned issues
## metadata$Col2<-revalue(metadata$Col2, 
##      c("Taraxacum officinale!!!!!"="Taraxacum officinale F.H. Wigg."))
## metadata$Col2<-revalue(metadata$Col2,
##      c("Polygonum mite (=Persicaria laxiflora)"="Persicaria mitis (Schrank) Assenov"))
## metadata$Col2<-revalue(metadata$Col2,
##      c("Fallopia convolvulus (L.) A. Löwe"="Fallopia convolvulus (L.) Á. Löve"))
## metadata$Col2<-revalue(metadata$Col2,
##      c("Setaria pumila (Poir.) Schult."="Setaria pumila (Poir.) Roem. & Schult."))
## metadata$Col2<-revalue(metadata$Col2,
##      c("Phleum pratense agg."="Phleum pratense L."))


###################################################
### code chunk number 6: taxize_2 (eval = FALSE)
###################################################
## check_names<-tnrs(metadata$Col2,source="iPlant_TNRS")
## issues<-with(check_names,check_names[score!="1",])
## issues[,c("submittedname","acceptedname","score")]


###################################################
### code chunk number 7: two (eval = FALSE)
###################################################
## final_dataframe<-merge(metadata,check_names,
##          by.x = "Col2",by.y="submittedname")


###################################################
### code chunk number 8: three (eval = FALSE)
###################################################
## final_dataframe<-final_dataframe[
##         !final_dataframe$Col2%in%issues$submittedname,]  


###################################################
### code chunk number 9: tr8_ex (eval = FALSE)
###################################################
## species_names<-final_dataframe$acceptedname
## my_traits<-c("h_max","le_area","leaf_mass","li_form_B","strategy")
## retrieved_traits<-tr8(species_list = species_names,download_list = my_traits)


###################################################
### code chunk number 10: cirsium (eval = FALSE)
###################################################
## ## we extract the data from the object returned by tr8()
## traits<-extract_traits(retrieved_traits)
## ## first I convert the column to character
## traits$h_max<-as.character(traits$h_max)
## traits$h_max[which(row.names(traits)=="Convolvulus arvensis")]<-"42.5"


###################################################
### code chunk number 11: convert (eval = FALSE)
###################################################
## traits$h_max<-as.numeric(traits$h_max)


###################################################
### code chunk number 12: leArea (eval = FALSE)
###################################################
## traits$le_area<-revalue(traits$le_area,
##       c("0.1-1"=0.55,
##         "1-10"=5.5,
##         "10-100"=55,
##         "100-1000"=550,
##         "1-10;0.1-1"=1,
##         "10-100;1-10"=10,
##         "100-1000;10-100"=100,
##         "10-100;100-1000"=100))
## ## and convert them to numeric
## traits$le_area<-as.numeric(as.character(traits$le_area))


###################################################
### code chunk number 13: liform (eval = FALSE)
###################################################
## 
## 
## traits$li_form_B<-revalue(traits$li_form_B,
##   c("C (Chamaephyte) - H (Hemicryptophyte)"="C - H",
##     "G (Geophyte)"="G",
##     "G (Geophyte) - H (Hemicryptophyte)"="G - H",
##     "H (Hemicryptophyte)"="H",
##     "H (Hemicryptophyte) - T (Therophyte)"="H - T",
##     "M (Macrophanerophyte)"="M",
##     "M (Macrophanerophyte) - N (Nanophanerophyte)"="M - N",
##     "T (Therophyte)"="T"))
## ## convert it to factor
## traits$li_form_B<-as.factor(traits$li_form_B)


###################################################
### code chunk number 14: strategy (eval = FALSE)
###################################################
## traits$strategy<-revalue(traits$strategy,c("c (competitors)"="c",
##  "cr (competitors/ruderals)"="cr",
##  "cs (competitors/stress-tolerators)"="cs",
##  "csr (competitors/stress-tolerators/ruderals)"="csr",
##  "r (ruderals)"="r"))
## traits$strategy<-as.factor(traits$strategy)


###################################################
### code chunk number 15: a (eval = FALSE)
###################################################
## row.names(traits)<-mapvalues(row.names(traits),
##    from=final_dataframe$acceptedname,to=final_dataframe$Col1)


###################################################
### code chunk number 16: b (eval = FALSE)
###################################################
## traits<-traits[complete.cases(traits),]


###################################################
### code chunk number 17: c (eval = FALSE)
###################################################
## vegetation<-veg_data[,names(veg_data)%in%row.names(traits)]


###################################################
### code chunk number 18: d (eval = FALSE)
###################################################
## library(ade4)
## coa<-dudi.coa(vegetation,scannf=F)


###################################################
### code chunk number 19: e (eval = FALSE)
###################################################
## hil.traits<-dudi.hillsmith(traits,row.w=coa$cw,scannf = FALSE)


###################################################
### code chunk number 20: f (eval = FALSE)
###################################################
## ##select which columns have at least one non-zero value
## selection<-colSums(vegetation)>0
## ## and now we choose only those columns
## vegetation<-vegetation[,selection]


###################################################
### code chunk number 21: g (eval = FALSE)
###################################################
## traits<-traits[row.names(traits)%in%names(vegetation),]


###################################################
### code chunk number 22: hh (eval = FALSE)
###################################################
## vegetation<- vegetation[,order(names(vegetation))]
## traits<-traits[order(row.names(traits)),]


###################################################
### code chunk number 23: h (eval = FALSE)
###################################################
## coa<-dudi.coa(vegetation,scannf=F)
## traits.hill<-dudi.hillsmith(traits,row.w=coa$cw,scannf = F)


###################################################
### code chunk number 24: i (eval = FALSE)
###################################################
## env.hill<-dudi.hillsmith(env_data,row.w=coa$lw,scannf = FALSE)


###################################################
### code chunk number 25: l (eval = FALSE)
###################################################
## env_data$Treat<-as.factor(env_data$Treat)


###################################################
### code chunk number 26: i (eval = FALSE)
###################################################
## env.hill<-dudi.hillsmith(env_data,row.w=coa$lw,scannf = FALSE)


###################################################
### code chunk number 27: l (eval = FALSE)
###################################################
## rlq_tr8<-rlq(env.hill,coa,traits.hill,scannf = F)


###################################################
### code chunk number 28: m (eval = FALSE)
###################################################
## plot(rlq_tr8)


###################################################
### code chunk number 29: m (eval = FALSE)
###################################################
## clust<-hclust(dist(rlq_tr8$lQ),method="ward.D2")
## plot(clust,sub="Ward minimum variance clustering",xlab="TR8 tutorial")


###################################################
### code chunk number 30: o (eval = FALSE)
###################################################
## rect.hclust(clust,k=6)


###################################################
### code chunk number 31: p (eval = FALSE)
###################################################
## cuts<-cutree(clust,6)


###################################################
### code chunk number 32: q (eval = FALSE)
###################################################
## s.class(rlq_tr8$lQ,as.factor(cuts),col=1:6)
## s.arrow(rlq_tr8$c1,add.plot = TRUE)


###################################################
### code chunk number 33: aa (eval = FALSE)
###################################################
## par(mfrow=c(3,2))
## plot(traits$h_max~as.factor(cuts),main="Maxim height",
##      ylab="max height",border = 1:6,xlab="Group number")
## plot(traits$le_area~as.factor(cuts),main="Leaf area",
##      ylab="leaf area",border = 1:6,xlab="Group number")
## plot(traits$leaf_mass~as.factor(cuts),main="Leaf mass",
##      ylab="leaf mass",border = 1:6,xlab="Group number")
## plot(table(cuts,traits$strategy),main="CSR strategy",
##      ylab="strategy",border = 1:6,xlab="Group number")
## plot(table(cuts,traits$li_form_B),main="Life form",
##      ylab="life form",border = 1:6,xlab="Group number")
## par(mfrow=c(1,1))


