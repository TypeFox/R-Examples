summary.simfam <- function(object, ...){


design<-match(attr(object,"design"), c("pop","pop+","cli","cli+","twostage"))
design.name <- c("Population-based study with affected probands",
"Population-based study with affected and mutation carrier probands",
"Clinic-based study with affected probands",
"Clinic-based study with affected and mutation carrier probands",
"Two-stage design")[design]

cat("Study design:     ", design.name,"\n")

k <- ifelse(design==5,2,1)
for(i in 1:k){
if(design==5){
	 if(i==1){ 
	 	data <- object[object$weight==1,]
		cat(" ## High risk families", "\n")
	 }
	 else{
	 	data <- object[object$weight!=1,]
		cat(" ## Low risk families", "\n")
	 }
}
else data<- object	 

numfam <- length(data$famID[data$proband==1])
avgnumaffec <- round(sum(data$status)/numfam,2)
avgnumcarr <- round(sum(data$mgene,na.rm=T)/numfam,2)
avgfamsize<-round(mean(data$fsize[data$proband==1]),2)
avgageonset<-round(mean(data$ageonset[data$status==1]),2)
wt <- round(data$weight[1],2)
#ans<-data.frame(numfam, totalnumaffec,avgnumaffec, avgfamsize,avgage,avgageonset)
ans<-list(design=design.name, numfam=numfam, avgnumaffec=avgnumaffec, avgnumcarr=avgnumcarr, avgfamsize=avgfamsize,avgageonset=avgageonset, weight=wt)

nn<-c(
"Number of families:                    ", 
"Average number of affected per family: ", 
"Average number of carriers per family: ", 
"Average family size:                   ",
"Average age of onset for affected:     ",
"Sampling weight:                       ")

for(j in 2:7)  cat(nn[j-1], ans[[j]], "\n")
cat("\n")
}

}
