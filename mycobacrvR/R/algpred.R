algpred <-
function(data,organism,ginumber,prediction)
{xa<- read.table(data,header=F,sep="\t");
	if(missing(ginumber)==TRUE){ginumber="ALL"}
	if(missing(prediction)==TRUE){prediction="Non Allergen"}
	x=which(as.character(xa[,2])==organism)
	res1=xa[x,]
	if(ginumber[1] != "ALL")
	 {
	   res2=res1[which(as.numeric(res1[,1])==ginumber),]
	   res3=res2[which(as.character(res2[,3])==prediction),]
	   colnames(res3)=c("ginumber","species","overall pred","IgE pred","IgE epitope","Seq matched","position","Perc identity","MAST pred","Svm pred","Svn score","Svm thold","Svm ppv","Svm npv","Svm dipep pred","Svm dipep score","Svm dipep thold","Svm dipep ppv","Svm dipep npv","Blast pred","Hits_with_ARPs_database")
	   write.table(res3,file="filtered_algpred.txt",sep="\t",row.names=F)
	  }
	else 
	 { res2=res1
	   res3=res2[which(as.character(res2[,3])==prediction),]
	    colnames(res3)=c("ginumber","species","overall pred","IgE pred","IgE epitope","Seq matched","position","Perc identity","MAST pred","Svm pred","Svn score","Svm thold","Svm ppv","Svm npv","Svm dipep pred","Svm dipep score","Svm dipep thold","Svm dipep ppv","Svm dipep npv","Blast pred","Hits_with_ARPs_database")
	   write.table(res3,file="filtered_algpred.txt",sep="\t",row.names=F)
	   }
  return(res3)
	 }
