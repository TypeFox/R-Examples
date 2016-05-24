allermatch <-
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
	   colnames(res3)=c("ginumber","species","prediction","hit no","Db","Allermatch id","best hit iden","No_hits_ident_g35","Perc_hits_gt35","prec_ident","Seq_len_fasta_aligned","external_link","link_db","Genus_name","Spc_name")
	   write.table(res3,file="filtered_allermatch.txt",sep="\t",row.names=F)
	  }
	else 
	 {res2=res1
	 res3=res2[which(as.character(res2[,3])==prediction),]
	    colnames(res3)=c("ginumber","species","prediction","hit no","Db","Allermatch id","best hit iden","No_hits_ident_g35","Perc_hits_gt35","prec_ident","Seq_len_fasta_aligned","external_link","link_db","Genus_name","Spc_name")
	   write.table(res3,file="filtered_allermatch.txt",sep="\t",row.names=F)
	   }
  return(res3)
  #write.table(res6,file="filtered_firstlayer.txt")
	 }
