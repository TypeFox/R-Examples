filter.firstlayer <-
function(data,organism,ginumber,spaanscore,subcelllocal,tmhelices,Hrefhits)
{xa<- read.table(data,header=F,sep="\t");
	if(missing(ginumber)==TRUE){ginumber="ALL"}
	if(missing(spaanscore)==TRUE){spaanscore=0.6}
	if(missing(subcelllocal)==TRUE){subcelllocal= c("Extracellular","Cellwall")}
	if(missing(tmhelices)==TRUE){tmhelices=2}
	if(missing(Hrefhits)==TRUE){Hrefhits="No Hits found"}
	x=which(as.character(xa[,3])==organism)
	res1=xa[x,]
	if(ginumber[1] != "ALL")
	 {
	   res2=res1[which(as.numeric(res1[,1])==ginumber),]
	   res3=res2[which(as.numeric(res2[,5])>spaanscore),]
	   res4=res3[union(which(as.character(res3[,13])==subcelllocal[1]),which(as.character(res3[,13])==subcelllocal[2])),]
	   res5=res4[which(as.numeric(res4[,14])<tmhelices),]
	   res6=res5[which(as.character(res5[,17])==Hrefhits),]	
	   colnames(res6)=c("ginumber","annotation","organism","length","SPAAN_score","paralogs","orthologs_inter_sps","orthologs_intra_sps","COG","signalp","is_signalP","psortbscore","subcellular_localization_from_psortb","No_of_TMhelix","topology_of_TMhelix","betawrap","Hrefhits","cddhits")	   
	   write.table(res6,file="filtered_firstlayer.txt",sep="\t",row.names=F)
	  }
	else 
	 {res2=res1
	 	res3=res2[which(as.numeric(res2[,5])>spaanscore),]
	   res4=res3[union(which(as.character(res3[,13])==subcelllocal[1]),which(as.character(res3[,13])==subcelllocal[2])),]
	   res5=res4[which(as.numeric(res4[,14])<tmhelices),]
	   res6=res5[which(as.character(res5[,17])==Hrefhits),]
 colnames(res6)=c("ginumber","annotation","organism","length","SPAAN_score","paralogs","orthologs_inter_sps","orthologs_intra_sps","COG","signalp","is_signalP","psortbscore","subcellular_localization_from_psortb","No_of_TMhelix","topology_of_TMhelix","betawrap","Hrefhits","cddhits")
	   write.table(res6,file="filtered_firstlayer.txt",sep="\t",row.names=F)
	   }
  return(res6)
	 }
