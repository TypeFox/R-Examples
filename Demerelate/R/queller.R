queller <- function(pop1, pop2, allele.column, row.position, col.position, matrix.share, ref.pop)
	{ 

  
  pm1 <- table(rbind(ref.pop[,allele.column],ref.pop[,allele.column+1]))
  pm <- pm1/sum(pm1)
  

  # Calculations after Oliehoek et al 2006 (15) + (16) 
  # Share.RE for colnames/pop2
  pxma <- as.numeric(pm[which(names(pm)==pop1[which(row.names(matrix.share)[row.position]==pop1[,1]),allele.column])])
	pymb <- as.numeric(pm[which(names(pm)==pop1[which(row.names(matrix.share)[row.position]==pop1[,1]),allele.column+1])])
  het <- as.numeric(table(as.numeric(pop1[which(row.names(matrix.share)[row.position]==pop1[,1]),allele.column:(allele.column+1)]))[1])-1
  
  share.RE <- (sum(
  as.numeric(pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),allele.column]==pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),allele.column]),
  as.numeric(pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),allele.column]==pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),allele.column+1]),
  as.numeric(pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),allele.column+1]==pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),allele.column]),
	as.numeric(pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),allele.column+1]==pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),allele.column+1])
	)*0.5-pxma-pymb)/(1+het-pxma-pymb)
  
  
  # Share.RAT for pop2
  # het == 1 if Alleles identical otherwise het == 0
  pxma <- as.numeric(pm[which(names(pm)==pop2[which(colnames(matrix.share)[col.position]==pop2[,1]),allele.column])])
  pymb <- as.numeric(pm[which(names(pm)==pop2[which(colnames(matrix.share)[col.position]==pop2[,1]),allele.column+1])])
  het <- as.numeric(table(as.numeric(pop2[which(colnames(matrix.share)[col.position]==pop2[,1]),allele.column:(allele.column+1)]))[1])-1
    
  share.RAT <- (sum(
  as.numeric(pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),allele.column]==pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),allele.column]),
  as.numeric(pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),allele.column]==pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),allele.column+1]),
  as.numeric(pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),allele.column+1]==pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),allele.column]),
	as.numeric(pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),allele.column+1]==pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),allele.column+1])
	)*0.5-pxma-pymb)/(1+het-pxma-pymb)

  
  share <- (share.RE+share.RAT)/2		
  if (length(share)>0){return(share)}else{return(NA)}
  
	}
