genotypePlot <- function(snp,gene,eqtl,geneAnnot=NULL, ylab=NULL, xlab=NULL, mainlab=FALSE){
   # First bring the geno objects into the same order as the expression values are...
   rowGex <- rownames(eqtl$gex)
   takeThese <- is.element(as.character(eqtl$geno$fam[,1]),rowGex)
   eqtl$geno$genotypes <- eqtl$geno$genotypes[takeThese,]
   eqtl$geno$genotypes <- eqtl$geno$genotypes[order(rowGex),]
   snpCol <- which((eqtl$geno$map[,2]==snp)==TRUE)
   snpValues <- as.numeric(as.vector(eqtl$geno$genotypes[,snpCol]))
   nGroups <- unique(snpValues)
   grExpr <- list()
   for(i in 1:length(nGroups)){
     temp <- rownames(eqtl$geno$genotypes)[snpValues==i]
     grExpr[[i]] <- eqtl$gex[is.element(rowGex,temp)]
   }
   
   temp <- eqtl$eqtl[names(eqtl$eqtl)==gene]
   temp2 <- temp[2]$TestedSNP[,2]==snp
   
   snpP <- temp[3]$p.values[temp2]
   
   if(is.null(geneAnnot)) geneAnnot <- gene
   ifelse(mainlab==TRUE, mainTitle <- paste(snp," and ",geneAnnot,". (P-value:",snpP,")",sep=""), mainTitle <- "")
   
   boxplot(grExpr[[1]],xlim=c(0.5,length(grExpr)+0.5),ylim=c(min(as.vector(eqtl$gex)),max(as.vector(eqtl$gex))),main=mainTitle, ylab=ylab, xlab=xlab)
   for(i in 2:length(grExpr)){
     boxplot(grExpr[[i]],at=i,add=TRUE)
   }
   
   # Fill the x-axis
   tempRow <- eqtl$geno$map[eqtl$geno$map[,2]==snp,]
   al1 <- as.character(tempRow[5])
   al2 <- as.character(tempRow[6])
   axis(1,at=c(1,2,3),labels=c(paste(al1,"/",al1,sep=""),paste(al1,"/",al2,sep=""),paste(al2,"/",al2,sep="")))
}
