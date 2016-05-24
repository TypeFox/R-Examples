LBI <-
function(infile,name.M="M.norm",ind.array=1:2,graph=TRUE,graphout="FigM1M2") {

ind <- which(names(infile) %in% paste(name.M,ind.array,sep=""))
cond <- as.data.frame(colSums(apply(infile[,ind],1,is.na)))
nbgene1 <- table(cond)[unique(cond)==1]
condb <- as.data.frame(colSums(apply(infile[,ind],1,is.na)))
nbgene2 <- table(condb)[unique(condb)==2]
cat("Number of genes with only one observation :",nbgene1,"\n")
cat("Number of genes wihtout any observation :",max(nbgene2,0),"\n")

VarParGene2 <- rowVars(infile[,ind],na.rm=TRUE)
VarParGene1 <- rowVars(data.frame(infile[,ind[1]],-infile[,ind[2]]),na.rm=TRUE)

Var2 <- mean(VarParGene2,na.rm=TRUE)
Var1 <- mean(VarParGene1,na.rm=TRUE)
cat("Variance per gene (Mean Log-ratio) : ", Var1,"\n")
cat("Variance per gene (Dye bias) : ", Var2,"\n")
cat("LBI : ", Var1/Var2,"\n")

if (graph) {
png(file=paste(graphout,".png",sep=""))  
plot(infile[,ind[1]],infile[,ind[2]],xlim=c(-3,3),ylim=c(-3,3),xlab="log2(R1/G1)",ylab="log2(R2/G2)")
dev.off()
}
# (c) 2010 Institut National de la Recherche Agronomique
 
}

