# Fonction pour convertir un data frame contenant les genotypes (paires d'alleles) de multiples sujets (lignes) et SNPs (paires de colonnes)
# en un autre data frame contenant, dans le cas allelic, les (nombres d'alleles mineurs)/2 (valeurs 0, 0.5 ou 1; une colonne par SNP);
# dans le cas dominant, la valeur 1 si le genotype est constitue d'au moins un allele mineur, 0 sinon;
# dans le cas recessif, la valeur 1 si le genotype est constitue de 2 alleles mineurs, 0 sinon.

# correspond au fichier alleles2sums_v2.R dans le dossier "programmes"

# par Jordie Croteau
# 2 juin 2011

alleles2sums=function(geno.table,MA.vec,snp.names,mode="allelic")
{
# geno.table : data frame des genotypes (2 colonnes par SNP) 
# MA.vec :     vecteur des alleles mineurs
# snp.names :  vecteur des noms de SNPs
# mode :       mode de transmission ("allelic" (defaut), "recessive" ou "dominant")

num.snps=ncol(geno.table)/2
geno.table[geno.table==0]=NA 

sums=rep(NA,nrow(geno.table))
for(j in 1:num.snps){
cols.tmp=geno.table[,c(-1,0)+2*j]
cols.tmp.u=unique(as.vector(as.matrix(cols.tmp)))
if(length(cols.tmp.u[!is.na(cols.tmp.u)])>2) stop(paste("SNP",snp.names[j],"has more than two distinct alleles"))
sums=data.frame(sums,apply(cols.tmp,1,function(alleles,minor) sum(alleles==minor),minor=MA.vec[j]))
}
sums=as.matrix(sums[,-1])
colnames(sums)=snp.names

if(mode=="allelic") sums=sums/2

if(mode=="recessive"){
sums[sums==1]=0
sums[sums==2]=1
}

if(mode=="dominant") sums[sums==2]=1

sums
}



