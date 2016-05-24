`is.snp` <-
function(data,sep){
   j<-1
   if (!is.data.frame(data)) stop("Data is not a data.frame object. \n")
   if (any(apply(data,2,is.factor))) stop("Some data columns are factors. \n")
   if (ncol(data)==1) stop("There is only one SNP, is not possible to perform a haplotype analysis. \n")
   
   col.rm<-NULL
   for (j in 1:ncol(data))

        {

                if (sum(is.na(data[,j]))==nrow(data)) stop("At least one SNP values are all missings. \n")

                geno<-genotype(data[,j],sep=sep)
                a<-allele(geno) 
                if (length(attr(a,"allele.names"))==1){
                    warning("At least one SNP is homozygous for everyone.\n")
                    #col.rm<-c(col.rm,j)
                }                
                alleles<- sort(attr(a,"allele.names"))
                if (length(alleles)>2) stop("SNP data is not biallelic.\n")
                if ( sum(geno[!is.na(geno)]%in%c(paste(alleles[1],alleles[1],sep=sep),paste(alleles[1],alleles[2],sep=sep),
                   paste(alleles[2],alleles[1],sep=sep),paste(alleles[2],alleles[2],sep=sep)))!=nrow(data)-sum(is.na(geno)))
                        stop("SNP data is not correct.\n")

          } 
    #data<-data[,-col.rm]

}

