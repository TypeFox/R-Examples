library(R.utils)
## Data can be loaded at this place
## "http://cbbp.thep.lu.se/~markus/software/BAFsegmentation/"

data = read.table("CRL2324_dilutionSeries_TableExport", header=TRUE, sep="\t")

## Keep only interesting columns
## Use genomic DNA of breast carcinoma cells mixed with DNA from lymphoblastoid
## cells at known proportions (100, 79, 50, 34, 14).
## Genotype use are those of lymphoblastoid cells.

dat <- data[,c("Chr", "Position", "CRL2325.GType", "CRL2324.Log_R_Ratio",
               "CRL2324.B_Allele_Freq", "CRL2324_79pc_Tum.Log_R_Ratio",
               "CRL2324_79pc_Tum.B_Allele_Freq", "CRL2324_50pc_Tum.Log_R_Ratio",
               "CRL2324_50pc_Tum.B_Allele_Freq", "CRL2324_34pc_Tum.Log_R_Ratio",
               "CRL2324_34pc_Tum.B_Allele_Freq", "CRL2324_14pc_Tum.Log_R_Ratio",
               "CRL2324_14pc_Tum.B_Allele_Freq")]

## Rename columns
names(dat) <-  c("chromosome", "Position", "genotype","logR_100", "baf_100",
                 "logR_79", "baf_79","logR_50", "baf_50","logR_34", "baf_34",
                 "logR_14", "baf_14")

## Define genotype
dat$muN <- (as.numeric(dat$genotype)-1)/2
dat$muN[which(dat$muN==1.5)] <- NA

saveObject(dat,"CRL2324_dilutionSeries.xdr")
