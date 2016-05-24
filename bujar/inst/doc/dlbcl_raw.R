#Data source
#http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?token=rhojvaiwkcsaihq&acc=GSE10846
#Series Matrix
#Citation(s) 	
#    * Lenz G, Wright G, Dave SS, Xiao W et al. Stromal gene signatures in large-B-cell lymphomas. N Engl J Med 2008 Nov 27;359(22):2313-23. PMID: 19038878
### GSE10846_series_matrix-1.txt: CHOP and the first part of CHOP data
### GSE10846_series_matrix-2.txt: the second part of R-CHOP data
# Note: 1. It seems the original GSE10846_series_matrix.txt was split to two txt files as below.
#       2. The updated file from the above webpage can be different, for instance, number of variables
#       3. The last column 54722 below is empty with variable name !series_matrix_table_end, thus omitted

lymph1 <- t(read.delim("GSE10846_series_matrix-1.txt",sep="",header=FALSE,skip=33))[,-54722]
colnames(lymph1) <- lymph1[1,]
chop <- lymph1[2:182,]  ### CHOP patients
tmp <- substr(chop[,9], 34, 38)=="DEAD"
chop[tmp,9] <- 1
chop[!tmp,9] <- 0
colnames(chop)[9:10] <- c("status","survtime")
chop[,10] <- substr(chop[,10], 33, 36)
chop <- chop[,c(10,9,47:54721)]
storage.mode(chop) <- "numeric"

lymph2 <- t(read.delim("GSE10846_series_matrix-2.txt",sep="",header=FALSE,skip=33))[,-54722]
colnames(lymph2) <- lymph2[1,]
rchop <- rbind(lymph1, lymph2[-1,])[183:(183+232),] ### R-CHOP patients
tmp <- substr(rchop[,9], 34, 38)=="DEAD"
rchop[tmp,9] <- 1
rchop[!tmp,9] <- 0
colnames(rchop)[9:10] <- c("status","survtime")
rchop[,10] <- substr(rchop[,10], 33, 36)
rchop <- rchop[,c(10,9,47:54721)]
storage.mode(rchop) <- "numeric"
#begin with column 78, there are 54675 rows as genes

chop <- as.data.frame(chop)
rchop <- as.data.frame(rchop)
cat("\ncensoring rate in the full data is",sum(chop$status==0)/nrow(chop),"\n")
### unsupervised gene screening, cf: Additive risk survival model with microarray data
#Shuangge Ma corresponding author and Jian Huang, BMC Bioinformatics, 2007; 8:192
geoname <- read.delim("GPL570-39741.txt",sep="",header=TRUE,skip=16)[,1:15]
colnames(geoname) <- c(colnames(geoname)[-3],"")
### Lossos (2004), NEJM, Prediction of survival in diffuse large-B-cell lymphoma based on the expression of six genes
### 6 genes found (see below)
#tmp <- which(geoname$Gene.Symbol%in%c("LMO2","BCL6","FN1","CCND2","SCYA3","BCL2"))
#sixgene <- geoname$ID[tmp]
quar <- function(x) {
tmp <- quantile(x, c(0.25, 0.75))
if(tmp[2] - tmp[1] < tmp[1])
 return(FALSE)
else return(TRUE)
} ### inter quartile: difference between 3rd quartile and 1st one

lowvar <- function(x){
tmp <- quantile(x, 0.10)
if(var(x) < tmp)
return(FALSE)
else return(TRUE)
}

#res <- c("TRUE", "TRUE",apply(chop[, -(1:2)], 2, quar))
res <- c("TRUE", "TRUE",apply(chop[, -(1:2)], 2, lowvar))
table(res)
### there are 1871 - 2 genes passed the unsupervised screening
chop <- chop[,res=="TRUE"]

###selected genes to construct gene signatures in Lenz et al, 2008 NEJM, ref is in Lenz_data.txt, the first 5 lines
tmp <- read.delim(file="Lenz_data.txt",header=TRUE,sep="",skip=5)
n <- dim(tmp)[1]
res <- rep(NA,n)
for (i in 1:n){
 if(nchar(as.character(tmp[i,1])) == 7 && substr(tmp[i,1],1,1)==1)
 res[i] <- 1
}
tmp <- subset(tmp, res==1)

cat("How many the remaining genes match those gene signatures found in Lenz et al (2008)\n")
sum(colnames(chop)%in%tmp[,2])

