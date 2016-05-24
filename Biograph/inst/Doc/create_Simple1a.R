# test data
rm(list=ls())
library (Biograph)
# Example 1: three hypothetical individuals 
id <- c(1,2,3)
born <- c("1986-04-05","1986-08-08","1986-11-28")
interview <- rep("2019-05-09",3)
sex <- factor(c("F","M","F"))
educ <- factor(c("High","Medium","Medium"))

namstates <- c("H","A","C","M")
A <- c("2004-08-20","2011-09-15","2006-08-10")
C <- c("2011-12-1",NA,NA)
M <- c(NA,NA,"2012-03-16")
d <- data.frame(A=A,C=C,M=M,stringsAsFactors =FALSE)
nsample <- nrow(d)
 # d= character object
 dd<- apply(d,1,function(x) y=as.Date(x))
 dd <- data.frame(t(dd))  #  dd is numeric
dimnames(dd) <- dimnames(d)
f <- Sequences.ind.0(dd,namstates,absorb=NULL)
dates <- data.frame (f$d)
for (i in 1:3)
 {dates[,i] <- as.Date(dates[,i],origin="1970-01-01") 
 }
path <- as.character(f$path)

bio  <- data.frame (ID=id,born=born,start=born,end=interview,sex=sex,educ=educ,path=as.character(path),dates[,1:(max(nchar(path))-1)],stringsAsFactors=FALSE)
namtrans <- paste("Tr",1:ncol(f$d),sep="")
colnames(bio)[8:9] <- namtrans[1:2]

attr(bio,"format.date") <- "%Y-%m-%d"
attr(bio,"format.born") <- "%Y-%m-%d"
attr(bio,"param") <- Parameters (bio)

# Path to folder where bio should be stored
zzz <- "/Users/frans/Documents/R/0 0 MAC/Package/Biograph.TEST/Chapters/AnnexA/Simple 1"
setwd(zzz)
save (bio,file="simple1a.RData")

