# test data
rm(list=ls())
library (Biograph)
# Define state space
namstates <- c("H","A","C","M")
# Specify subject identification numbers
id <- c(1,2,3)
# Dates of birth
born <- c(1036,1040,1043)
# Observation window
start <- born   
interview <- rep(1433,3)
# Covariates
sex <- factor(c("F","M","F"))
educ <- factor(c("High","Medium","Medium"))

# Create a data frame of dates of transition
A <- c(1256,1341,1280)
C <- c(1344,NA,NA)
M <- c(NA,NA,1347)
d <- data.frame(A=A,C=C,M=M,stringsAsFactors =FALSE)
# nsample = sample size
nsample <- nrow(d)
# Sort and rearrange the dates and obtain state sequence (path)
f <- Sequences.ind.0(d,namstates)
dates <- data.frame (f$d)
colnames(dates) <- paste ("Tr",1:(length(namstates)-1),sep="")
# State sequence (path)
path <- as.character(f$path)
# Number of states (ns)

# Create Biograph object
bio  <- data.frame (ID=id,born=born,start=start,end=interview,sex=sex,educ=educ,path=as.character(path),dates[,1:(max(nchar(path))-1)],stringsAsFactors=FALSE)

attr(bio,"format.date") <- "CMC"
attr(bio,"format.born") <- "CMC"
param <- Parameters (bio)
attr(bio,"param") <- param
# ---------  Biograph object completed  --------------

