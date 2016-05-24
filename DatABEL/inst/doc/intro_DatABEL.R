### R code from vignette source 'intro_DatABEL.Rnw'

###################################################
### code chunk number 1: intro_DatABEL.Rnw:87-89
###################################################
unlink("*.fv?")
unlink("*.txt")


###################################################
### code chunk number 2: intro_DatABEL.Rnw:93-94
###################################################
library(DatABEL)


###################################################
### code chunk number 3: intro_DatABEL.Rnw:100-105
###################################################
matr <- matrix (c(1:12),ncol=3,nrow=4)
matr[3,2] <- NA
matr
dimnames(matr) <- list(paste("row",1:4,sep=""),paste("col",1:3,sep=""))
matr


###################################################
### code chunk number 4: intro_DatABEL.Rnw:116-119
###################################################
list.files(pattern="*.fv?")
dat1 <- as(matr,"databel")
list.files(pattern="*.fv?")


###################################################
### code chunk number 5: intro_DatABEL.Rnw:126-127
###################################################
dat1


###################################################
### code chunk number 6: intro_DatABEL.Rnw:136-138
###################################################
dat2 <- matrix2databel(matr, filename="matr",cachesizeMb=16, type="UNSIGNED_CHAR",readonly=FALSE)
dat2


###################################################
### code chunk number 7: intro_DatABEL.Rnw:142-143
###################################################
list.files(pattern="*.fv?")


###################################################
### code chunk number 8: intro_DatABEL.Rnw:151-153
###################################################
dat3 <- databel("matr")
dat3


###################################################
### code chunk number 9: intro_DatABEL.Rnw:158-159
###################################################
write.table(matr,"matr.txt",row.names=TRUE,col.names=TRUE,quote=FALSE)


###################################################
### code chunk number 10: intro_DatABEL.Rnw:162-164
###################################################
dat4 <- text2databel("matr.txt",outfile="matr1",R_matrix=TRUE,type="UNSIGNED_INT")
dat4


###################################################
### code chunk number 11: intro_DatABEL.Rnw:169-170
###################################################
dat5 <- dat4


###################################################
### code chunk number 12: intro_DatABEL.Rnw:173-175
###################################################
dat6 <- dat1[c("row1","row3"),c("col1","col2")]
dat6


###################################################
### code chunk number 13: intro_DatABEL.Rnw:188-189
###################################################
dat1[1,1] <- 321


###################################################
### code chunk number 14: intro_DatABEL.Rnw:195-196
###################################################
dat6


###################################################
### code chunk number 15: intro_DatABEL.Rnw:209-214
###################################################
dim(dat1)
length(dat1)
dimnames(dat1)
colnames(dat1)
rownames(dat1)


###################################################
### code chunk number 16: intro_DatABEL.Rnw:219-221
###################################################
dimnames(dat1) <- list(paste("ID",1:4,sep=""),paste("SNP",1:3,sep=""))
dimnames(dat1)


###################################################
### code chunk number 17: <
###################################################
backingfilename(dat1)


###################################################
### code chunk number 18: intro_DatABEL.Rnw:230-231
###################################################
cachesizeMb(dat1)


###################################################
### code chunk number 19: intro_DatABEL.Rnw:235-237
###################################################
cachesizeMb(dat1) <- 1
cachesizeMb(dat1)


###################################################
### code chunk number 20: intro_DatABEL.Rnw:246-247
###################################################
set_dimnames(dat1) <- list(dimnames(dat1)[[1]],c("duplicate","col2","duplicate"))


###################################################
### code chunk number 21: intro_DatABEL.Rnw:252-253
###################################################
dimnames(dat1)


###################################################
### code chunk number 22: intro_DatABEL.Rnw:256-257
###################################################
get_dimnames(dat1)


###################################################
### code chunk number 23: intro_DatABEL.Rnw:262-267
###################################################
disconnect(dat1)
setReadOnly(dat6) <- FALSE
dat6[1,1] <- 123
dat6
dat1


###################################################
### code chunk number 24: intro_DatABEL.Rnw:274-278
###################################################
newm <- as(dat2,"matrix")
class(newm)
class(newm[1,1])
newm


###################################################
### code chunk number 25: intro_DatABEL.Rnw:283-284
###################################################
databel2text(dat2,file="dat2.txt")


###################################################
### code chunk number 26: intro_DatABEL.Rnw:288-289
###################################################
read.table("dat2.txt")


###################################################
### code chunk number 27: intro_DatABEL.Rnw:298-300
###################################################
apply2dfo(SNP,dfodata=dat2,anFUN="sum",MAR=2)
apply2dfo(SNP,dfodata=dat2,anFUN="sum",MAR=1)


###################################################
### code chunk number 28: intro_DatABEL.Rnw:304-305
###################################################
apply2dfo(SNP^2,dfodata=dat2,anFUN="sum",MAR=2)


###################################################
### code chunk number 29: intro_DatABEL.Rnw:308-311
###################################################
Y <- rnorm(4)
apply2dfo(Y~SNP,dfodata=dat2,anFUN="lm",MAR=2)
apply2dfo(Y~SNP+I(SNP^2),dfodata=dat2,anFUN="lm",MAR=2)


###################################################
### code chunk number 30: intro_DatABEL.Rnw:322-326
###################################################
rm(list=ls())
gc()
unlink("*.fv?")
unlink("*.txt")


