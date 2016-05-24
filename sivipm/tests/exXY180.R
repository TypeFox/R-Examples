# Example exXY180: transformed data,  monomes description on a file.
# Qualitative variables split into 0-1 or -1/+1 variables.
# 2 responses. It is exX18Y2, with transformed data
library("sivipm")

# -------------------------------------
# READ DATA (transformed data)

#  XY180 <-read.table("../inst/extdata/XY180.txt",  header=TRUE)

datancol <- ncol(XY180)
nX=datancol -2
X180 = XY180[,1:nX]
Y180 = XY180[,(nX+1):datancol]
# -------------------------------------
#  MONOME DESCRIPTION: monomes coded by the inputs numbers
Pexp <- as.character(1:18)
# derniers: 17, 18
for (i in 1:13) {
  Pexp <- c(Pexp, paste(i,"*",i, sep=""))
}
# dernier 13*13
for (i in 1:13) {
    j <- Pexp[18+i]
     Pexp <-  c(Pexp, paste(j,"*",i, sep=""))
}
# derniers: 31*13
for (i in 1:12 ){
  for (j in (i+1):13) {
    Pexp <- c(Pexp, paste(i, "*", j,sep=""))
   }
}
# derniers: 12*13

for (i in 1:18) {
    if (i != 15)
    Pexp <- c(Pexp, paste("15*",i,sep=""))
  }
# derniers: 15*14, 15*16, 15*17, 15*18
for (i in 16:18) {
    for (j in 9:11) {
      Pexp <- c(Pexp, paste(i, "*", j,sep=""))
   }
}
Pexp <- c(Pexp, c("13*16"), c("13*17"), c("13*18"),
             c("14*16"), c("14*17"), c("14*18"),
            c("14*9"), c("14*10"), c("12*14"), c("2*14"))

# X18Y2 <-read.table("../inst/extdata/X18Y2.txt",  header=TRUE)
varnames <- colnames(X18Y2)[1:18]
zz <- vect2polyXT(varnames , X180, Pexp)


# -------------------------------------
# CALCULATIONS 
print("TSIVIP without alea on the first response")
nc <- 21# idem que avant-derniere page du raptech 2012
nc <- 2 # plus rapide
res <- sivipm(as.matrix(Y180[,1], ncol=1), zz,nc)
print(res, all=TRUE)
# idem que avant-derniere page du raptech 2012

# -------------------------------------
# COMPARE WITH exX18Y2 (raw data provided)
 nvar <- length(varnames) # number of "raw" variables
datancol <- ncol(X18Y2)
X = X18Y2[,1:nvar]
Y = X18Y2[,(nvar+1):datancol]

#X[X[,"MODEL"]==0, "MODEL"]=-1

zz=vect2polyX(X,Pexp)
res2 <- sivipm(Y[,1], zz, nc)
res@fo.isivip <- res@fo.isivip[1:nvar]
names(res2@fo.isivip) = names(res@fo.isivip)
print(all.equal(res, res2))

                 
