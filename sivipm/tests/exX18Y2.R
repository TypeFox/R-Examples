# Example exX18Y2: 18 variables including qualitative variables
# split into 0-1  or -1/+1 variables. No missing. 2 responses.
# Raw data are on file X18Y2.txt
# Corresponding transformed data  are on  file XY180.txt


library("sivipm")

# -------------------------------------
# READ DATA
# X18Y2 <-read.table("../inst/extdata/X18Y2.txt",  header=TRUE)
datancol <- ncol(X18Y2)
nvar=datancol -2
X = X18Y2[,1:nvar]
Y = X18Y2[,(nvar+1):datancol]
varnamesX18 <- colnames(X)

# -------------------------------------
#  CREATE POLYNOME
# monomes <- scan("../inst/extdata/mononamesXY180.txt", what=character(), sep="\n")
monomes <- scan(system.file("extdata", "mononamesXY180.txt",package="sivipm"),  what=character(), sep="\n")
P <- vect2polyX(X, monomes)

# -------------------------------------
# CALCULATIONS
print(" without alea")
nc=2
# too long: print(sivipm(Y, P,nc))
res1 <- sivipm(Y[,1], P,nc, output="Q2")
print(res1, all=TRUE)

# compute tsivip Y by Y
nc <- 2
print(apply(Y, 2, sivipm, P, nc))
# -------------------------------------
# COMPARISON: are the results the same when the data are provided
# in their transformed version?
# XY180 <-read.table("../inst/extdata/XY180.txt",  header=TRUE)

datancol <- ncol(XY180)
nX=datancol -2
dataexp = XY180[,1:nX]
Y = XY180[,(nX+1):datancol]
Pext <- vect2polyXT(varnamesX18, dataexp, monomes)
res2 <- sivipm(Y[,1], Pext,nc)
print(res2, all=TRUE)
res2@fo.isivip <- res2@fo.isivip[1:18]
names(res2@fo.isivip) = names(res1@fo.isivip)
print(all.equal(res1, res2))
