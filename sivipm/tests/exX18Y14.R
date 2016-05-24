# Example ex18Y14: the same 18 inputs as in X18Y2 but with 14 responses
library("sivipm")

# -------------------------------------
# READ DATA 

# X18Y14 <-read.table("../inst/extdata/X18Y14.txt", header=TRUE, na.strings =".", colClasses  = "numeric" )

X18Y14 <-read.table(system.file("extdata", "X18Y14.txt", package="sivipm"),  header=TRUE, na.strings =".", colClasses  = "numeric" )

datancol <- ncol(X18Y14)
nX=18
X = X18Y14[,1:nX]
Y = X18Y14[,(nX+1):datancol]
# -------------------------------------
# CALCULATIONS WITH POLYNOME DESCRIPTION NOT REQUIRED
nc=2
b <- new("polyX", dataX.exp=X)
A <-  sivipm(Y, b, nc=nc)
print( A, all=TRUE)

# -------------------------------------
#  MONOME DESCRIPTION: monomes coded by the inputs numbers
Pexp <- as.character(1:18)
for (i in 1:13) {
  Pexp <- c(Pexp, paste(i,"*",i, sep=""))
}
for (i in 1:13) {
  Pexp <- c(Pexp, paste(i,"*",i, "*",i, sep=""))
}
for (i in 1:12) {
  for (j in (i+1):13) {
     Pexp <-  c(Pexp, paste(i,"*",j, sep=""))
   }
}
for (i in 1:18) {
    if (i != 15)
    Pexp <- c(Pexp, paste("15*",i,sep=""))
  }
for (i in 16:18) {
    for (j in 9:11) {
      Pexp <- c(Pexp, paste(i, "*", j,sep=""))
   }
}
Pexp <- c(Pexp, c("13*16"), c("13*17"), c("13*18"),
             c("14*16"), c("14*17"), c("14*18"),
            c("14*9"), c("14*10"), c("12*14"), c("2*14"))



zz <- vect2polyX(X, Pexp)

# -------------------------------------
# CALCULATIONS
res1 <- sivipm(Y[,1], zz,nc)
print(res1, all=TRUE)
