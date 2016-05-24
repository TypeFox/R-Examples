`design.graeco` <-
function (trt1, trt2, serie = 2, seed = 0, kinds = "Super-Duper",randomization=TRUE)
{
number<-10
if(serie>0) number<-10^serie
r <- length(trt1)
if (floor(r/2) * 2 == r) {
if (r == 6 | r > 13) {
cat("not implemented design ", r, "x", r, ", see help(design.graeco)\n")
return()
}
}
if (seed == 0) {
genera<-runif(1)
seed <-.Random.seed[3]
}
set.seed(seed,kinds)
parameters<-list(design="graeco",trt1=trt1,trt2=trt2,r=r,serie=serie,seed=seed,kinds=kinds,randomization)
col <- rep(gl(r, 1), r)
fila <- gl(r, r)
fila <- as.character(fila)
fila <- as.numeric(fila)
plots <- fila*number+(1:r)
C1 <- data.frame(plots, row=factor(fila), col)
if ((r == 4) | (r == 8) | (r == 10) | (r == 12)) {
if (r == 4) {
c1 <- c(1,2,3,4,2,1,4,3,3,4,1,2,4,3,2,1)
c2 <- c(1,4,2,3,2,3,1,4,3,2,4,1,4,1,3,2)
}
if (r == 8) {
c1<-c(1,5,2,3,7,4,8,6,2,8,1,7,3,6,5,4,3,4,7,1,2,5,6,8,4,3,6,5,8,1,7,2,5,1,8,4,6,3,2,7,6,7,4,8,5,2,3,1,7,6,3,2,1,8,4,5,8,2,5,6,4,7,1,3)
c2<- c(1,2,3,4,5,6,7,8,2,1,7,6,8,4,3,5,3,7,1,5,4,8,2,6,4,6,5,1,3,2,8,7,5,8,4,3,1,7,6,2,6,4,8,2,7,1,5,3,7,3,2,8,6,5,1,4,8,5,6,7,2,3,4,1)
}
if (r == 10) {
c1<-c(1,5,2,8,3,10,9,4,7,6,9,2,6,3,8,4,10,5,1,7,10,9,3,7,4,8,5,6,2,1,6,10,9,4,1,5,8,7,3,2,8,7,10,9,5,2,6,1,4,3,7,8,1,10,9,6,3,2,5,4,4,1,8,2,10,9,7,3,6,5,2,3,4,5,6,7,1,8,9,10,3,4,5,6,7,1,2,9,10,8,5,6,7,1,2,3,4,10,8,9)
c2<-c(1,8,9,7,10,4,6,5,2,3,7,2,8,9,1,10,5,6,3,4,6,1,3,8,9,2,10,7,4,5,10,7,2,4,8,9,3,1,5,6,4,10,1,3,5,8,9,2,6,7,9,5,10,2,4,6,8,3,7,1,8,9,6,10,3,5,7,4,1,2,5,6,7,1,2,3,4,8,9,10,2,3,4,5,6,7,1,10,8,9,3,4,5,6,7,1,2,9,10,8)
}

if (r == 12) {
c1<-c(1,12,6,7,5,4,10,11,9,8,2,3,2,11,5,8,6,3,9,12,10,7,1,4,3,10,8,5,7,2,12,9,11,6,4,1,4,9,7,6,8,1,11,10,12,5,3,2,5,4,10,11,9,8,2,3,1,12,6,7,6,3,9,12,10,7,1,4,2,11,5,8,7,2,12,9,11,6,4,1,3,10,8,5,8,1,11,10,12,5,3,2,4,9,7,6,9,8,2,3,1,12,6,7,5,4,10,11,10,7,1,4,2,11,5,8,6,3,9,12,11,6,4,1,3,10,8,5,7,2,12,9,12,5,3,2,4,9,7,6,8,1,11,10)
c2<-c(1,2,3,4,9,10,11,12,5,6,7,8,2,1,4,3,10,9,12,11,6,5,8,7,3,4,1,2,11,12,9,10,7,8,5,6,4,3,2,1,12,11,10,9,8,7,6,5,5,6,7,8,1,2,3,4,9,10,11,12,6,5,8,7,2,1,4,3,10,9,12,11,7,8,5,6,3,4,1,2,11,1,9,10,8,7,6,5,4,3,2,1,12,11,10,9,9,10,11,12,5,6,7,8,1,2,3,4,10,9,12,11,6,5,8,7,2,1,4,3,11,12,9,10,7,8,5,6,3,4,1,2,12,11,10,9,8,7,6,5,4,3,2,1)

}
t1<-trt1
if(randomization)t1 <- sample(trt1, r, replace = FALSE)
t2<-trt2
if(randomization)t2 <- sample(trt2, r, replace = FALSE)
t1 <- t1[c1]
t2 <- t2[c2]
C1 <- data.frame(C1[, 1:3], t1,t2)
C1[, 4] <- trt1[C1[, 4]]
C1[, 5] <- trt2[C1[, 5]]
}
else
{
C2 <- C1
a <- 1:(r * r)
dim(a) <- c(r, r)
for (i in 1:r) {
for (j in 1:r) {
k <- i + j - 1
if (k > r) 
k <- i + j - r - 1
a[i, j] <- k
}
}
m<-trt1
if(randomization)m <- sample(trt1, r)
C1 <- data.frame(C1, m[a])
m<-trt2
if(randomization)m <- sample(trt2, r)
C2 <- data.frame(C2, m[a])
ntr <- length(trt1)
C1 <- data.frame(C1, B = 0)
for (k in 1:r) {
x <- C1[k, 4]
i <- 1
for (j in 1:(r^2)) {
y <- C2[(k - 1) * r + i, 4]
if (C1[j, 4] == x) {
C1[j, 5] <- y
i <- i + 1
}
}
}
C1[, 5] <- trt2[C1[, 5]]
}
C1[, 4] <- as.factor(C1[, 4])
C1[, 5] <- as.factor(C1[, 5])
names(C1)[4] <- c(paste(deparse(substitute(trt1))))
names(C1)[5] <- c(paste(deparse(substitute(trt2))))
outdesign<-list(parameters=parameters,sketch=matrix(paste(C1[,4], C1[,5]), byrow = TRUE, ncol = r),book=C1)
return(outdesign)
}
