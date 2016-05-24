## ----setup, include=FALSE, cache=FALSE--------------------------------------------------
library(knitr)

# set global chunk options
opts_chunk$set(fig.path='figure/minimal-', fig.align='center', fig.show='hold')
options(replace.assign=TRUE,width=90)

## ---------------------------------------------------------------------------------------
library("fgpt")

## ---------------------------------------------------------------------------------------
x.coor <- 1:10
y.coor <- rep(0,10)

x.coor
y.coor

## ---------------------------------------------------------------------------------------
set.seed(23)
apply(array(x.coor, dim=c(10,5)),2,sample, size=10)

## ----out.width="300pt", out.height="300pt"----------------------------------------------
set.seed(45)
rand.sets1 <- apply(array(x.coor, dim=c(10,1000)),2,sample, size=10)
plot(as.vector(rand.sets1)+runif(10000, -0.4,0.4), rep(1:10,1000)+runif(10000,
    -0.4,0.4), col="#FF000030", pch=16, xlab="original location", 
    ylab="assigned location", main="non spatial")

## ---------------------------------------------------------------------------------------
cor.test(rep(1:10,1000),as.vector(rand.sets1))

## ---------------------------------------------------------------------------------------
set.seed(72)
xy <- cbind(x.coor,y.coor)
fgperm(xy, z=x.coor, scale=3, iter=5, as.matrix=TRUE)

## ----out.width="300pt", out.height="300pt"----------------------------------------------
set.seed(24)
rand.sets2 <- fgperm(xy, z=x.coor,scale=3, iter=1000) 
plot(unlist(rand.sets2)+runif(10000, -0.4,0.4), rep(1:10,1000)+runif(10000,
    -0.4,0.4), col="#FF000030", pch=16, xlab="original location",
    ylab="assigned location", main="FGPT: scale = 3")

## ---------------------------------------------------------------------------------------
cor.test(rep(1:10,1000),unlist(rand.sets2))

## ----out.width="200pt", out.height="200pt", fig.width=4, fig.height=4-------------------
hist(unlist(rand.sets2)[which(rep(1:10,1000)==1)], breaks=0:10+0.5,
     xlab="assigned location", main="observed location 1")

## ----out.width="200pt", out.height="200pt", fig.width=4, fig.height=4-------------------
hist(unlist(rand.sets2)[which(rep(1:10,1000)==5)], breaks=0:10+0.5,
     xlab="assigned location", main="observed location 5")

## ---------------------------------------------------------------------------------------
set.seed(7)
xy <- cbind(x.coor,y.coor)
fgperm(xy, z=x.coor, scale=6, iter=5, as.matrix=TRUE)

## ----out.width="300pt", out.height="300pt"----------------------------------------------
rand.sets3 <- fgperm(xy, z=x.coor,scale=6, iter=1000, add.obs=FALSE) 
plot(unlist(rand.sets3)+runif(10000, -0.4, 0.4), 
     rep(1:10,1000)+runif(10000, -0.4, 0.4), col="#FF000030", pch=16, 
     xlab="original location", ylab="assigned location", main="FGPT: scale = 6")

cor.test(rep(1:10,1000),unlist(rand.sets3))

## ----out.width="300pt", out.height="300pt"----------------------------------------------
set.seed(6)
xy <- cbind(x.coor,y.coor)
fgperm(xy, z=x.coor, scale=1000, iter=5, as.matrix=TRUE)

rand.sets4 <- fgperm(xy, z=x.coor,scale=1000, iter=1000, add.obs=FALSE) 
plot(unlist(rand.sets4)+runif(10000, -0.4,0.4),
     rep(1:10,1000)+runif(10000, -0.4,0.4), col="#FF000030",
     pch=16, xlab="original location", ylab="assigned location", main="FGPT: scale = 1000")
cor.test(rep(1:10,1000),unlist(rand.sets4))

## ---------------------------------------------------------------------------------------
x.coor <- rep(1:10, times=10)
y.coor <- rep(1:10, each=10)

## ----out.width="200pt", out.height="200pt", fig.width=4, fig.height=4-------------------
set.seed(29)
mark1 <- x.coor+rnorm(100,10,2)
mark2 <- y.coor+rnorm(100,10,2)
plot(x.coor, y.coor, cex=(mark1-5)/10, pch=16, col="red",main="mark 1")
plot(x.coor, y.coor, cex=(mark2-5)/10, pch=16, col="red",main="mark 2")

## ----out.width="200pt", out.height="200pt", fig.width=4, fig.height=4-------------------
set.seed(9)
mark3 <- mark1+rnorm(100,0,0.5)+11-2*x.coor
plot(x.coor, y.coor, cex=(mark3-3)/10, pch=16, col="red",main="mark 3")
points(x.coor, y.coor, cex=(mark1-3)/10)
plot(mark1,mark3, pch=16, col="red", main="Mark 1 versus mark 3")

## ---------------------------------------------------------------------------------------
set.seed(435)
xy <- cbind(x.coor,y.coor)
rand1 <- fgperm(xy=xy, scale=3)
rand2 <- fgperm(xy=xy, z=mark1, scale=3)

rand1[[1]][1:10]
rand2[[1]][1:10]

## ---------------------------------------------------------------------------------------
set.seed(5)
replace <- function(x){x[sample.int(length(x),replace=TRUE)]}
rand3 <- fgperm(xy=xy, scale=3, FUN=replace)

rand3[[1]][1:10]

## ---------------------------------------------------------------------------------------
set.seed(7)
rand4 <- fgploc(xy=xy, scale=3, marks=mark1, FUN.mani=scale)

rand4[[1]][1:10]

## ---------------------------------------------------------------------------------------
round(tapply(rand4[[1]],x.coor, mean),2)
round(tapply(rand4[[1]],x.coor, mean, na.rm=TRUE),2)
round(tapply(rand2[[1]],x.coor, mean),2)

## ----out.width="200pt", out.height="200pt", fig.width=4, fig.height=4-------------------
rand1 <- fgperm(xy=xy, z=mark1, scale=6, add.obs=TRUE)

calc1 <- fgstat(rand1)
hist(calc1, xlab="mean mark1", main="scale = 6")
abline(v=calc1[1], col="red", lwd=2)
calc1[1]
quantile(calc1, probs=c(0.025,0.975))

## ---------------------------------------------------------------------------------------
dist.mat <- as.matrix(dist(xy))
diag(dist.mat) <- NA
w.mat <- ifelse(dist.mat<=max(apply(dist.mat,1,min, na.rm=TRUE)),1,0)
diag(w.mat) <- 0

w.mat[1:4,1:4]

## ---------------------------------------------------------------------------------------
mori <- function(x,w.mat){
  x <- x-mean(x)
  length(x)*sum((x %*% t(x))*w.mat,na.rm=TRUE)/(sum(w.mat,na.rm=TRUE)*sum(x^2))
  }

## ----out.width="200pt", out.height="200pt", fig.width=4, fig.height=4-------------------
calc2 <- fgstat(rand1, FUN=mori, w.mat=w.mat)
hist(calc2, xlab="Moran's I for mark1", main="scale = 6")
abline(v=calc2[1], col="red", lwd=2)
calc2[1]
quantile(calc2, probs=c(0.025,0.975))

## ----out.width="200pt", out.height="200pt", fig.width=4, fig.height=4-------------------
set.seed(34)
rand2 <- fgperm(xy=xy, scale=3, add.obs=TRUE)

calc3 <- fgstat(rand2,cbind(mark1,mark2), FUN=cor)
hist(calc3, xlab="corr. mark1 and mark2", main="scale = 3")
abline(v=calc3[1], col="red", lwd=2)
calc3[1]
quantile(calc3, probs=c(0.025,0.975))

## ----out.width="200pt", out.height="200pt", fig.width=4, fig.height=4-------------------
calc3 <- fgstat(rand2,cbind(mark2,mark3), FUN=cor)
calc4 <- fgstat(rand2,cbind(mark1,mark3), FUN=cor)
hist(calc3, xlab="corr. mark2 and mark3", main="scale = 3")
abline(v=calc3[1], col="red", lwd=2)
calc3[1]
quantile(calc3, probs=c(0.025,0.975))

hist(calc4, xlab="corr. mark1 and mark3", main="scale = 3")
abline(v=calc4[1], col="red", lwd=2)
calc4[1]
quantile(calc4, probs=c(0.025,0.975))

## ----out.width="200pt", out.height="200pt", fig.width=4, fig.height=4-------------------
rand5 <- fgperm(xy=xy, scale=20, add.obs=TRUE)
calc5 <- fgstat(rand5,cbind(mark1,mark3), FUN=cor)
hist(calc5, xlab="corr. mark1 and mark3", main="scale = 20")
abline(v=calc5[1], col="red", lwd=2)
calc5[1]
quantile(calc5, probs=c(0.025,0.975))

## ----out.width="200pt", out.height="200pt", fig.width=4, fig.height=4-------------------
set.seed(10)
spatial.dist <- as.matrix(dist(xy, diag=TRUE, upper=TRUE)) 
rel.mat <- array(rnorm(rep(1,10000),0.8-(as.vector(spatial.dist)/1000),0.02),
                 dim=c(100,100))

calc6 <- fgstat(rand1,rel.mat)
hist(calc6, xlab="relatedness", main="scale = 3")
abline(v=calc6[1], col="red", lwd=2)
calc6[1]
quantile(calc6, probs=c(0.025,0.975))

calc7 <- fgstat(rand5,rel.mat)
hist(calc7, xlab="relatedness", main="scale = 20")
abline(v=calc7[1], col="red", lwd=2)
calc7[1]
quantile(calc7, probs=c(0.025,0.975))

## ---------------------------------------------------------------------------------------
data("Gpulex")

str(Gp.xy)
str(f.pheno)
str(m.pheno)

## ---------------------------------------------------------------------------------------
set.seed(12)

fg1a <- fgeasy(xy=Gp.xy, marks=m.pheno, iter=99)  

summary(fg1a)

## ----out.width="300pt", out.height="200pt", fig.width=6, fig.height=4-------------------
plot(fg1a)

## ---------------------------------------------------------------------------------------
set.seed(9)
remove.f <- sample(200,5)
remove.m <- sample(200,5)

remove.f
remove.m

f.pheno[remove.f] <- NA
m.pheno[remove.m] <- NA

## ----out.width="300pt", out.height="200pt", fig.width=6, fig.height=4-------------------
set.seed(456)

fg1b <- fgeasy(xy=Gp.xy, marks=cbind(f.pheno,m.pheno), pairwise=TRUE, 
               correlate="pearson", iter=99)  

summary(fg1b)
plot(fg1b)

## ---------------------------------------------------------------------------------------
size.diff <- matrix(m.pheno, nrow=200, ncol=200) - 
                matrix(f.pheno, nrow=200, ncol=200,byrow = TRUE)

set.seed(99)
fg1c <- fgeasy(xy=Gp.xy, marks=size.diff, iter=99, pairwise=TRUE, correlate=FALSE)  

summary(fg1c)

## ----out.width="300pt", out.height="200pt", fig.width=6, fig.height=4-------------------
plot(fg1c)

## ---------------------------------------------------------------------------------------
data("Pmajor")

str(xy)
str(rel2)

## ----out.width="300pt", out.height="200pt", fig.width=6, fig.height=4-------------------
set.seed(78)
fg2a <- fgeasy(xy=xy, marks=rel2, iter=99, pairwise=TRUE, correlate=FALSE)  

summary(fg2a)
plot(fg2a)

## ----out.width="300pt", out.height="200pt", fig.width=6, fig.height=4-------------------
set.seed(11)
fg2b <- fgeasy(xy=xy, marks=rel2, scale.seq=c(1,100,200,300,400,500,800),
               iter=99, pairwise=TRUE, correlate=FALSE)  

summary(fg2b)
plot(fg2b, plane=1)

