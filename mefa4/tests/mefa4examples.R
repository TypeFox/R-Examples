## examples

## load library
library(mefa)
library(mefa4)
## create input data
x <- data.frame(
    sample = paste("Sample", c(1,1,2,2,3,4), sep="."),
    species = c(paste("Species", c(1,1,1,2,3), sep="."),  "zero.pseudo"),
    count = c(1,2,10,3,4,0),
    segment = letters[c(6,13,6,13,6,6)])
x
samp <- data.frame(samples=levels(x$sample), var1=1:2)
taxa <- data.frame(specnames=levels(x$species), var2=c("b","a"))
rownames(samp) <- samp$samples
rownames(taxa) <- taxa$specnames
samp
taxa
## Xtab class, counts by repetitions in RHS
(x0 <- Xtab(~ sample + species, x))
## counts by LHS and repetitions in RHS
(x1 <- Xtab(count ~ sample + species, x))
## drop all empty rows
(x2 <- Xtab(count ~ sample + species, x, cdrop=FALSE,rdrop=TRUE))
## drop all empty columns
Xtab(count ~ sample + species, x, cdrop=TRUE,rdrop=FALSE)
## drop specific columns by placeholder
Xtab(count ~ sample + species, x, cdrop="zero.pseudo")

## 3-way crosstab
x33 <- Xtab(count ~ sample + species + segment, x)

## Mefa class, standard
(x3 <- Mefa(x1, samp, taxa))
unclass(x3)
## effects of left join, NULL taxa slot, xtab is (not sparse) matrix
(x4 <- Mefa(as.matrix(x1), samp[1:2,]))
unclass(x4)
## effects of inner join (intersect)
(x5 <- Mefa(x2, samp, taxa, join="inner"))
unclass(x5)
unclass(Mefa(x1, samp[1:2,], join="inner"))
## xtab only Mefa
(x6 <- Mefa(x1))
## creating new Mefa object without Mefa()
new("Mefa", xtab=x1, samp=samp, taxa=taxa,join="left")

## 0 and 1 row/col Mefa object
x3[c(FALSE,FALSE,FALSE,FALSE),c(FALSE,FALSE,FALSE,FALSE)]
x3[c(TRUE,FALSE,FALSE,FALSE),c(FALSE,FALSE,FALSE,FALSE)]
x3[c(FALSE,FALSE,FALSE,FALSE),c(TRUE,FALSE,FALSE,FALSE)]
x3[c(TRUE,FALSE,FALSE,FALSE),c(TRUE,FALSE,FALSE,FALSE)]

## Melt
x0 <- Xtab(count ~ sample + species, x)
x33 <- Xtab(count ~ sample + species + segment, x)
(M1 <- Melt(x0))
Melt(as.matrix(x0))
(M2 <- Melt(x33))
stopifnot(identical(Xtab(value ~ rows + cols, M1), x0))
stopifnot(identical(Xtab(value ~ rows + cols + segm, M2), x33))

## stack
stack(x3)

## accessing the xtab slot
xtab(x3)
## replacing the slot value
x1[3,1] <- 999
xtab(x3) <- x1
xtab(x3)

## accessing and replacing the samp slot
samp(x3)
samp(x3) <- NULL
samp(x3)
samp(x3) <- samp[1:3,]
samp(x3)

## accessing and replacing the taxa slot
taxa(x3)
taxa(x3) <- NULL
taxa(x3)
taxa(x3) <- taxa[1:3,]
taxa(x3)

## subsetting
unclass(x3[3:2, 1:2])
unclass(x3[3:2,])
unclass(x3[,1:2])

## simple methods, dim, dimnames
dim(x5)
dimnames(x5)
dn <- list(paste("S", 1:dim(x5)[1], sep=""),
    paste("SPP", 1:dim(x5)[2], sep=""))
dimnames(x5) <- dn
unclass(x5)
dimnames(x5)[[1]] <- paste("S", 1:dim(x5)[1], sep="_")
unclass(x5)
dimnames(x5)[[2]] <- paste("SPP", 1:dim(x5)[2], sep="_")
unclass(x5)

## transpose
x5
t(x5)
unclass(x5)
unclass(t(x5))

## aggregation, sums
groupSums(as.matrix(x2), 1, c(1,1,2))
groupSums(as.matrix(x2), 2, c(1,1,2,2))
groupSums(x2, 1, c(1,1,2))
groupSums(x2, 2, c(1,1,2,2))
groupSums(x5, 1, c(1,1,2))
groupSums(x5, 2, c(1,1,2,2))
## aggregation, means
groupMeans(as.matrix(x2), 1, c(1,1,2))
groupMeans(as.matrix(x2), 2, c(1,1,2,2))
groupMeans(x2, 1, c(1,1,2))
groupMeans(x2, 2, c(1,1,2,2))
groupMeans(x5, 1, c(1,1,2))
groupMeans(x5, 2, c(1,1,2,2))

## back and foth -- this is important to keep in the vignette
as.stcs(x1)
as.mefa(x1)
as.stcs(x3)
(a <- as.mefa(x3))
xtab(a)
samp(a)
taxa(a)
segm(a)
segm(x3)
as.Mefa(a)
as.Xtab(a)
(s <- melt(a))
as.Xtab(s)
as.Mefa(s)
melt(x1)
melt(x3)
## sparse matrix list
as.mefa(x33)

## mbind
x=matrix(1:4,2,2)
rownames(x) <- c("a","b")
colnames(x) <- c("A","B")
y=matrix(11:14,2,2)
rownames(y) <- c("b","c")
colnames(y) <- c("B","C")

sampx <- data.frame(x1=1:2, x2=2:1)
rownames(sampx) <- rownames(x)
sampy <- data.frame(x1=3:4, x3=10:11)
rownames(sampy) <- rownames(y)
taxay <- data.frame(x1=1:2, x2=2:1)
rownames(taxay) <- colnames(y)
taxax <- NULL

mbind(x,y)
mbind(as(x,"sparseMatrix"),as(y,"sparseMatrix"))
unclass(mbind(Mefa(x,sampx),Mefa(y,sampy,taxay)))

## 1x1 cases
unclass(Mefa(y,sampy,taxay)[1,1])
unclass(Mefa(y,sampy[,1,drop=FALSE],taxay[,1,drop=FALSE]))
unclass(Mefa(y,sampy[,1,drop=FALSE],taxay[,1,drop=FALSE])[1,1])

z0 <- Mefa(y)
samp(z0) <- sampy[,1,drop=FALSE]
taxa(z0) <- taxay[,1,drop=FALSE]
unclass(z0)
xtab(z0) <- y[1,1,drop=FALSE]
unclass(z0)

