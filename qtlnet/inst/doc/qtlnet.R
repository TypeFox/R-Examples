### R code from vignette source 'qtlnet.Rnw'

###################################################
### code chunk number 1: qtlnet.Rnw:44-45
###################################################
library(qtlnet)


###################################################
### code chunk number 2: qtlnet.Rnw:50-51
###################################################
example(acyclic)


###################################################
### code chunk number 3: qtlnet.Rnw:56-57
###################################################
example(cyclica)


###################################################
### code chunk number 4: qtlnet.Rnw:62-63
###################################################
example(cyclicb)


###################################################
### code chunk number 5: qtlnet.Rnw:68-69
###################################################
example(cyclicc)


###################################################
### code chunk number 6: qtlnet.Rnw:74-75
###################################################
example(glxnet)


###################################################
### code chunk number 7: qtlnet.Rnw:82-83
###################################################
library(qtlnet)


###################################################
### code chunk number 8: qtlnet.Rnw:85-92
###################################################
# Make width of chunks 60.
options(width=60)
if(!file.exists("qdgPDF")) {
  dir.create("qdgPDF")
  warning(paste("Creating Sweave directory qdgPDF"),
    call. = FALSE, immediate. = TRUE)
}


###################################################
### code chunk number 9: qtlnet.Rnw:98-99
###################################################
mymap <- sim.map(len=rep(100,20), n.mar=10, eq.spacing=FALSE, include.x=FALSE)


###################################################
### code chunk number 10: qtlnet.Rnw:104-106
###################################################
n.ind <- 200
mycross <- sim.cross(map=mymap, n.ind=n.ind, type="f2")


###################################################
### code chunk number 11: qtlnet.Rnw:111-112
###################################################
mycross <- sim.geno(mycross,n.draws=1)


###################################################
### code chunk number 12: qtlnet.Rnw:117-135
###################################################
genotypes <- pull.geno(mycross)
geno.names <- dimnames(genotypes)[[2]]
m1 <- sample(geno.names,2,replace=FALSE)
m2 <- sample(geno.names,2,replace=FALSE)
m3 <- sample(geno.names,2,replace=FALSE)
m4 <- sample(geno.names,2,replace=FALSE)

## get marker genotypes
g11 <- genotypes[,m1[1]]; g12 <- genotypes[,m1[2]]
g21 <- genotypes[,m2[1]]; g22 <- genotypes[,m2[2]]
g31 <- genotypes[,m3[1]]; g32 <- genotypes[,m3[2]]
g41 <- genotypes[,m4[1]]; g42 <- genotypes[,m4[2]]

## generate phenotypes
y1 <- runif(3,0.5,1)[g11] + runif(3,0.5,1)[g12] + rnorm(n.ind)
y2 <- runif(3,0.5,1)[g21] + runif(3,0.5,1)[g22] + rnorm(n.ind)
y3 <- runif(1,0.5,1) * y1 +  runif(1,0.5,1) * y2 + runif(3,0.5,1)[g31] + runif(3,0.5,1)[g32] + rnorm(n.ind)
y4 <- runif(1,0.5,1) * y3 + runif(3,0.5,1)[g41] + runif(3,0.5,1)[g42] + rnorm(n.ind)


###################################################
### code chunk number 13: qtlnet.Rnw:140-141
###################################################
mycross$pheno <- data.frame(y1,y2,y3,y4)


###################################################
### code chunk number 14: qtlnet.Rnw:146-148
###################################################
markers <- list(m1,m2,m3,m4)
names(markers) <- c("y1","y2","y3","y4")


###################################################
### code chunk number 15: qtlnet.Rnw:153-163
###################################################
allqtls <- list()
m1.pos <- find.markerpos(mycross, m1)
allqtls[[1]] <- makeqtl(mycross, chr = m1.pos[,"chr"], pos = m1.pos[,"pos"])
m2.pos <- find.markerpos(mycross, m2)
allqtls[[2]] <- makeqtl(mycross, chr = m2.pos[,"chr"], pos = m2.pos[,"pos"])
m3.pos <- find.markerpos(mycross, m3)
allqtls[[3]] <- makeqtl(mycross, chr = m3.pos[,"chr"], pos = m3.pos[,"pos"])
m4.pos <- find.markerpos(mycross, m4)
allqtls[[4]] <- makeqtl(mycross, chr = m4.pos[,"chr"], pos = m4.pos[,"pos"])
names(allqtls) <- c("y1","y2","y3","y4")


###################################################
### code chunk number 16: qtlnet.Rnw:168-177
###################################################
out <- qdg(cross=mycross, 
           phenotype.names = c("y1","y2","y3","y4"), 
           marker.names = markers, 
           QTL = allqtls, 
           alpha = 0.005, 
           n.qdg.random.starts=10, 
           skel.method="pcskel")

out


###################################################
### code chunk number 17: qtlnet.Rnw:182-184
###################################################
graph <- graph.qdg(out)
plot(graph)


