### R code from vignette source 'clusterCrit.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: clusterCrit.Rnw:1342-1359
###################################################
    library(clusterCrit)
    dataPath <- file.path(path.package(package="clusterCrit"),"unitTests","data","testsInternal_400_4.Rdata")
    load(file=dataPath, envir=.GlobalEnv)
    traj <- traj_400_4
    part <- part_400_4
	critName <- "Ball_Hall"
    vals <- vector(mode="numeric",length=6)
    for (k in 2:7) {
        idx <- intCriteria(traj,part[[k]],critName)
        vals[k-1] <- idx
    }
	par(mfrow=c(1,2))
	plot(traj[,1],traj[,2],col=part[[4]],xlab="",ylab="",asp=1)
    plot(2:7,vals,type='o',main=critName,xlab="# of clusters",ylab="Index",asp=1)
    # bestidx <- bestCriterion(vals,critName)
    # points(bestidx,vals[bestidx+1],pch=8,col="green")
	par(mfrow=c(1,1))


###################################################
### code chunk number 2: clusterCrit.Rnw:1885-1886
###################################################
library(clusterCrit)


###################################################
### code chunk number 3: clusterCrit.Rnw:1890-1893
###################################################
x <- rbind(matrix(rnorm(100, mean = 0, sd = 0.5), ncol = 2),
           matrix(rnorm(100, mean = 1, sd = 0.5), ncol = 2),
           matrix(rnorm(100, mean = 2, sd = 0.5), ncol = 2))


###################################################
### code chunk number 4: clusterCrit.Rnw:1900-1901
###################################################
cl <- kmeans(x, 3)


###################################################
### code chunk number 5: clusterCrit.Rnw:1905-1906
###################################################
getCriteriaNames(TRUE)


###################################################
### code chunk number 6: clusterCrit.Rnw:1911-1914
###################################################
intIdx <- intCriteria(x,cl$cluster,"all")
length(intIdx)
intIdx[["trace_w"]]


###################################################
### code chunk number 7: clusterCrit.Rnw:1918-1919
###################################################
intCriteria(x,cl$cluster,c("C_index","Calinski_Harabasz","Dunn"))


###################################################
### code chunk number 8: clusterCrit.Rnw:1924-1925
###################################################
intCriteria(x,cl$cluster,c("det","cal","dav"))


###################################################
### code chunk number 9: clusterCrit.Rnw:1931-1933
###################################################
part1<-sample(1:3,150,replace=TRUE)
part2<-sample(1:5,150,replace=TRUE)


###################################################
### code chunk number 10: clusterCrit.Rnw:1937-1938
###################################################
getCriteriaNames(FALSE)


###################################################
### code chunk number 11: clusterCrit.Rnw:1942-1945
###################################################
extIdx <- extCriteria(part1,part2,"all")
length(extIdx)
extIdx[["jaccard"]]


###################################################
### code chunk number 12: clusterCrit.Rnw:1949-1950
###################################################
extCriteria(part1,part2,c("Rand","Folkes"))


###################################################
### code chunk number 13: clusterCrit.Rnw:1954-1955
###################################################
extCriteria(part1,part2,c("ra","fo"))


