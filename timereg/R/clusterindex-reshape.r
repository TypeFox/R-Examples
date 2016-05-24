cluster.index <- function(clusters,index.type=FALSE,num=NULL,Rindex=0)
{ ## {{{
antpers <- length(clusters)

if (index.type==FALSE)  {
	if (is.numeric(clusters)) clusters <-  sindex.prodlim(unique(clusters),clusters)-1 else  {
	   max.clust <- length(unique(clusters))
	   clusters <- as.integer(factor(clusters, labels = 1:max.clust))-1
	}
}

 nclust <- .C("nclusters",
	as.integer(antpers), as.integer(clusters), as.integer(rep(0,antpers)), 
	as.integer(0), as.integer(0), package="timereg")
  maxclust <- nclust[[5]]
  antclust <- nclust[[4]]
  cluster.size <- nclust[[3]][1:antclust]

if ((!is.null(num))) { ### different types in different columns
   mednum <- 1
if (is.numeric(num)) numnum <-  sindex.prodlim(unique(num),num)-1
else numnum <- as.integer(factor(num, labels = 1:maxclust)) -1
maxclust <- max(numnum)+1; 
} else { numnum <- 0; mednum <- 0; }

init <- -1*Rindex
clustud <- .C("clusterindex",
	      as.integer(clusters), as.integer(antclust),
	      as.integer(antpers), as.integer(rep(init,antclust*maxclust)),
	      as.integer(rep(0,antclust)), as.integer(mednum), 
	      as.integer(numnum),as.integer(rep(0,antclust)), package="timereg")

if (Rindex==1) idclust  <- matrix(clustud[[4]],antclust,maxclust)+1
else idclust <- matrix(clustud[[4]],antclust,maxclust)
if(Rindex==1) idclust[idclust==0] <- NA 
if (Rindex==1) firstclustid <- clustud[[8]]+1 else firstclustid <- clustud[[8]]
out <- list(clusters=clusters,maxclust=maxclust,antclust=antclust,idclust=idclust,
	    cluster.size=cluster.size,firstclustid=firstclustid)
} ## }}}

###faster.reshape <- function(data,clusters,index.type=FALSE,num=NULL,Rindex=1)
###{ ## {{{
###data <- as.matrix(data)
###if (NCOL(data)==1) data <- cbind(data)
###
###antpers <- length(clusters)
###if (index.type==FALSE)  {
###	max.clust <- length(unique(clusters))
###	clusters <- as.integer(factor(clusters, labels = 1:max.clust))-1 
###}
###
### nclust <- .C("nclusters",
###	as.integer(antpers), as.integer(clusters), as.integer(rep(0,antpers)), 
###	as.integer(0), as.integer(0), package="timereg")
###  maxclust <- nclust[[5]]
###  antclust <- nclust[[4]]
###  cluster.size <- nclust[[3]][1:antclust]
###
###if ((!is.null(num)) && (Rindex==1)) { ### different types in different columns
###   mednum <- 1
###   numnum <- as.integer(factor(num, labels = 1:maxclust)) -1
###} else { numnum <- 0; mednum <- 0; }
###
###data[is.na(data)] <- nan 
###p <- ncol(data); 
###init <- -1*Rindex;
###clustud <- .C("clusterindexdata",
###	        as.integer(clusters), as.integer(antclust),as.integer(antpers),
###                as.integer(rep(init,antclust*maxclust)),as.integer(rep(0,antclust)), as.integer(mednum), 
###		as.integer(numnum), as.double(c(data)), 
###		as.integer(p), as.double(rep(init*1.0,antclust*maxclust*p)), package="timereg")
###idclust <- matrix(clustud[[4]],antclust,maxclust)
###xny <- matrix(clustud[[10]],antclust,maxclust*p)
###if(Rindex==1) xny[idclust==-1] <- NA 
###if(Rindex==1) xny[idclust==-1] <- NA 
###if(Rindex==1) idclust[idclust==-1] <- NA 
###  mnames <- c()
###print(maxclust)
###  for (i in 1:maxclust) {
###     mnames <- c(mnames,paste(names(data),".",i,sep=""))
###  }
###  xny <- data.frame(xny)
###  names(xny) <- mnames
###out <- xny; 
###} ## }}}
###
###fast.reshape <- function(data,id=id,num=NULL) {  ## {{{
###   if (NCOL(data)==1) data <- cbind(data)
###   cud <- cluster.index(id,num=num,Rindex=1) ## NA for index not there, index starts at 0 for use in C
###   dataw <- c()
###   mnames <- c()
###  for (i in 1:cud$maxclust) {
###     if (i==1) dataw <- data[cud$idclust[,i]+1,]
###     else dataw <- cbind(dataw,data[cud$idclust[,i]+1,])
###     mnames <- c(mnames,paste(names(data),".",i,sep=""))
###  }
###  names(dataw) <- mnames
###  return(dataw)
###}  ## }}}
