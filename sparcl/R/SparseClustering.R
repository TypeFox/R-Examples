# this file contains the hidden functions for the sparcl package
 
GetWCSS <- function(x, Cs, ws=NULL){
  wcss.perfeature <- numeric(ncol(x))
  for(k in unique(Cs)){
    whichers <- (Cs==k)
    if(sum(whichers)>1) wcss.perfeature <- wcss.perfeature + apply(scale(x[whichers,],center=TRUE, scale=FALSE)^2, 2, sum)
  }
  bcss.perfeature <- apply(scale(x, center=TRUE, scale=FALSE)^2, 2, sum)-wcss.perfeature
  if(!is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), wcss.ws=sum(wcss.perfeature*ws),
                               bcss.perfeature=bcss.perfeature))
  if(is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), bcss.perfeature=bcss.perfeature))
}
 
UpdateCs <- function(x, K, ws, Cs){
  x <- x[,ws!=0]
  z <- sweep(x, 2, sqrt(ws[ws!=0]), "*")
  nrowz <- nrow(z)
  mus <- NULL
  if(!is.null(Cs)){
    for(k in unique(Cs)){
      if(sum(Cs==k)>1) mus <- rbind(mus, apply(z[Cs==k,],2,mean))
      if(sum(Cs==k)==1) mus <- rbind(mus, z[Cs==k,])
    }
  }
  if(is.null(mus)){
    km <- kmeans(z, centers=K, nstart=10)
  } else {
    distmat <- as.matrix(dist(rbind(z, mus)))[1:nrowz, (nrowz+1):(nrowz+K)]
    nearest <- apply(distmat, 1, which.min)
    if(length(unique(nearest))==K){
      km <- kmeans(z, centers=mus)
    } else {
      km <- kmeans(z, centers=K, nstart=10)
    }
  }
  return(km$cluster)
}

#distmat <- function(x){
#  return(dist(x))
#}

BinarySearch <- function(argu,sumabs){
  if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu))-1e-5
  iter <- 1
  while(iter<=15 && (lam2-lam1)>(1e-4)){
    su <- soft(argu,(lam1+lam2)/2)
    if(sum(abs(su/l2n(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    iter <- iter+1
  }
  return((lam1+lam2)/2)
}

soft <- function(x,d){
  return(sign(x)*pmax(0, abs(x)-d))
}

l2n <- function(vec){
  return(sqrt(sum(vec^2)))
}

GetUW <- function(ds, wbound,niter,uorth,silent){
  # uorth would be a $n^2 x k$ matrix containing $k$ previous
  # dissimilarity matrices found, if we want to do sparse comp clust
  # after having already done $k$ of these things
  # Example:
  # out <- HierarchicalSparseCluster(x,wbound=5)
  # out2 <- HierarchicalSparseCluster(x,wbound=5, uorth=out$u)
  # Then out2 contains a sparse complementary clustering
  p <- ncol(ds)
  w <- rep(1/p, p)*wbound
  iter <- 1
  if(!is.null(uorth)){
    if(sum(abs(uorth-t(uorth)))>1e-10) stop("Uorth must be symmetric!!!!!!!!!!")
    uorth <- matrix(uorth[lower.tri(uorth)],ncol=1)
    uorth <- uorth/sqrt(sum(uorth^2))
  }
  u <- rnorm(nrow(ds))
  w <- rep(1, ncol(ds))
  w.old <- rnorm(ncol(ds))
#  u.old <- rnorm(nrow(ds))
  while(iter<=niter && (sum(abs(w.old-w))/sum(abs(w.old)))>1e-4){#(sum(abs(u.old-u))/sum(abs(u.old)))>1e-4){ # was 1e-3 and involved ws until 11.12.09
    if(!silent) cat(iter,fill=FALSE)
#    u.old <- u
    if(iter>1) u <- ds[,argw>=lam]%*%matrix(w[argw>=lam],ncol=1) 
    if(iter==1) u <- ds%*%matrix(w,ncol=1)
    if(!is.null(uorth)) u <- u - uorth%*%(t(uorth)%*%u)
    iter <- iter+1
    u <- u/l2n(u)
    w.old <- w
    argw <- matrix(pmax(matrix(u,nrow=1)%*%ds,0),ncol=1) 
    lam <- BinarySearch(argw,wbound)
    w <- soft(argw,lam) # Don't need to normalize b/c this will be done soon
    w <- w/l2n(w)
  } 
  u <-  ds[,argw>=lam]%*%matrix(w[argw>=lam],ncol=1)/sum(w)
  if(!is.null(uorth)) u <- u - uorth%*%(t(uorth)%*%u)
  u <- u/l2n(u)
  w <- w/l2n(w)
  crit <- sum(u*(ds%*%matrix(w,ncol=1)))
  u2 <- matrix(0,nrow=ceiling(sqrt(2*length(u))),ncol=ceiling(sqrt(2*length(u))))
  u2[lower.tri(u2)] <- u
  u <- as.matrix(as.dist(u2))/sqrt(2)
  return(list(u=u, w=w, crit=crit))
}

 
UpdateWs <- function(x, Cs, l1bound){
  wcss.perfeature <- GetWCSS(x, Cs)$wcss.perfeature
  tss.perfeature <- GetWCSS(x, rep(1, nrow(x)))$wcss.perfeature
  lam <- BinarySearch(-wcss.perfeature+tss.perfeature, l1bound)
  ws.unscaled <- soft(-wcss.perfeature+tss.perfeature,lam)
  return(ws.unscaled/l2n(ws.unscaled))
}


output.cluster.files.fun <- function(x,out,outputfile.prefix,genenames=NULL,genedesc=NULL){
         p=ncol(x)
         n=nrow(x)
         geneids=dimnames(x)[[2]]
         samplenames=dimnames(x)[[1]]
         if(is.null(geneids)) geneids <- paste("Gene", 1:ncol(x))
         if(is.null(samplenames)) samplenames <- paste("Sample",1:nrow(x))
         if(is.null(genenames)){genenames=geneids}
         if(is.null(genedesc)){genedesc <- rep("", ncol(x))}

         xx=x[,out$ws!=0]
         geneids.subset=geneids[out$ws!=0]
         genenames.subset=genenames[out$ws!=0]
         genedesc.subset=genedesc[out$ws!=0]
         pp=ncol(xx)
         sample.order=out$hc$order
          samplenames.o=samplenames[sample.order]
         arrynames=paste("ARRY",as.character(sample.order),"X",sep="")
         feature.order=1:pp
   if(!is.null(out$hc.features)){feature.order=out$hc.features$order}
    xx.o=xx[sample.order,feature.order]
    geneids.subset.o=geneids.subset[feature.order]
    genenames.subset.o=genenames.subset[feature.order]
    genedesc.subset.o=genedesc.subset[feature.order]
    genex=paste("GENE",as.character(1:pp),"X",sep="")
    genex.o=genex[feature.order]
    arrynames.o=arrynames[sample.order]

# output cdt
	file=paste(outputfile.prefix,".cdt",sep="")
         xprefix <- rbind(c("GID","UID","NAME","GWEIGHT",samplenames.o),c("AID","","","",arrynames))
         xbig <- matrix(NA, ncol=(n+4),nrow=pp)
         for(i in 1:pp){
           xbig[i,] <- c(genex.o[i],geneids.subset.o[i],genedesc.subset.o[i],"1",xx.o[,i]) # was xx
         }
         xbig <- rbind(xprefix,xbig)
         write.table(file=file,xbig,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
	
	#output atr file
	atr=out$hc$merge
	atr.new=atr
	for(i in 1:nrow(atr)){
	 for(j in 1:2){
	   if(atr[i,j]<0){atr.new[i,j]=paste("ARRY",as.character(abs(atr[i,j])),"X",sep="")}
	   if(atr[i,j]>0){atr.new[i,j]=paste("NODE",as.character(abs(atr[i,j])),"X",sep="")}
	}}
        col1=paste("NODE",as.character(1:nrow(atr.new)),"X",sep="")
	atr.new=cbind(col1,atr.new,1-out$hc$height/2)
	output.matrix(atr.new, paste(outputfile.prefix,".atr",sep=""))
	
	if(!is.null(out$hc.features)){
	#output gtr file
	gtr=out$hc.features$merge
	gtr.new=gtr
	for(i in 1:nrow(gtr)){
	 for(j in 1:2){
	   if(gtr[i,j]<0){gtr.new[i,j]=paste("GENE",as.character(abs(gtr[i,j])),"X",sep="")}
	   if(gtr[i,j]>0){gtr.new[i,j]=paste("NODE",as.character(abs(gtr[i,j])),"X",sep="")}
	}}
        col1=paste("NODE",as.character(1:nrow(gtr.new)),"X",sep="")
	gtr.new=cbind(col1,gtr.new,1-out$hc.features$height/2)
	output.matrix(gtr.new, paste(outputfile.prefix,".gtr",sep=""))
}

return()
}


output.matrix <- function(x,file){
  write.table(file=file,x,quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}



read.gct <- function(file) {
        if (is.character(file))
        if (file == "")
            file <- stdin()
        else {
            file <- file(file, "r")
            on.exit(close(file))
        }
        if (!inherits(file, "connection"))
        stop("argument `file' must be a character string or connection")

   # line 1 version
        version <- readLines(file, n=1)

        # line 2 dimensions
        dimensions <- scan(file, what=list("integer", "integer"), nmax=1, quiet=TRUE)
        rows <- dimensions[[1]]
        columns <- dimensions[[2]]
        # line 3 Name\tDescription\tSample names...
        column.names <- read.table(file, header=FALSE, quote='', nrows=1, sep="\t", fill=FALSE, comment.char='')
        column.names <-column.names[3:length(column.names)]


        if(length(column.names)!=columns) {
                stop(paste("Number of sample names", length(column.names), "not equal to the number of columns", columns, "."))
        }

        colClasses <- c(rep(c("character"), 2), rep(c("double"), columns))

        x <- read.table(file, header=FALSE, quote="", row.names=NULL, comment.char="", sep="\t", colClasses=colClasses, fill=FALSE)
        row.descriptions <- as.character(x[,2])
        data <- as.matrix(x[seq(from=3, to=dim(x)[2], by=1)])

        column.names <- column.names[!is.na(column.names)]

        colnames(data) <- column.names
        row.names(data) <- x[,1]
        return(list(row.descriptions=row.descriptions, data=data))
}


extract.prefix <- function(file){
#  i=0
#  while(substring(file,i,i)!="." & (i <nchar(file))) {i=i+1}
#  if(i==nchar(file)){stop('Error in file name')}
#  pre=substring(file,1,i-1)
#  return(pre)
  tmp <- strsplit(file,"\\.")[[1]]
  return(paste(tmp[-length(tmp)],collapse="."))
}
