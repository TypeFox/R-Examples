GWAS <- function(pheno,geno,fixed=NULL,K=NULL,n.PC=0,min.MAF=0.05,n.core=1,P3D=TRUE,plot=TRUE) {

qvalue <- function(p) {
	smooth.df = 3
	
    if(min(p)<0 || max(p)>1) {
      print("ERROR: p-values not in valid range.")
      return(0)
    }
    
    lambda=seq(0,0.90,0.05)
    m <- length(p)

    pi0 <- rep(0,length(lambda))
    for(i in 1:length(lambda)) {pi0[i] <- mean(p >= lambda[i])/(1-lambda[i])}

    spi0 <- smooth.spline(lambda,pi0,df=smooth.df)
    pi0 <- predict(spi0,x=max(lambda))$y
    pi0 <- min(pi0,1)

    if(pi0 <= 0) {
      print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
      return(0)
    }

	#The estimated q-values calculated here
    u <- order(p)

    # ranking function which returns number of observations less than or equal
    qvalue.rank <- function(x) {
      idx <- sort.list(x)

      fc <- factor(x)
      nl <- length(levels(fc))
      bin <- as.integer(fc)
      tbl <- tabulate(bin)
      cs <- cumsum(tbl)
 
      tbl <- rep(cs, tbl)
      tbl[idx] <- tbl

      return(tbl)
    }

    v <- qvalue.rank(p)
    
    qvalue <- pi0*m*p/v
    qvalue[u[m]] <- min(qvalue[u[m]],1)
    for(i in (m-1):1) {qvalue[u[i]] <- min(qvalue[u[i]],qvalue[u[i+1]],1)}

    return(qvalue)
}

manhattan <- function(input,fdr.level=0.05) {
	#first column is marker name
	#second is chromosome
	#third is map position
	#fourth is score
	input <- input[order(input[,2],input[,3]),]
	
	chroms <- unique(input[,2])	
	n.chrom <- length(chroms)
	chrom.start <- rep(0,n.chrom)
	chrom.mid <- rep(0,n.chrom)

	if (n.chrom > 1) {
		for (i in 1:(n.chrom-1)) {chrom.start[i+1] <- chrom.start[i]+max(input[which(input[,2]==chroms[i]),3])+1}
	}
	x.max <- chrom.start[n.chrom]+max(input[which(input[,2]==chroms[n.chrom]),3])
	plot(0,0,type="n",xlim=c(0,x.max),ylim=c(0,max(input[,4])+1),ylab="-log(p)",xlab="Chromosome",xaxt="n")

	for (i in seq(1,n.chrom,by=2)) {
		ix <- which(input[,2]==chroms[i])
		chrom.mid[i] <- median(chrom.start[i]+input[ix,3])
		points(chrom.start[i]+input[ix,3],input[ix,4],col="dark blue",pch=16)
	}	

	if (n.chrom > 1){
	for (i in seq(2,n.chrom,by=2)) {
		ix <- which(input[,2]==chroms[i])
		chrom.mid[i] <- median(chrom.start[i]+input[ix,3])
		points(chrom.start[i]+input[ix,3],input[ix,4],col="cornflowerblue",pch=16)
	}	
	}

	q.ans <- qvalue(10^-input[,4])
	temp <- cbind(q.ans,input[,4])
	temp <- temp[order(temp[,1]),]	
	if (temp[1,1]<fdr.level) {
		temp2 <- tapply(temp[,2],temp[,1],mean)
		qvals <- as.numeric(rownames(temp2))
		x <- which.min(abs(qvals-fdr.level))
		first <- max(1,x-2)
		last <- min(x+2,length(qvals))
		if ((last-first)<4) {last <- first + 3}
		splin <- smooth.spline(x=qvals[first:last],y=temp2[first:last],df=3)
		lines(x=c(0,x.max),y=rep(predict(splin,x=fdr.level)$y,2),lty=2)
	}
	axis(side=1,at=chrom.mid,labels=chroms)
}

qq <- function(scores) {
	remove <- which(scores==0)
	if (length(remove)>0) {
		x <- sort(scores[-remove],decreasing=TRUE)
	} else {
		x <- sort(scores,decreasing=TRUE)
	}
	n <- length(x)
	unif.p <- -log10(ppoints(n))
	plot(unif.p,x,pch=16,xlab="Expected -log(p)",ylab="Observed -log(p)")
	lines(c(0,max(unif.p)),c(0,max(unif.p)),lty=2)	
}

score.calc <- function(M) {
  m <- ncol(M)
  scores <- array(0,m)

for (i in 1:m) {

  	Mi <- M[,i]
	freq <- mean(Mi+1,na.rm=TRUE)/2
	MAF <- min(freq,1-freq)
    if (MAF >= min.MAF) {
		not.NA.gid <- which(!is.na(Mi))
		temp <- rep(1,length(Mi))
		temp[not.NA.gid] <- 0
		not.NA.obs <- which(Z%*%temp!=1)
		n2 <- length(not.NA.obs)
		y2 <- matrix(y[not.NA.obs],n2,1)
		Z2 <- Z[not.NA.obs,not.NA.gid]
		X3 <- cbind(X2[not.NA.obs,],Z2%*%Mi[not.NA.gid])
		p <- ncol(X3)
		v1 <- 1
		v2 <- n2-p

		if (!P3D) {	  	
			H2inv <- mixed.solve(y=y2,X=X3,Z=Z2,K=K2[not.NA.gid,not.NA.gid],return.Hinv=TRUE)$Hinv
		} else {
			H2inv <- Hinv[not.NA.obs,not.NA.obs]
		}

		W <- crossprod(X3,H2inv%*%X3)
		Winv <- try(solve(W),silent=TRUE)
		if (class(Winv)!="try-error") {
			beta <- Winv %*% crossprod(X3,H2inv%*%y2)
			resid <- y2 - X3 %*% beta
			s2 <- as.double(crossprod(resid,H2inv%*%resid))/v2
			CovBeta <- s2*Winv
			Fstat <- beta[p]^2/CovBeta[p,p]
			x <- v2/(v2+v1*Fstat)
			scores[i] <- -log10(pbeta(x,v2/2,v1/2))
		} #if try
	} #MAF
  } #for i
      
  return(scores)
} #end score.calc

make.full <- function(X) {
	svd.X <- svd(X)
	r <- max(which(svd.X$d>1e-8))
	return(as.matrix(svd.X$u[,1:r]))
}


n <- nrow(pheno)
pheno.ix <- 2:ncol(pheno)
names <- colnames(pheno)
X <- matrix(1,n,1)		
if (!is.null(fixed)) {
	p <- length(fixed)
	for (i in 1:p) {
		xpos <- match(fixed[i],names)
		pheno.ix <- setdiff(pheno.ix,xpos)
		xx <- factor(pheno[,xpos])	
		if (length(unique(xx)) > 1) {X <- cbind(X,model.matrix(~x-1,data.frame(x=xx)))}
	}
}

geno <- geno[order(geno[,2],geno[,3]),] 

M <- t(geno[,-c(1:3)])  #first column is marker name, second is chrom, third is map position
map <- geno[,1:3]
m <- ncol(M)  # number of markers

geno.gid <- colnames(geno)[-c(1:3)]
rownames(M) <- geno.gid

if (is.null(K)) {
	K <- A.mat(M,shrink=FALSE)
}

if (n.PC > 0) {
	eig.vec <- eigen(K)$vectors
}

if (length(which(rownames(K)!=geno.gid))>0) {
	stop("Line names in K and genotype file do not match.")
}

n.phenos <- length(pheno.ix)
all.scores <- matrix(0,m,n.phenos)
trait.names <- colnames(pheno)[pheno.ix]
colnames(all.scores) <- trait.names

if (plot) {
	p1 <- floor(sqrt(n.phenos))
	p2 <- ceiling(n.phenos/p1)
	par(mfrow=c(p1,p2))
}

if (n.phenos==0) {
	stop("No phenotypes.")
}
for (i in 1:n.phenos) {
	print(paste("GWAS for trait:",trait.names[i]))
	y <- pheno[,pheno.ix[i]]
	not.miss <- which(!is.na(y))
	y <- y[not.miss]
	n <- length(y)

	pheno.gid <- unique(pheno[not.miss,1])
	n.gid <- length(pheno.gid)
	ix.pheno <- match(pheno.gid,geno.gid)
	miss.pheno.gid <- which(is.na(ix.pheno))
	if (length(miss.pheno.gid)>0) {
		stop(paste("The following lines have phenotypes but no genotypes:",paste(unique(pheno.gid[miss.pheno.gid]),collapse=" ")))
	}

	Z <- matrix(0,n,length(pheno.gid))
	Z[cbind(1:n,match(pheno[not.miss,1],pheno.gid))] <- 1
	K2 <- K[ix.pheno,ix.pheno]
	if (n.PC > 0) {
		X2 <- make.full(cbind(X[not.miss,],Z%*%eig.vec[ix.pheno,1:n.PC]))
	} else {
		X2 <- make.full(X[not.miss,])
	}

	if (P3D) {
		Hinv <- mixed.solve(y,X=X2,Z=Z,K=K2,return.Hinv=TRUE)$Hinv  
		print("Variance components estimated. Testing markers.")
	}

if ((n.core > 1) & requireNamespace("parallel",quietly=TRUE)) {
    	it <- split(1:m,factor(cut(1:m,n.core,labels=FALSE)))
    	scores <- unlist(parallel::mclapply(it,function(markers){score.calc(M[ix.pheno,markers])},mc.cores=n.core))
	 } else {
    	scores <- score.calc(M[ix.pheno,])
	 }
	 if (plot) {
	 	qq(scores)
	 	title(main=trait.names[i])
	 }      
	 all.scores[,i] <- scores
}

if (plot) {
	if (length(grep("RStudio",names(dev.cur())))==0) { 
		if (dev.cur()==dev.next()) {
    	  dev.new()   
	    } else {
    	  dev.set(dev.next())
	    }
	}
	par(mfrow=c(p1,p2))
	for (j in 1:n.phenos) {
		manhattan(cbind(map,all.scores[,j]))
		title(main=trait.names[j])
	}
}

return(data.frame(map,all.scores))
  
} #end function

