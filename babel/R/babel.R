#library(edgeR)
#library(parallel)

#Turn one-sided p-values into two-sided p-values
one.to.two <- function(p)
  {
    2*min(p,1-p)
  }

#group is a vector of group labels
#lib.size is the effective library size

estDisp <- function(counts,group,lib.size)
  {
    estimateCommonDisp(DGEList(counts=counts,group=group,lib.size=lib.size))$common.dispersion
  }    

rpDisp <- function(x,y,nbins=20,min.bin=10,trim.x=0.1,fixdisp=0.2,fixalpha=0.01)
  {
#Sort for later    
    order.x <- order(x)
    x <- x[order.x]
    y <- y[order.x]
#Lower quantile of x
    x.lower <- quantile(x,trim.x)
#Upper quantile of x
    x.upper <- quantile(x,1-trim.x)
    keepers <- which(x>=x.lower & x<=x.upper)
    x <- x[keepers]
    y <- y[keepers]
#Now fit by ls
    fit <- lsfit(x,y,intercept=FALSE)
#    fit <- rlm(x,y)    
    fitted.values <- x*fit$coefficients
#Do two runs, one with a guess of alpha, the other with the estimated alpha
    j <- 1
    while(j<3)
      {
        if(j==1) pvalues <- pnbinom(y,mu=fitted.values,size=1/fixdisp)
        else if (j==2) pvalues <- pnbinom(y,mu=fitted.values,size=1/disp)
        pvalues <- sapply(pvalues,one.to.two)
        keepers <- rep(TRUE,length(pvalues))
        keepers[which(pvalues<fixalpha)] <- FALSE
        x.keepers <- x[keepers]
        y.keepers <- y[keepers]
#    lowessfit <- rlm(x,y,f=f)    
#Starting and ending of bins based on quantiles so same number in every bin
        x.bins <- quantile(x.keepers,seq(0,1,length.out=nbins+1))
#    x.bins <- quantile(x[x>=x.lower & x<=x.upper],seq(0,1,length.out=nbins+1))    
        v <- rep(NA_real_,nbins)
        lambda <- rep(NA_real_,nbins)
        count.bins <- rep(NA_integer_,nbins)
        for(i in seq_len(nbins))
          {
#Which x are in the current bin
            which.bin <- which(x.keepers>=x.bins[i] & x.keepers<x.bins[i+1])
#How many are in the current bin
            count.bins[i] <- length(which.bin)
#Unnecessary check for enough in bin
            if(length(which.bin)>=min.bin)
              {
#Which ys are in bin            
                current.y <- y.keepers[which.bin]
##Variance of trimmed y value
#            v[i] <- mad(current.y)^2
                v[i] <- var(current.y)
#            v[i] <- mad(current.y)^2            
#            v[i] <- var(trimmed.y)            
#Which value of sorted x is in middle of bin
                x.middle <- median(x.keepers[which.bin])
#            lowessind <- findIndices(lowessfit$x,(x.bins[i]+x.bins[i+1])/2)            
#What is the y value at the middle bin position, that is the intensity
                lambda[i] <- x.middle*fit$coefficients
#            lambda[i] <- lowessfit$y[lowessind]            
              }
          }
        lambda <- lambda[!is.na(lambda)]
        v <- v[!is.na(v)]
        term1 <- sum(v*lambda*(1+lambda))
        term2 <- sum(lambda^3)
        term3 <- sum(lambda^4)
        disp <- (term1-term2)/term3
        j <- j+1
      }
    list(disp=disp,lambda=lambda,v=v,x.bins=x.bins,count.bins=count.bins,pvalues=sort(pvalues),fitted.values=fitted.values,keepers=keepers,x=x,y=y)
  }

fitter <- function(rna.vector, rp.vector, method, trim.x=NULL, trim.y=FALSE, disp=0.2, alpha=0.001)
  {
    if(method=="rlm") fit <- lsfit(rna.vector,rp.vector)
#    if(method=="rlm") fit <- rlm(rna.vector,rp.vector)
    else
      {
        if(!is.null(trim.x))
          {
            quantile.low <- quantile(rna.vector,trim.x)
            quantile.high <- quantile(rna.vector,1-trim.x)
            keepers.x <- which(rna.vector>=quantile.low & rna.vector<=quantile.high)
            rna.vector <- rna.vector[keepers.x]
            rp.vector <- rp.vector[keepers.x]
          }
        if(trim.y)
          {
            fit1 <- lsfit(rna.vector,rp.vector,intercept=FALSE)
            fitted.values <- rna.vector*fit1$coefficients
            pvalue.vector <- pnbinom(rp.vector,mu=fitted.values,size=1/disp)
            keepers.y <- which(pvalue.vector>alpha & pvalue.vector<(1-alpha))
            rna.vector <- rna.vector[keepers.y]
            rp.vector <- rp.vector[keepers.y]
          }
        if(method=="ls")
          {
            fit <- lsfit(rna.vector,rp.vector,intercept=FALSE)
          }
        else if(method=="blue")
          {
            fit1 <- lsfit(rna.vector,rp.vector,intercept=FALSE)
            fit1.values <- rna.vector*fit1$coefficients
            fit1.phi <- rpDisp(rna.vector,rp.vector)$disp
            fit1.variances <- fit1.values*(1+fit1.values*fit1.phi)
            fit <- lsfit(rna.vector,rp.vector,wt=1/fit1.variances,intercept=FALSE)
          }
      }
    fit
  }

combined3.p <- function(vec)
  {
    n <- length(vec)
    n.fac <- gamma(n+1)
    sum.vec <- sum(vec)
    floor.vec <- floor(sum.vec)
    output <- 0
    for(k in 0:floor.vec)
      {
        output <- output+((-1)^k)*(n.fac)*((sum.vec-k)^n)/(gamma(n-k+1)*gamma(k+
1))
#        print(output)
    }
    output <- output/n.fac
    output
  }

doBetween <- function(pmat,group,nSD=3,type=c("one-sided","two-sided"))
  {
    zmat <- qnorm(as.matrix(pmat))
    group <- unlist(group)
    unique.group <- unique(group)
    n <- length(unique.group)
    m <- choose(n,2)
    pval <- matrix(NA_real_,nrow(pmat),m)
    direction <- matrix(1,nrow(pmat),m)
    colnames(pval) <- 1:m
    counter <- 1
    for(i in seq_len(n-1))
      {
        which.i <- which(group==unique.group[i])
        zmat.i <- as.matrix(zmat[,which.i])
        num1 <- rowSums(zmat.i)
        for(j in (i+1):n)
          {
            which.j <- which(group==unique.group[j])
            zmat.j <- as.matrix(zmat[,which.j])
            num2 <- rowSums(zmat.j)
            num <- num1-num2
            p <- length(which.i)+length(which.j)
            which.del <- which(abs(num)>(nSD*sqrt(p)))
            zmat.ij <- cbind(zmat.i,zmat.j)
            if(length(which.del)>0) zmat.ij <- zmat.ij[-which.del,]
            covmat <- cov(zmat.ij)
            var.ij <- sum(diag(covmat))
            for(ii in seq_len(ncol(zmat.ij)-1))
                {
                  for(jj in (ii+1):ncol(zmat.ij))
                    {
                      if(ii<=length(which.i) & jj>length(which.i))
                        {
                          var.ij <- var.ij-2*covmat[ii,jj]
                        }
                      else
                        {
                          var.ij <- var.ij+2*covmat[ii,jj]
                        }
                    }
                }
            se.ij <- sqrt(var.ij)
            tstat <- num/se.ij
            if(type=="one-sided") pval[,counter] <- pnorm(tstat)
            else if (type=="two-sided")
              {
                pval[,counter] <- 2*(1-pnorm(abs(tstat)))
                which.switch <- which(tstat>0)
                direction[which.switch,counter] <- (-1)
              }
            colnames(pval)[counter] <- paste(unique.group[i],".vs.",unique.group[j],sep="")
            counter <- counter+1
          }
      }
    list(pval=pval,direction=direction)
  }

doCount <- function(rnv.pos,rnadisp,fit,rpdisp,match.rnv.urnv,rpv)
  {
    length.pos <- length(rnv.pos)
    len <- 10000*length.pos
    nbsample <- rnbinom(len,mu=rnv.pos,size=1/rnadisp)
    mu <- (nbsample*fit$coefficients)
    ref <- (matrix(rnbinom(len,mu=mu,size=1/rpdisp),nrow=length.pos,ncol=10000))[match.rnv.urnv,]
    count <- rowSums(ref>rpv)
    count
  }
    
doWithin <- function(rna,rp,group,nreps,keeper.genes,trim.x=0.1,trim.var=0.1,rnadisp=NULL)
{
  n <- ncol(rna)
  rna <- rna[keeper.genes,]
  rp <- rp[keeper.genes,]
  smat <- matrix(0,nrow(rna),n)
  rownames(smat) <- rownames(rna)
  colnames(smat) <- colnames(rna)
  if(!is.null(rnadisp)) rnadisp <- rnadisp
  else
    {
      unique.group <- unique(group)
      count.group <- rep(NA_integer_,length(unique.group))
      for(i in seq_len(length(unique.group))) count.group[i] <- length(which(group==unique.group[i]))
      keeper.groups <- which(count.group>1)
      if(length(keeper.groups)<2)
         {
           warning("Must be at least two sanmples in at least two groups to estimate rna dispersion.  Fixing at 0.1.  May want to select value for rnadisp in argument to babel")
           rnadisp <- 0.1
         }
      else
         {
           keeper.samples <- which(!is.na(match(group,unique.group[keeper.groups])))
           rnadisp <- estDisp(rna[,keeper.samples],group[keeper.samples],NULL)
         }
    }
#
# Permute
#
  loops <- nreps%/%10000
  for(i in seq_len(n))
    {
      print(paste("Running",colnames(rna)[i]))
      rnv <- rna[,i]; # This sample's RNA
      rpv <- rp[,i]; # This sample's Ribosome
      rpdisp <- rpDisp(rnv,rpv)$disp; # Defaults for rpDisp currently match our desired input
      urnv <- unique(rnv); # Unique values of RNA read counts
      match.rnv.urnv <- match(rnv,urnv)
      pos <- match(urnv,rnv); # Position of unique RNA read count values
      fit <- fitter(rnv,rpv,method="ls",trim.x=trim.x,trim.y=FALSE)
      rnv.pos <- rnv[pos]
      counts <- mclapply(1:loops, FUN=function(j) {
        doCount(rnv.pos=rnv.pos, rnadisp=rnadisp, fit=fit,                   
                rpdisp=rpdisp, match.rnv.urnv=match.rnv.urnv, rpv=rpv)
      })
          count <- Reduce(`+`,counts)
          smat[,i] <- count
    }
  pmat <- (smat+1)/(nreps+2)
  pmat
}

doCombined <- function(pmat,group)
{
  unique.group <- unique(group)
  n <- length(unique.group)
  cmat <- matrix(NA_real_,nrow(pmat),n)
  for(i in seq_len(n))
    {
      which.i <- which(group==unique.group[i])
      p.i <- as.matrix(pmat[,which.i])
      cmat[,i] <- apply(p.i,1,combined3.p)
    }
  colnames(cmat) <- paste("combined.",unique.group,sep="")
  cmat
}

formatWithin <- function(within,method,p,rnames,cnames,keeper.genes,n)
  {
    twosided <- 2*within
    alt <- 2*(1-within)
    which.switch <- which(alt<twosided)
    twosided[which.switch] <- alt[which.switch]
    qvalues <- apply(twosided,2,p.adjust,method=method)
    direction <- matrix(1,nrow(within),ncol(within))
    direction[which(within>0.5)] <- (-1)
    output.within <- vector("list",p)
    names(output.within) <- cnames
    for(i in seq_len(p))
      {
        new.direction <- new.within <- new.twosided <- new.qvalues <- rep(NA_real_,n)
        new.direction[keeper.genes] <- direction[,i]
        new.within[keeper.genes] <- within[,i]
        new.twosided[keeper.genes] <- twosided[,i]
        new.qvalues[keeper.genes] <- qvalues[,i]
        output.within[[i]] <- cbind.data.frame(rnames,new.direction,new.within,new.twosided,new.qvalues)
#        output.within[[i]] <- cbind.data.frame(rnames,direction[,i],within[,i],twosided[,i],qvalues[,i])        
        rownames(output.within[[i]]) <- NULL
        colnames(output.within[[i]]) <- c("Gene","Direction","P-value (one-sided)","P-value (two-sided)","FDR")        
      }
    output.within
  }

formatCombined <- function(combined,group,method,rnames,keeper.genes,n)
{
  twosided <- 2*combined
  alt <- 2*(1-combined)
  which.switch <- which(alt<twosided)
  twosided[which.switch] <- alt[which.switch]
  qvalues <- apply(twosided,2,p.adjust,method=method)
  direction <- matrix(1,nrow(combined),ncol(combined))
  direction[which(combined>0.5)] <- (-1)
  ugroup <- unique(group)
  lgroup <- length(ugroup)
  output.combined <- vector("list",lgroup)
  names(output.combined) <- ugroup
  for(i in seq_len(lgroup))
    {
        new.direction <- new.twosided <- new.qvalues <- rep(NA_real_,n)
        new.direction[keeper.genes] <- direction[,i]
        new.twosided[keeper.genes] <- twosided[,i]
        new.qvalues[keeper.genes] <- qvalues[,i]
        output.combined[[i]] <- cbind.data.frame(rnames,new.direction,new.twosided,new.qvalues)
#        output.combined[[i]] <- cbind.data.frame(rnames,direction[,i],twosided[,i],qvalues[,i])
        rownames(output.combined[[i]]) <- NULL
        colnames(output.combined[[i]]) <- c("Gene","Direction","P-value","FDR")        
    }
  output.combined
}

formatBetween <- function(between,rna,group,method,keeper.genes=keeper.genes,n=n)
  {
    min.group <- min(table(group))
    ugroup <- unique(group)
    lgroup <- length(ugroup)
    m <- choose(lgroup,2)
    output.between <- vector("list",m)
    counter <- 1
    for(i in seq_len(lgroup-1))
      {
        for(j in (i+1):lgroup)
          {
            names(output.between)[counter] <- paste(ugroup[i],".vs.",ugroup[j],sep="")
            qvalues <- p.adjust(between$pval[,counter],method=method)
            if(min.group<2)
              {
                new.between <- new.qvalues <- new.direction <- rep(NA_real_,n)
                new.between[keeper.genes] <- between$pval[,counter]
                new.qvalues[keeper.genes] <- qvalues
                new.direction[keeper.genes] <- between$direction[,counter]
                output.between[[counter]] <-cbind.data.frame(rownames(rna),new.between,new.qvalues,new.direction)
#                output.between[[counter]] <-cbind.data.frame(rownames(rna),between$pval[,counter],qvalues,between$direction[,counter])
                colnames(output.between[[counter]]) <- c("Gene","P-value","FDR","Direction")
              }
            else
              {
                which.ij <- which(group==ugroup[i]|group==ugroup[j])
                dge <- suppressMessages(DGEList(counts=rna[,which.ij],group=group[which.ij]))
                dsp <- estimateCommonDisp(dge)
                dsp <- estimateTagwiseDisp(dsp)
                res <- topTags(exactTest(dsp,pair=c(ugroup[j],ugroup[i])),n=nrow(rna))$table
                match.genes <- match(rownames(rna),rownames(res))
                mrna.qval <- res$FDR[match.genes]
                mrna.logfc <- res$logFC[match.genes]
                change.type <- rep("both",nrow(rna))
                change.type[mrna.qval>0.05 & abs(mrna.logfc)<1.5] <- "translational_only"
                new.between <- new.qvalues <- new.direction <- rep(NA_real_,n)
                new.between[keeper.genes] <- between$pval[,counter]
                new.qvalues[keeper.genes] <- qvalues
                new.direction[keeper.genes] <- between$direction[,counter]
                output.between[[counter]] <-cbind.data.frame(rownames(rna),mrna.logfc,mrna.qval,change.type,new.between,new.qvalues,new.direction)
#                output.between[[counter]] <-cbind.data.frame(rownames(rna),mrna.logfc,mrna.qval,change.type,between$pval[,counter],qvalues,between$direction[,counter])
                colnames(output.between[[counter]]) <- c("Gene","mRNA_logFC","mRNA_FDR","Change_type","P-value","FDR","Direction")
              }
            counter <- counter+1
          }
      }
    output.between
  }
                          
babel <- function(rna,rp,group,nreps,method.adjust="BH",min.rna=10,...)
  {
    if(sum(dim(rna)==dim(rp))<2) stop("rna and rp are different sizes")
    n <- nrow(rna)
    p <- ncol(rp)
    if(length(group)!=p) stop("group must have same length as the number of columns of rna")
    if(is.null(rownames(rna))) rownames(rna) <- 1:n
    if(is.null(rownames(rp))) rownames(rp) <- 1:n
    if(is.null(colnames(rna))) colnames(rna) <- 1:p
    if(is.null(colnames(rp))) colnames(rp) <- 1:p
    if(sum(rownames(rna)==rownames(rp))<n) stop("rownames of rna and rp must match")
    if(sum(colnames(rna)==colnames(rp))<p) stop("colnames of rna and rp must match")
    if((nreps%%10000)!=0) stop("nreps must be divisible by 10000")
    if(nreps<100000) stop("nreps must at least 100,000")
    if(min.rna<1) stop("min.rna needs to be at least 1")
#    if(length(unique(group))!=2) stop("There must be exactly two groups")
    mins.rna <- apply(rna,1,min)
    keeper.genes <- which(mins.rna>=min.rna)
    within <- doWithin(rna=rna,rp=rp,group=group,nreps=nreps,keeper.genes=keeper.genes,...)
    output.within <- formatWithin(within,method=method.adjust,p=p,rnames=rownames(rna),cnames=colnames(rna),keeper.genes=keeper.genes,n=n)
    combined <- doCombined(within,group)
    output.combined <- formatCombined(combined,group=group,method=method.adjust,rnames=rownames(rna),keeper.genes=keeper.genes,n=n)
    between <- doBetween(within,group,type="two-sided")
    output.between <- formatBetween(between,rna=rna,group=group,method=method.adjust,keeper.genes=keeper.genes,n=n)
    list(within=output.within,combined=output.combined,between=output.between)
  }

