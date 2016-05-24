vario <- function (n.bins=20, size.bins=NULL, extent=0.5, data, data2=NULL, 
                   is.latlon=TRUE, is.centered=FALSE, nrands=0,
                   type=c("semivar", "cov", 
                          "pearson", "spearman", "kendall", "moran", "geary"),
                   alternative=c("one.tailed", "two.tailed"),
                   mult.test.corr=c("none", "holm", "hochberg", "bonferroni"),
                  quiet = FALSE) {
  
  tails=c("one.tailed", "two.tailed")
  alternative=match.arg(tolower(alternative), tails)
  
  types=c("semivar", "cov", "pearson", "spearman", "kendall", "moran", "geary")
  type=match.arg(tolower(type), types)
  
  mults=c("none", "holm", "hochberg", "bonferroni")
  mult.test.corr=match.arg(tolower(mult.test.corr), mults)
  
  n.cols=NCOL(data)
  if (n.cols > 3)
    is.multivar=TRUE
  else
    is.multivar=FALSE
  
  if (!is.null(size.bins))
    n.bins=NULL
  
  results=vario.aux (n.bins=n.bins, size.bins=size.bins, extent=extent, data=data, data2=data2, 
                     is.latlon=is.latlon, is.centered=is.centered, 
                     is.multivar=is.multivar, type=type)
  
  if (nrands > 0) {
    rands=matrix(NA, nrow=nrands+1, ncol=length(results$bins))
    if (!quiet)
      prog.bar=txtProgressBar(min = 0, max = nrands, style = 3)
    
    if (is.null(data2))
      data2=data
    for (i in 1:nrands) {
      s=sample(results$grpdata)
      
      if (is.multivar) {
        vals=results$vals
        rands[i,]=tapply(vals, s, FUN=mean, na.rm=TRUE)
      }
      else {
        for (j in 1:(length(results$bins))) {        
          tmp=results$all.combs[s==j,]
          tmp=tmp[complete.cases(tmp),]
          rands[i,j]=vario.func(data[tmp[,1], 3], data2[tmp[,2], 3], results$glob.mean, 
                                results$glob.sd, results$glob.N, is.multivar, type)
        }
      }
      if (!quiet)
        setTxtProgressBar(prog.bar, i)
    }
    rands=rands-rowMeans(rands, na.rm=TRUE)
    crit.val=0
    if (!is.centered) {
      rands=rands+results$regional.mean
      crit.val=results$regional.mean
    }
    rands[nrands+1,]=results$vario
    
    if (alternative == "two.tailed") {
      pvals=apply(rands, MARGIN=2, 
                  function (x) {
                    sum(abs(x - crit.val) >= abs(x[nrands+1] - crit.val))/(nrands+1)
                  })
    }
    else {
      pvals=apply(rands, MARGIN=2, function (x) {
        ifelse (x[nrands+1] > crit.val, 
                sum(x >= x[nrands+1])/(nrands+1), 
                sum(x <= x[nrands+1])/(nrands+1))})
    }

    if (mult.test.corr != "none") {
      pvals=p.adjust(pvals, method=mult.test.corr[1])
    }
    
    colnames(rands)=names(results$bins)
    names(pvals)=names(results$bins)
    results$pvals=pvals
    results$rands=rands
    results$alternative=alternative
    results$mult.test.corr=mult.test.corr[1]
  }
  
  ## Remove extraneous elements
  results[8:13]=NULL
  results$is.multivar=is.multivar
  
  class(results)="vario"
  return(results)
}

vario.aux <- function (n.bins=20, size.bins=NULL, extent=0.5, data, data2=NULL, is.latlon=TRUE, 
                       is.centered=FALSE, is.multivar=FALSE,
                       type=c("semivar", "cov", "pearson", "spearman", "kendall", 
                              "moran", "geary")) {
  
  n.cols=NCOL(data)
  all.dists=coord2dist(data[, 1:2], is.latlon)
  ## Compute maximum distance
  max.dist=max(all.dists)
  max.extent=max.dist*extent
  
  if (!is.null(data2))
    include.lag0=TRUE
  else
    include.lag0=FALSE
  
  types=c("semivar", "cov", "pearson", "spearman", "kendall", "moran", "geary")
  type=match.arg(tolower(type), types)
  if (is.null(size.bins)) {
    bins=seq(0, max.extent, length.out=n.bins+1)
    grpdata <-cut(all.dists, breaks=bins, labels=1:(length(bins)-1), right=TRUE)
  }
  else {
    bins=seq(0, max.extent+size.bins, by=size.bins)
    grpdata <-cut(all.dists, breaks=bins, labels=1:(length(bins)-1), right=FALSE)
  }
  if (is.multivar) {
    glob.mean=NA
    glob.sd=NA
    glob.N=NA
    all.combs=NA
    
    if (is.null(data2)) {
      if (type=="cov")
        vals=cov(t(data[, 3:n.cols]))
      else {
        vals=suppressWarnings(cor(t(data[, 3:n.cols]), method=type, use = "pairwise.complete.obs"))
      }
      vals=vals[lower.tri(vals)]
    }
    else {
      if (type=="cov")
        vals=suppressWarnings(cov(x=t(data[, 3:n.cols]), y=t(data2[, 3:n.cols]), use = "pairwise.complete.obs"))
      else
        vals=suppressWarnings(cor(x=t(data[, 3:n.cols]), y=t(data2[, 3:n.cols]), method=type, use = "pairwise.complete.obs"))
      # vals=vals[row(vals)!=col(vals)]
      vals=vals[lower.tri(vals)]
    }
    
    regional.mean=mean(vals, na.rm=TRUE)
    if (is.centered) {
      vals=vals-regional.mean
    }
    vario=tapply(vals, grpdata, mean, na.rm=TRUE)
    npoints=tapply(vals, grpdata, FUN=function (x) {length(na.omit(x))})
    bin.dist=tapply(all.dists, grpdata, FUN=mean, na.rm=TRUE)
  }
  else {
    vals=NA
    bin.dist=numeric(length(bins)-1)*NA
    vario=numeric(length(bins)-1)*NA
    npoints=numeric(length(bins)-1)*NA
    
    all.combs=t(combn(NROW(data), 2))
    if (!is.null(data2)) {
      glob.mean=c(mean(data[,3], na.rm=TRUE), mean(data2[,3], na.rm=TRUE))
      glob.sd=c(sd(data[,3], na.rm=TRUE), sd(data2[,3], na.rm=TRUE))
      glob.N=NROW(data[,3])
      all.combs=cbind(all.combs, all.combs[, c(2, 1)])
      ## Control for the fact that including lag0 means pairs of sites are
      ## counted twice
      denom.N=2
    }
    else {
      data2=data
      glob.mean=rep(mean(data[,3], na.rm=TRUE), 2)
      glob.sd=rep(sd(data[,3], na.rm=TRUE), 2)
      glob.N=NROW(data[,3])
      denom.N=1  
    }
    for (i in 1:(length(bins)-1)) {
      if (include.lag0) {
        tmp=rbind(all.combs[grpdata==i, 1:2], all.combs[grpdata==i, c(3, 4)])
      }
      else
        tmp=all.combs[grpdata==i, 1:2]
      tmp=tmp[complete.cases(tmp),]
      # produce NAs when there are no pairs of points in a bin
      if (!is.null(dim(tmp)) & length(tmp)>0) {
        x=data[tmp[,1], 3:n.cols]
        y=data2[tmp[,2], 3:n.cols]    
        npoints[i]=NROW(x)/denom.N
        vario[i]=vario.func(x, y, glob.mean, glob.sd, glob.N, is.multivar, type=type)        
        bin.dist[i]=mean(all.dists[grpdata==i], na.rm=T)
      }
      else {
        vario[i]=NA
        npoints[i]=length(tmp)/denom.N
        bin.dist[i]=NA
      }
    }
    regional.mean=mean(vario, na.rm=TRUE)
    if (is.centered)
      vario=vario-regional.mean
  }
  
  bins=bins[1:(length(bins))-1]
  col.names=1:length(bins)
  names(bins)=col.names
  names(bin.dist)=col.names
  names(vario)=col.names
  names(npoints)=col.names
  
  return (list(bins=bins, mean.bin.dist=bin.dist,
               vario=vario, npoints=npoints, metric=type, is.centered=is.centered,
               regional.mean=regional.mean, all.combs=all.combs, grpdata=grpdata, 
               glob.mean=glob.mean, glob.sd=glob.sd, glob.N=glob.N, vals=vals))
}
