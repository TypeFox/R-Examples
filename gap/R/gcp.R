gcp <- function(y,cc,g,handle.miss=1,miss.val=0,n.sim=0,locus.label=NULL,quietly=FALSE)
{
# autosomal data
  tmp <- geno.recode(g,miss.val=miss.val)
  g <- tmp$grec
  nloci <- dim(g)[2]/2
  loci <- rep(0,nloci)
  for(i in 1:nloci) loci[i] <- length(tmp$alist[[i]]$allele)
  g[is.na(g)] <- 0
  control <- gc.control(handle.miss=handle.miss,verbose=F)
# marker-marker association
  if (cc==0)
  {
    x <- 0
    n.obs <- length(g[,1])
    n.loci <- length(g[1,])/2
    g <- g
    g.gc <- genecounting(g,loci=loci,control=control)
    n.haps <- length(g.gc$h)
    x2obs <- 2*(g.gc$l1-g.gc$l0)
    pobs <- 1 - pchisq(x2obs, n.haps-1)
    zobs <- z <- ph <- rep(0,n.haps)
    h0 <- g.gc$h0
    h1 <- g.gc$h
    m0 <- (h0 > 0)
    lh <- length(m0)
    e <- 2.0 * h0 * n.obs
    o <- 2.0 * h1 * n.obs
    zobs[m0] <- sqrt(o) + sqrt(o + 1) - sqrt(4 * e + 1)
    psim <- 0
    for (i in 1:n.sim)
    {
    # all association
      for (j in 1:(n.loci-1))
      {
          rand.ord <- order(runif(n.obs))
          g[,(2*j-1)] <- g[,(2*j-1)][rand.ord]
          g[,(2*j)] <- g[,(2*j)][rand.ord]
      }
      g.gc <- genecounting(g,loci=loci,control=control)
      x2sim <- 2*(g.gc$l1-g.gc$l0)
      if (x2sim >= x2obs) x <- x + 1
      h0 <- g.gc$h0
      h1 <- g.gc$h
      m1 <- (h0 > 0)
      e <- 2.0 * h0 * n.obs
      o <- 2.0 * h1 * n.obs
      z[m1] <- sqrt(o[m1]) + sqrt(o[m1] + 1) - sqrt(4 * e[m1] + 1)
      m2 <- (abs(z[m1]) >= abs(zobs[m1]))
      ph[m2] <- ph[m2] + 1
    }
  } else {
#   marker-disease association
    is.case <- (y==1)
    caco.obs <- genecounting(g,loci=loci,control=control)
    ca.obs <- genecounting(g[is.case,],loci=loci,control=control)
    co.obs <- genecounting(g[!is.case,],loci=loci,control=control)
    yobs <- y
    n.haps <- length(caco.obs$h)
    x2obs <- 2*(ca.obs$l1+co.obs$l1-caco.obs$l1)
    pobs <- 1 - pchisq(x2obs, n.haps-1)
    zobs <- z <- ph <- rep(0,n.haps)
    n.caco <- length(y)
    n.ca <- length(y[y==1])
    n.co <- n.caco - n.ca
    h0 <- co.obs$h
    h1 <- ca.obs$h
    m <- (n.ca*h1+n.co*h0)/(n.ca+n.co)
    m0 <- (m > 0)
    lh <- length(m[m0])
    zobs[m0] <- (h1[m0]-h0[m0])/sqrt((1/n.ca+1/n.co)*m[m0]*(1-m[m0]))
    psim <- 0
    for (i in 1:n.sim)
    {
       y <- yobs[order(runif(n.ca+n.co))]
       is.case <- (y==1)
       ca.j <- genecounting(g[is.case,],loci=loci,control=control)
       co.j <- genecounting(g[!is.case,],loci=loci,control=control)
       x2 <- 2*(ca.j$l1+co.j$l1-caco.obs$l1)
       if (x2 >= x2obs) psim <- psim + 1
       h0 <- co.j$h
       h1 <- ca.j$h
       mh <- (n.ca*h1+n.co*h0)/(n.ca+n.co)
       m1 <- ((m > 0) & (mh > 0))
       z[m1] <- (h1[m1]-h0[m1])/sqrt((1/n.ca+1/n.co)*mh[m1]*(1-mh[m1]))
       m2 <- (abs(z[m1]) >= abs(zobs[m1]))
       ph[m2] <- ph[m2] + 1
    }
  }
  m <- (1:n.haps)[m0]
  hap <- revhap(loci,m)
  t1 <- zobs[m0]
  if (is.null(locus.label)) locus.label <- paste("locus-",1:nloci,sep="")
  if (n.sim > 0)
  {
    if (!quietly)
    {
      cat("\nLRT = ",round(x2obs,3), ", p = ",pobs, ", sim p = ",psim/n.sim, "\n")
      cat("\nz-test of individual haplotypes\n")
      t2 <- ph[m0]/n.sim
      tbl <- data.frame(m,round(t1,3),round(t2,4),grec2g(hap,nloci,tmp))
      names(tbl) <- c("hap.id","z","sim p",locus.label)
      print(tbl, quote=FALSE)
    }
  # likelihood-based LD statistics to be added
    invisible(list(x2obs=x2obs,pobs=pobs,psim=psim/n.sim,zobs=zobs,ph=ph/n.sim))
  }
  else
  {
    if (!quietly)
    {
      cat("\nLRT = ",round(x2obs,3), ", p = ",pobs, "\n")
      cat("\nz-test of individual haplotypes\n")
      tbl <- data.frame(m,round(t1,3),grec2g(hap,nloci,tmp))
      names(tbl) <- c("hap.id","z",locus.label)
      print(tbl, quote=FALSE)
    }
    invisible(list(x2obs=x2obs,pobs=pobs,zobs=zobs))
  }
}

# 1-8-2004 start implementation
# 24-9-2004 incorporate code from haplo.score
