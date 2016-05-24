allele.recode <- function (a1, a2, miss.val = NA)
{
    n <- length(a1)
    if (is.factor(a1))
        a1 <- as.character(a1)
    if (is.factor(a2))
        a2 <- as.character(a2)
    is.ch <- is.character(a1) | is.character(a2)
    if (is.ch) {
        t <- factor(c(a1, a2), exclude = miss.val)
    }
    if (!is.ch) {
        lev <- sort(unique(c(a1, a2)))
        t <- factor(c(a1, a2), levels = lev, exclude = miss.val)
    }
    allele.label <- levels(t)
    t <- as.numeric(t)
    a1 <- t[1:n]
    a2 <- t[(n + 1):(2 * n)]
    return(list(a1 = a1, a2 = a2, allele.label = allele.label))
}

geno.recode <- function (geno, miss.val = 0)
{
    n.loci <- ncol(geno)/2
    alist <- vector("list", n.loci)
    grec <- NULL
    for (i in 1:n.loci) {
        t <- (i - 1) * 2 + 1
        tmp <- allele.recode(geno[, t], geno[, t + 1], miss.val = miss.val)
        grec <- cbind(grec, tmp$a1, tmp$a2)
        alist[[i]] <- list(allele = tmp$allele.label)
    }
    return(list(grec = grec, alist = alist))
}

a2g <- function(a1,a2)
{
  i <- ifelse(a1 < a2,a2,a1)
  j <- ifelse(a1 < a2,a1,a2)
  genocode <- ifelse (j==0, 0, i*(i-1)/2 + j)
  invisible (genocode)

}

g2a.c <- function (g)
{
    d <- 1 + 8 * (g - 1)
    u <- 1 + ((1 + (sqrt(d) - 1) - 1) / 2)
    u <- ceiling(u)
    l = g - u * (u - 1) / 2
    list (l=l,u=u)
}

g2a <- function (g)
{
    i <- 1 + floor((sqrt(8 * g + 1) - 1)/2)
    j <- g - i * (i - 1)/2
    i <- ifelse(j == 0, i - 1, i)
    j <- ifelse(j == 0, i, j)
    return(cbind(j,i))
}

is.miss <- function(data,data.int,miss.val=0)
{
   if (data.int==2) # genotype
   {
      id <- array(FALSE, length(data))
      for (i in 1:length(miss.val))
          id <- id | data==miss.val[i]
   }
   else # allele
   {
      id <- array(FALSE, length(data[,1]))  
      for (i in 1:length(miss.val))
          id <- id | apply(data==miss.val[i],1,any)
   }
   return (id)  
}

# adapted from gcp.c on 24/9/2004
# 17-6-2004 JH Zhao

revhap <- function(loci,hapid)
{
   nloci <- length(loci)
   nalleles <- vector("numeric",nloci)
   nalleles[nloci] <- 1
   m <- nloci
   for(k in nloci:1)
   {
      m <- m - 1
      nalleles[m] <- nalleles[m+1] * loci[k]
   }
  n <- length(hapid)
  hap <- matrix(0,n,nloci)
  for (i in 1:n)
  {
    l <- hapid[i] - 1
    for (j in 1:nloci)
    {
      hap[i,j] <- floor(l/nalleles[j])
      if (j==nloci) hap[i,j] <- l
      else l <- l%%nalleles[j]
      hap[i,j] <- hap[i,j] + 1
    }
  }
  invisible(hap)

}

revhap.i <- function(loci,hapid)
{
   nloci <- length(loci)
   nalleles <- vector("numeric",nloci)
   nalleles[nloci] <- 1
   m <- nloci
   for(k in nloci:1)
   {
      m <- m - 1
      nalleles[m] <- nalleles[m+1] * loci[k]
   }

  hap <- vector("numeric",nloci)
  l <- hapid - 1
  for (j in 1:nloci)
  {
    hap[j] <- floor(l/nalleles[j])
    if (j==nloci) hap[j] <- l
    else l <- l%%nalleles[j]
    hap[j] <- hap[j] + 1
  }
  invisible(hap)

}

gcode <- function(a1,a2) {

  i <- ifelse(a1 < a2,a2,a1)
  j <- ifelse(a1 < a2,a1,a2)
  genocode <- i*(i-1)/2 + j
  return(genocode)

}

ungcode <- function(g) {

  i <- 1 + floor((sqrt(8*g+1)-1)/2)
  j <- g - i*(i-1)/2
  
  # the following 2 lines were added as a patch to make this ungcode work:
  i <- ifelse(j==0,i-1,i)
  j <- ifelse(j==0,i,j)

  return(cbind(j,i))
}

grec2g <- function (h, n, t)
{
  hh <- h
  for (i in 1:n)
  {
    hh[,i] <- t$alist[[i]]$allele[h[,i]]
  }
  invisible(hh)
}

# a simple scheme to represent SNPs
# similar to a2g() and not unlike Mike Weale's Matlab function
# 10-2-2005 JHZ
m2plem <- function(a1,a2)
{
  a <- vector("numeric")
# 1=A, 2=a
  for (i in 1:length(a1))
  {
      if(a1[i]==1&a2[i]==2) {
         a[i] <- 0
      }
      else if (a1[i]==1&a2[i]==1) {
         a[i] <- 1
      }
      else if (a1[i]==2&a2[i]==2) {
         a[i] <- 2
      }
  }
  invisible(a)
}

plem2m <- function(a)
{
  a1 <- vector("numeric")
  a2 <- vector("numeric")
# 1=A, 2=a
  for (i in 1:length(a))
  {
      if(a[i]==0) {
         a1[i] <- 1
         a2[i] <- 2
      }
      else if (a[i]==1) {
         a1[i] <- 1
         a2[i] <- 1
      }
      else if (a[i]==2) {
         a1[i] <- 2
         a2[i] <- 2
      }
  }
  invisible(list(a1,a2))
}

micombine <- function (est, std.err, confidence = 0.95)
{
    qstar <- est[[1]]
    for (i in 2:length(est)) {
        qstar <- cbind(qstar, est[[i]])
    }
    qbar <- apply(qstar, 1, mean)
    u <- std.err[[1]]
    for (i in 2:length(std.err)) {
        u <- cbind(u, std.err[[i]])
    }
    if (!is.null(dimnames(qstar)[[1]]))
        dimnames(u)[[1]] <- dimnames(qstar)[[1]]
    u <- u^2
    ubar <- apply(u, 1, mean)
    bm <- apply(qstar, 1, var)
    m <- dim(qstar)[2]
    tm <- ubar + ((1 + (1/m)) * bm)
    rem <- (1 + (1/m)) * bm/ubar
    nu <- (m - 1) * (1 + (1/rem))^2
    alpha <- 1 - (1 - confidence)/2
    low <- qbar - qt(alpha, nu) * sqrt(tm)
    up <- qbar + qt(alpha, nu) * sqrt(tm)
    pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
    fminf <- (rem + 2/(nu + 3))/(rem + 1)
    result <- list(est = qbar, std.err = sqrt(tm), df = nu, signif = pval,
        lower = low, upper = up, r = rem, fminf = fminf)
    result
}

VR <- function(v1,vv1,v2,vv2,c12)
{
  nv <- v2^2*vv1 + v1^2*vv2 - 2*v1*v2*c12
  dv <- v2^4
  nv/dv
}

h2G <- function(V,VCOV,verbose=TRUE)
{
  VG <- V[1]
  Ve <- V[2]
  Vp <- VG + Ve
  VVG <- VCOV[1,1]
  VVe <- VCOV[2,2]
  cVGVe <- VCOV[2,1]
  h2G <- VG / Vp
  VVp <- VVG + VVe + 2 * cVGVe
  cVpVG <- VVG + cVGVe
  Varh2G <- VR(VG,VVG,Vp,VVp,cVpVG)
  if (verbose) {
     cat("Vp =",Vp,"SE =",sqrt(VVp),"\n")
     cat("h2G =",h2G,"SE =",sqrt(Varh2G),"\n")
  }
  invisible(list(Vp=Vp,VVp=VVp,h2G=h2G,Varh2G=Varh2G))
}

h2GE <- function(V,VCOV,verbose=TRUE)
{
  VG <- V[1]
  VGE <- V[2]
  Ve <- V[3]
  Vp <- VG + VGE + Ve
  VVG <- VCOV[1,1]
  VVGE <- VCOV[2,2]
  VVe <- VCOV[3,3]
  cVGVGE <- VCOV[2,1]
  cVGVe <- VCOV[3,1]
  cVGEVe <- VCOV[3,2]
  s <- VVG + VVGE + VVe
  s12 <- 2*cVGVGE
  s13 <- 2*cVGVe
  s23 <- 2*cVGEVe
  VVp <- s + s12 + s13 + s23
  cVpVG <- VVG + cVGVGE + cVGVe
  cVpVGE <- cVGVGE + VVGE + cVGEVe
  h2G <- VG / Vp
  Varh2G <- VR(VG, VVG, Vp, VVp, cVpVG)
  h2GE <- VGE / Vp
  Varh2GE <- VR(VGE,VVGE,Vp,VVp, cVpVGE)
  if (verbose) {
     cat("Vp =",Vp,"SE =",sqrt(VVp),"\n")
     cat("h2G =",h2G,"SE =",sqrt(Varh2G),"h2GE =",h2GE,"SE =",sqrt(Varh2GE),"\n")
  }
  invisible(list(Vp=Vp,VVp=VVp,h2G=h2G,Varh2G=Varh2G,h2GE=h2GE,Varh2GE=Varh2GE))
}

h2l <- function(K=0.05,P=0.5,h2,se,verbose=TRUE)
{
  x <- qnorm(1-K)
  z <- dnorm(x)
  1/sqrt(2*pi)*exp(-x^2/2)
  fK <- (K*(1-K)/z)^2
  fP <- P*(1-P)
  f <- fK/fP
  h2l <- f*h2
  sel <- f*se
  z2 <- K^2*(1-K)^2/(f*fP)
  x2 <- -log(2*pi*z2)
  if (verbose) {
     cat("K = ", K, "P = ", P, "\n")
     cat("h2 =",h2,"SE =",se,"h2l =",h2l,"SE =",sel,"\n")
  }
  invisible(list(h2=h2,se=se,h2l=h2l,sel=sel,cc=f,z=sqrt(x2)))
}

# R script to read the GRM binary file

ReadGRMBin <- function(prefix, AllN=FALSE, size=4)
{
  BinFileName <- paste(prefix,".grm.bin",sep="")
  NFileName <- paste(prefix,".grm.N.bin",sep="")
  IDFileName <- paste(prefix,".grm.id",sep="")
  id <- read.table(IDFileName)
  n <- dim(id)[1]
  BinFile <- file(BinFileName, "rb");
  grm <- readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  close(BinFile)
  NFile <- file(NFileName, "rb");
  if(AllN) N <- readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  else N <- readBin(NFile, n=1, what=numeric(0), size=size)
  close(NFile)
  i <- sapply(1:n, function(i) i*(i+1)/2)
  GRM <- matrix(NA,n,n)
  GRM[upper.tri(GRM,diag=TRUE)] <- grm
  GRM[lower.tri(GRM)] <- t(GRM)[lower.tri(GRM)]
  invisible(list(grm=grm, id=id, N=N, GRM=GRM))
}

WriteGRMBin <- function(prefix, grm, N, id, size=4)
{
  BinFileName <- paste(prefix,".grm.bin",sep="")
  NFileName <- paste(prefix,".grm.N.bin",sep="")
  IDFileName <- paste(prefix,".grm.id",sep="")
  grm.bin <- file(BinFileName, "wb")
  writeBin(grm,grm.bin,size=size)
  close(grm.bin)
  grm.N.bin <- file(NFileName, "wb")
  writeBin(N,grm.N.bin,size=size)
  close(grm.N.bin)
  write.table(id,IDFileName,col.names=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
}

ReadGRM <- function(prefix=51)
{
  idfile <- paste(prefix,".grm.id",sep="")
  id <- read.delim(idfile,header=FALSE,sep="\t")
  N <- dim(id)[1]
  L <- N * (N + 1) /2
  grmfile <- paste(prefix,".grm.gz",sep="")
  gz <- gzfile(grmfile)
  grm_lines <- readLines(gz)
  close(gz)
  GRM_list <- sapply(grm_lines,strsplit,"\t")
  M <- as.integer(lapply(GRM_list,"[[",3))
  grm <- as.numeric(lapply(GRM_list,"[[",4))
  GRM <- matrix(NA,N,N)
  GRM[upper.tri(GRM,diag=TRUE)] <- grm
  GRM[lower.tri(GRM)] <- t(GRM)[lower.tri(GRM)]
  invisible(list(GRM=GRM,N=M,id=id,grm=grm))
}

WriteGRM <- function(prefix=51,id,N,GRM)
{
  M <- N
  N <- dim(GRM)[1]
  k2l <- matrix(NA,N*(N+1)/2,4)
  L <- 1
  for(i in 1:N)
  {
    for(j in 1:i)
    {
      k2l[L,] <- c(i,j,M[L],GRM[i,j])
      L <- L + 1
    }
  }
  idfile <- paste(prefix,".grm.id",sep="")
  write(t(id),file=idfile,2,sep="\t")
  grmfile <- paste(prefix,".grm.gz",sep="")
  gz <- gzfile(grmfile,"w")
  write.table(k2l,file=gz,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  close(gz)
}

ReadGRMPLINK <- function(prefix,diag=1)
{
  f <- paste(prefix,".genome",sep="")
  g <- read.table(f,header=TRUE,as.is=TRUE)
  L <- dim(g)[1]
  N <- (1+sqrt(1+8*L))/2
  pi_hat <- g["PI_HAT"]
  PIHAT <- matrix(NA,N,N)
  diag(PIHAT) <- diag
  PIHAT[lower.tri(PIHAT)] <- pi_hat[,1]
  PIHAT[upper.tri(PIHAT)] <- t(PIHAT)[upper.tri(PIHAT)]
  idplink <- function (N)
  {
    id <- function(N, diag=TRUE)
    {
      irep <- function(i) rep(i,i)
      iter <- function(i) 1:i
      i <- unlist(lapply(1:N,irep))
      j <- unlist(lapply(lapply(1:N,iter),cbind))
      ij <- cbind(i,j)
      if (diag) ij
      else subset(ij,ij[,1]!=ij[,2])
    }
    ij <- id(N, diag=FALSE)
    ij[order(ij[,2]),2:1]
  }
  fid1 <- g[1,"FID1"]
  iid1 <- g[1,"IID1"]
  fid <- c(fid1,g[1:(N-1),"FID2"])
  iid <- c(iid1,g[1:(N-1),"IID2"])
  ji <- idplink(N)
  idn <- 1:N
  ii <- cbind(i=idn,j=idn)
  ij <- rbind(ji,ii)
  idx <- ij[order(ij[,2]),2:1]
  invisible(list(pihat=PIHAT[lower.tri(PIHAT,diag=TRUE)],PIHAT=PIHAT,id=cbind(fid,iid),idx=idx))
}

WriteGRMSAS <- function(grmlist,outfile="gwas")
{
  for(p in c("foreign")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!require(p, quietly = TRUE, character.only=TRUE))
        warning(paste("WriteGRMSAS needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  with(grmlist,
  {
    parm <- 1
    colnames(GRM) <- row <- paste("COL",1:dim(id)[1],sep="")
    names(id) <- c("idn","id")
    write.dta(data.frame(id,parm,row,GRM),paste(outfile,".dta",sep=""))
  })
}

ReadGRMPCA <- function(prefix)
{
  values <- scan(paste(prefix,".eigenval",sep=""))
  vectors <- read.table(paste(prefix,".eigenvec",sep=""))
  N <- dim(vectors)[1]
  id <- vectors[,1:2]
  vectors <- as.matrix(vectors[,-(1:2)])
  s <- (values>0)
  posGRM <- vectors[,s]%*%diag(values[s])%*%t(vectors[,s])
  s <- (values<0)
  negGRM <- vectors[,s]%*%diag(values[s])%*%t(vectors[,s])
  invisible(list(N=N,posGRM=posGRM,negGRM=negGRM,id=id))
}
