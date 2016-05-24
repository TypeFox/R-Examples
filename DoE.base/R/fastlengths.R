### utility functions for fast WLP calculations in general case

## functions by Hongquan Xu adapted by Ulrike Groemping
Choose <- function(n, k) if(n < k) 0 else choose(n,k)

Kraw<-function(k,x,n,q)
{ # the Krawtchouk polynomials p<-k(x;n,q)
    kraw <-0;
    for(j in seq(0,k,2))
	kraw <- kraw + (q-1)^(k-j)* Choose(x,j)*Choose(n-x,k-j);
  ## made it work for k=0, May 10 13, UG
    if (k>0) for(j in seq(1,k,2))
	kraw <- kraw - (q-1)^(k-j)* Choose(x,j)*Choose(n-x,k-j);
    kraw
}

# ham: return Hamming distance of two codes
ham<-function(c1, c2)  sum(c1 != c2)

levels.no <- function(xx)
{
# return the number of levels for design xx
# apply(xx, 2, function(v) {length(table(v))} )
# changed table to unique, much faster UG 10 May 13
   apply(xx, 2, function(v) {length(unique(v))} )
}

## mixed level stuff
levelmix <- function(xx)
{ # return column sets for levels
  ss <- levels.no(xx)
  tabss  <- table(ss) # names contains the sj, content the nj
  sso <- as.numeric(names(tabss))
  tabss <- c(tabss)
  names(tabss) <- NULL
  colnos <- lapply(sso, function(obj) which(ss==obj))
  names(colnos) <- sso
  list(ss = sso, ns = tabss, colnos=colnos)
  ## nasty design for trying levelmix(L144.2.75.3.3.4.1.6.6.12.1)
}

distDistmix <-function(code, levm)
{
  # codeDist: return the distance distribution of code
  # i.e. the B
  # levm is a list created with levelmix
  code<-as.matrix(code);
  nRow <- nrow(code);
  colnos <- levm$colnos
  ss <- levm$ss
  ns <- levm$ns
  difflevs <- length(ss)
  dh <- matrix(0, 1, difflevs)  ## initialize
  #dists <- rep(NA, choose(nRow, 2));
  dists <- array(0, dim=ns+1)
    dn <- lapply(ns, function(obj) 0:obj)
    names(dn) <- ss
    dimnames(dists) <- dn
  dists[1] <- nRow  ## separate distances
  dist <- rep(0,sum(ns)+1)
  dist[1] <- nRow  ## overall distance
  for (i in 2:nRow){
    for (j in 1:(i-1)){
      for (k in 1:difflevs)
        dh[k] <-  ham(code[i,colnos[[k]]], code[j,colnos[[k]]]);
        dists[dh+1] <- dists[dh+1] + 2
        dist[sum(dh) + 1] <- dist[sum(dh) + 1] + 2
    }}
  return( list(BSep=dists / nRow, B=dist/nRow ));
}

Bprime <- function(dists, nmax=5){
  ## dists is the BSep element from an object created with distDistmix
  N <- sum(dists)
  ns <- dim(dists)-1
  ss <- as.numeric(names(dimnames(dists)))
  difflevs <- length(ss)
  Bprime <- 0*dists   ## create structure
  combs <- as.matrix(expand.grid(lapply(ns, function(obj) 0:obj)))
                     ## index set for the structure
  combsBd <- combs[rowSums(combs)<=nmax, ,drop=FALSE]

  ## can the necessary summands be reduced?
  for (j in 1:nrow(combsBd)){ ## Bprime-entry, i.e. the j from p.1072
    selj <- combsBd[j,,drop=FALSE]
    zwischen <- 0
    for (i in 1:nrow(combs)){ ## the dists entry
    sel <- combs[i,,drop=FALSE]
    hilf <- dists[sel+1]
    if (hilf>0)
    zwischen <- zwischen + hilf*prod(mapply(Kraw,selj,sel,ns,ss))
    }
    Bprime[selj+1]<-zwischen
  }
  Bprime/N
}

dualDistmix<-function(Bprime, nmax=5)
{
  # dual distance distribution
  # obtains A from Bprime
  ns <- dim(Bprime)-1
  difflevs <- length(ns)

  dual <- rep(0,nmax+1);  dual[1] <- 1;
  for (k in 1:nmax){
    jjj <- xsimplex(difflevs, k)
    if (length(jjj) > 1)
    jjj <- jjj[,apply(jjj,2, function(obj) all(obj<=ns,k)), drop=FALSE]
    else jjj <- matrix(jjj,1,1)
    J <- ncol(jjj)
    for (j in 1:J) dual[k+1] <- dual[k+1] + Bprime[t(jjj[,j,drop=FALSE])+1]
  }
  dual
}

GWLP <- function(design, ...) UseMethod("GWLP")
GWLP.design <- function(design, kmax=design.info(design)$nfactors, attrib.out=FALSE, with.blocks = FALSE, digits=NULL, ...){
    if (!"design" %in% class(design)) stop("GWLP.design is for class design objects only")
    if (with.blocks)
    GWLP.default(design[,c(design.info(design)$block.name,names(factor.names(design)))], kmax=kmax, attrib.out=attrib.out, digits=digits, ...)
    else
    GWLP.default(design[,names(factor.names(design))], kmax=kmax, attrib.out=attrib.out, digits=digits, ...)
}
GWLP.default <- function(design, kmax=ncol(design), attrib.out=FALSE, digits=NULL, ...){
  if (!is.null(digits)) if (!digits %% 1 ==0) stop("digits must be integer")
  levmix <- levelmix(design)
  if (max(levmix$ss) > 15) warning("at least one factor has more than 15 levels\nSomething wrong?")
  hilf <- distDistmix(design, levmix)
  Bd <- Bprime(hilf$BSep, nmax=kmax)
  A <- dualDistmix(Bd, nmax=kmax)
  names(A) <- 0:kmax
  if (!is.null(digits)) A <- round(A, digits)
  if (attrib.out){
     level.info <- rbind(nlevels=levmix$ss, nfactors=levmix$ns)
     colnames(level.info) <- rep("",ncol(level.info))
     attr(A, "B") <- hilf$B
     attr(A, "level.info") <- level.info
  }
  A
}