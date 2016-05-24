# Force loading of dynamic library 

# library.dynam("tdthap")
# Take pedfile input and convert to transmitted and untransmitted haplotypes
# Output is a dataframe

"hap.transmit" <- function(pedfile, markers=1:((ncol(pedfile)-6)/2),
                         multiple.cases=0, use.affected=TRUE) {
  if (!is.loaded("hap_transmit")) stop("hap_transmit C-function not loaded")
  n <- nrow(pedfile)
  m <- length(markers)
  select <- rep(2*markers, rep(2,m)) + rep(5:6, m)
  muse <- as.matrix(pedfile[,select])
  storage.mode(muse) <- "integer"
  res <- .C("hap_transmit", n=as.integer(n), ped=as.integer(pedfile[,1]),
            id=as.integer(pedfile[,2]), father=as.integer(pedfile[,3]),
            mother=as.integer(pedfile[,4]), sex=as.integer(pedfile[,5]),
            aff=as.integer(pedfile[,6]), if.qt=as.integer(FALSE), qt=double(),
            m=as.integer(m), mark=as.integer(muse),
            opt1=as.integer(multiple.cases), opt2=as.integer(use.affected),
            file ="hap-transmit-temporary-file",PACKAGE="tdthap")
  n <- res$n
  if (n==0) stop("No haplotypes could be created")
  nm <- n*m
  res <- .C("hap_read", n=as.integer(n), ped=integer(n), id=integer(n),
            father=integer(n), mother=integer(n),
            as.integer(FALSE), double(), m=as.integer(m),
            f.tr=matrix(integer(nm), nrow=n),
            f.un=matrix(integer(nm), nrow=n),
            m.tr=matrix(integer(nm), nrow=n),
            m.un=matrix(integer(nm), nrow=n),
            file=res$file,PACKAGE="tdthap")
  data.frame(ped=res$ped, id=res$id, father=res$father, mother=res$mother,
             f.tr=res$f.tr, f.un=res$f.un,
             m.tr=res$m.tr, m.un=res$m.un)
}

# Select a subset of markers, find unique haplotypes, and compute transmissions
# informative for the tdt

"tdt.select" <-  function(hap.data, markers=1:((ncol(hap.data)-4)/4),
                        complete=TRUE) {
  n <- nrow(hap.data)
  max.mark <- (ncol(hap.data)-4)/4
  rnames <- paste(hap.data[,1], hap.data[,2], hap.data[,3], hap.data[,4],
                  sep=".")
  for (hap in 1:4) {
    first <- TRUE
    hapcomp <- rep(TRUE, n)
    for (m in markers) {
      col <- 4 + m + (hap-1)*max.mark
      if (first) {
        first <- FALSE
        short.hap <- as.character(hap.data[,col])
      }
      else {
        short.hap <- paste(short.hap, hap.data[,col], sep=".")
      }
      hapcomp <- hapcomp & (hap.data[,col]!=0)
    }
    if (complete) {
      short.hap <- ifelse(hapcomp, short.hap, NA)
    }
    switch(hap,
           f.tr <- short.hap,
           f.un <- short.hap,
           m.tr <- short.hap,
           m.un <- short.hap)
  }
  all <- c(f.tr, f.un, m.tr, m.un)
  all <- all[!is.na(all)]
#  all <- all[all!=NA] Doesn't work in R
  hap.unique <- sort(unique(all))
  i.f.tr <- match(f.tr, hap.unique)
  i.f.un <- match(f.un, hap.unique)
  i.m.tr <- match(m.tr, hap.unique)
  i.m.un <- match(m.un, hap.unique)
  names(i.f.tr) <- paste(rnames,"p", sep=".")
  names(i.f.un) <- paste(rnames,"p", sep=".")
  names(i.m.tr) <- paste(rnames,"m", sep=".")
  names(i.m.un) <- paste(rnames,"m", sep=".")
  f.inf <- !is.na(i.f.tr) & !is.na(i.f.un) & (i.f.tr!=i.f.un)
  m.inf <- !is.na(i.m.tr) & !is.na(i.m.un) & (i.m.tr!=i.m.un)
  nhap <- length(hap.unique)
  tr <- factor(c(i.f.tr[f.inf], i.m.tr[m.inf]),
               levels=1:nhap)#, labels=hap.unique)
  un <- factor(c(i.f.un[f.inf], i.m.un[m.inf]),
               levels=1:nhap)#, labels=hap.unique)
  list(trans=tr, untrans=un)
}

# Function to do Bayes estimate of RR

"tdt.rr" <- function(hap, prior.df=c(0.5,0.5), prob=c(0.05, 0.95)){
  add <- 0
  names(add) <- "Prior"
  Trans <- c(add, table(hap$trans))
  Untrans <- c(add, table(hap$untrans))
  n <- Trans+Untrans
  r <- ifelse(n==0, NA, ifelse(Trans<Untrans, Untrans, Trans))
  p1 <- ifelse(n>0, 1-pf(r/(n-r+1), 2*(n-r+1), 2*r), NA)
  p.value <- ifelse(p1<0.5, 2*p1, 1)
  RR.est <- (prior.df[1]+Trans)/(prior.df[2]+Untrans)
  pi <- qbeta(prob[1], prior.df[1]+Trans, prior.df[2]+Untrans)
  RR.low <- pi/(1-pi)
  pi <- qbeta(prob[2], prior.df[1]+Trans, prior.df[2]+Untrans)
  RR.high <- pi/(1-pi)
  data.frame(Trans, Untrans, p.value, RR.est, RR.low, RR.high)
}

# Interface function to do the quadratic form test

"tdt.quad" <- function(hap, nsim=5000, funct=FALSE, keep=TRUE, seeds=c(0,0,0)) {
  if (!is.loaded("tdt_quad"))
    stop("tdt_quad C-function not loaded")
  if (funct && !is.loaded("tdt_similarity"))
    stop("tdt_similarity C-function not loaded")
  if (keep)
    lres <- nsim + 1
  else
    lres <- 2
  if (seeds[1]==0 && seeds[1]==0 && seeds[2]==0)
    seeds <- runif(3,1,30000)
  haplotypes <- levels(hap$trans)
  nhap <- length(haplotypes)
  ntr <- length(hap$tr)
  c.res <- .C("tdt_quad", as.integer(nhap), as.integer(ntr),
              as.character(haplotypes),
              as.integer((hap$trans)), as.integer((hap$untrans)),
              as.integer(nsim), as.integer(funct), as.logical(keep),
              as.integer(seeds), double(lres), PACKAGE="tdthap")
  obs <- c.res[[10]][1]
  if (is.na(obs)||obs < 0.0)
    stop("No test computed")
  mcs <- c.res[[10]][2:lres]
  if (keep) {
    ep <- sum(mcs>=obs)/nsim
    return(list(n.hap=nhap, n.trans=ntr, test=obs, p.value=ep, sim=mcs))
  }
  else {
    return(list(n.hap=nhap, n.trans=ntr, test=obs, p.value=mcs))
  }
}

"set.similarity" <- function(nloci=1, spacing=rep(1,nloci+1), focus=1, power=1){
  if (!is.loaded("set_tdt_similarity"))
    stop("set_tdt_similarity C-function not loaded")
  if (length(spacing)!=nloci+1 || any(spacing<0))
    stop("Illegal vector of spacings")
  if (power<=0)
    stop("Power must be positive")
  res <- .C("set_tdt_similarity", nloci=as.integer(nloci),
            spacing=as.double(spacing), focus=as.integer(focus),
            power=as.double(power),PACKAGE="tdthap")
  res
}

"get.similarity" <- function(nloci=1) {
  if (!is.loaded("get_tdt_similarity"))
    stop("get_tdt_similarity C-function not loaded")
  res <- .C("get_tdt_similarity", nloci=as.integer(1), 
            spacing=as.double(rep(1,nloci+1)), focus=as.integer(1),
            power=as.double(1),PACKAGE="tdthap")
  res
}
