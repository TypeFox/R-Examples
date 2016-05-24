dicedist <- function(regmat) 
{
    nart <- ncol(regmat)
    jdist <- rep(0, nart * nart)
    dim(jdist) <- c(nart, nart)
    for (i in 1:(nart - 1)) {
        for (j in (i + 1):nart) {
            uij <- sum(regmat[, i] + regmat[, j] >= 0.99)
            srij <- sum(regmat[, i] + regmat[, j] >= 1.99)
            jdist[j, i] <- jdist[i, j] <- 1 - 2 * srij/(2*srij+uij)
            if (is.na(jdist[i, j])) 
                cat("Warning! NA at i=", i, ", j=", j, "\n")
        }
    }
    jdist
}


lociplots <- function(indclust,locclust,locprab,lcluster,
                      symbols=NULL,brightest.grey=0.8,darkest.grey=0,
                      mdsdim=1:2){
  if (is.null(symbols))
    symbols <- indclust$symbols
  locfreq <- greyvector <- numeric(0)
  clustloci <- rep(FALSE,length(locprab$nonempty.species))
  clustloci[locprab$nonempty.species] <- locclust==lcluster
  locn <- sum(clustloci)
#  print(clustloci)
  locfreqmin <- locfreqmax <- locfreqmean <- numeric(0)
  for (i in 1:locprab$n.regions)
    locfreq[i] <- sum(clustloci[locprab$nonempty.species] &
                      as.logical(locprab$prab[i,]))/locn
#  print(locfreq)
  ic <- indclust$clustering
  if (min(indclust$clustering)==0)
    ic <- ic+1
  for (i in 1:max(ic)){
    ifreqs <- locfreq[ic==i]
#    print(i)
#    print(ifreqs)
    locfreqmin[i] <- min(ifreqs)
    locfreqmax[i] <- max(ifreqs)
    locfreqmean[i] <- mean(ifreqs)
  }
  lmin <- min(locfreqmin)
  lmax <- max(locfreqmax)
#  print(lmin)
#  print(lmax)
  for (i in 1:locprab$n.regions){
    if (lmax-lmin>0)
      greyvector[i] <- darkest.grey+
        (brightest.grey-darkest.grey)*(lmax-locfreq[i])/(lmax-lmin)
    else
      greyvector[i] <- rep(1,length(locfreq))
  }
#  print(greyvector)
  plot(indclust$points[,mdsdim],pch=symbols,col=grey(greyvector),
       xlab=paste("MDS dimension ",mdsdim[1]),
       ylab=paste("MDS dimension ",mdsdim[2]),
       main=paste("Loci cluster ",lcluster, " max.loc.freq.=",
         format(lmax,digits=3)))
  out <- list(locfreq=locfreq,locfreqmin=locfreqmin,locfreqmax=locfreqmax,
              locfreqmean=locfreqmean)
  out
}



# set.seed(5678)
# indprab <- prabinit(file="MartinezOrtega04AFLP.dat", distance = "jaccard")
# indclust <- prabclust(indprab, mdsmethod="kruskal")
# 
# locprab <- prabinit(file="MartinezOrtega04AFLP.dat", rows.are.species = FALSE)
# locclust <- prabclust(locprab, mdsmethod="kruskal") 
# 
# lp1 <- lociplots(indclust,locclust,indprab,locprab,lcluster=1)
# lp2 <- lociplots(indclust,locclust,indprab,locprab,lcluster=2)
# lp3 <- lociplots(indclust,locclust,indprab,locprab,lcluster=3)
# 
# symvec <- 1:10
# lp3 <- lociplots(indclust,locclust,indprab,locprab,lcluster=3,
#                  symbols=symvec[indclust$clustering])
# 
# 
# (etc.) 

# test
# nb <- list()
# data(kykladspecreg)
# for (i in 1:80)
#    nb <- c(nb,list(scan(file="kykladnbcol.dat",
#                    skip=i-1,nlines=1)))
# kycol <- prabinitcols(prabmatrix=kykladspecreg, neighborhood=nb)
# kytest <- prabtestcols(kycol, times=50, pd=0.95)
# summary(kytest)
# kytest2 <- prabtestcols(kycol, times=50, ignore.richness=TRUE)
# summary(kytest2)

# Hier ist der Kram aus Marshal06.R

allelepaircomp <- function(allelepair1,allelepair2,method="sum"){
#  if (nachar %in% allelepair1 | nachar %in% allelepair2) out <- NA
  out <- sum(allelepair1==allelepair2)
  if (!is.na(out) & out==0)
    out <- sum(allelepair1==allelepair2[c(2,1)])
  if (!is.na(out) & method=="geometrical"){
    if (out==1 & ((allelepair1[1]==allelepair1[2])|
          (allelepair2[1]==allelepair2[2])))
      out <- sqrt(0.5)
    else
      out <- out/2
  }    
  out
}

# Expects list of paired vectors
alleledist <- function(allelelist,ni,np,count=FALSE){
#  ni <- nrow(alleleframe)  individuals
#  np <- ncol(alleleframe)
  jdist <- matrix(0,nrow=ni,ncol=ni)
  for (i in 1:(ni-1)){
    if (count) cat("Dissimilarity: processing row ",i,"\n")
    for (j in (i+1):ni){
      agreements <- numeric(0)
      for (k in 1:np)
        agreements[k] <- allelepaircomp(allelelist[[k]][[i]],allelelist[[k]][[j]])
      sagreements <- sum(agreements,na.rm=TRUE)
      snna <- sum(!is.na(agreements))
      if (snna>0)
        jdist[i,j] <- jdist[j,i] <- 1-sagreements/(2*snna)
      else{
        jdist[i,j] <- jdist[j,i] <- 1
        cat("\n i=",i," j=",j,"\n")
        warning("Too many NAs, distance cannot be computed, set to 1.")
      }
    }
  }
  jdist
}

# # Expects list of paired vectors
# chorddist <- function(allelelist,ni,np,count=FALSE){
# #  ni <- nrow(alleleframe)  individuals
# #  np <- ncol(alleleframe)
#   jdist <- matrix(0,nrow=ni,ncol=ni)
#   for (i in 1:(ni-1)){
#     if (count) cat("Dissimilarity: processing row ",i,"\n")
#     for (j in (i+1):ni){
#       agreements <- numeric(0)
#       for (k in 1:np)
#         agreements[k] <- allelepaircomp(allelelist[[k]][[i]],allelelist[[k]][[j]],method="geometrical")
#       sagreements <- sum(sqrt(2*(1-agreements)),na.rm=TRUE)
#       snna <- sum(!is.na(agreements))
#       if (snna>0)
#         jdist[i,j] <- jdist[j,i] <- 2*sagreements/(pi*snna)
#       else{
#         jdist[i,j] <- jdist[j,i] <- 1
#         cat("\n i=",i," j=",j,"\n")
#         warning("Too many NAs, distance cannot be computed, set to 1.")
#       }
#     }
#   }
#   jdist
# }
# 
# # Expects list of paired vectors
# neidist <- function(allelelist,ni,np,count=FALSE){
# #  ni <- nrow(alleleframe)  individuals
# #  np <- ncol(alleleframe)
#   jdist <- matrix(0,nrow=ni,ncol=ni)
#   for (i in 1:(ni-1)){
#     if (count) cat("Dissimilarity: processing row ",i,"\n")
#     for (j in (i+1):ni){
#       agreements <- numeric(0)
#       for (k in 1:np)
#         agreements[k] <- allelepaircomp(allelelist[[k]][[i]],allelelist[[k]][[j]],method="geometrical")
#       sagreements <- sum(agreements,na.rm=TRUE)
#       snna <- sum(!is.na(agreements))
#       if (snna>0)
#         jdist[i,j] <- jdist[j,i] <- 1-sagreements/snna
#       else{
#         jdist[i,j] <- jdist[j,i] <- 1
#         cat("\n i=",i," j=",j,"\n")
#         warning("Too many NAs, distance cannot be computed, set to 1.")
#       }
#     }
#   }
#   jdist
# }
# 

build.ext.nblist <- function(neighbors,n.individuals=length(neighbors))
{
    ext.nblist <- list()
    for (i in seq(1,2*n.individuals-1,2)){
      ext.nblist[[i]] <- ext.nblist[[i+1]] <-
        c(2*neighbors[[floor((i+1)/2)]],2*neighbors[[floor((i+1)/2)]]-1)
      ext.nblist[[i]] <- c(ext.nblist[[i]],i+1)
      ext.nblist[[i+1]] <- c(ext.nblist[[i+1]],i)
    }
    ext.nblist
}
      
build.charmatrix <- function(allelelist,n.individuals,n.variables)
{
  cm <- matrix(nrow=2*n.individuals,ncol=n.variables)
  for (i in 1:n.individuals)
    for (k in 1:n.variables)
      cm[2*i-2+(1:2),k] <- allelelist[[k]][[i]]
  cm
}

unbuild.charmatrix <- function(charmatrix,n.individuals,n.variables)
{
  allelelist <- list()
  for (i in 1:n.variables){
    allelelist[[i]] <- list()
    for (k in 1:n.individuals){
      allelelist[[i]][[k]] <- charmatrix[2*k-2+(1:2),i]
      if (is.na(charmatrix[2*k-1,i])) allelelist[[i]][[k]] <- NA
    }
  }
  allelelist
}

nastats <- function(amatrix, nastr="--"){
  nano <- function(avector, nastr)
    sum(avector==nastr)
  
  narow <- apply(amatrix, 1, nano, nastr=nastr)
  nacol <- apply(amatrix, 2, nano, nastr=nastr)
  out <- list(narow=narow,nacol=nacol)
  out
}

# alleleconvert IS DOCUMENTED
# if file ("path/filename") is specified, data is read from file
# otherwise R-object strmatrix is read
# format.in can be "genepop" or "structure"
# format.out can be "prabclus" or "structurama"
# alength: length of allele entry in original genepop file
# orig.nachar: how are NAs specified in original data
# new.nachar: how are NAs specified in output
# firstcolname: ist first column of input individual name? (If so,
# will be added to output file)
# aletters: symbols for alleles in "prabclus" format
# outfile: filename to write output; if NULL, no file is written
alleleconvert <- function(file=NULL,strmatrix=NULL, format.in="genepop",
                          format.out="prabclus",
                          alength=3,orig.nachar="000",new.nachar="-",
                          rows.are.individuals=TRUE, firstcolname=FALSE,
                          aletters=intToUtf8(c(65:90,97:122),multiple=TRUE),
                          outfile=NULL){
  if (is.null(file))
    m1 <- strmatrix
  else
    m1 <- read.table(file,colClasses="character")
  if (firstcolname){
    if (format.in=="structure")
      colnamesv <- m1[seq(1,nrow(m1)-1,by=2),1]
    if (format.in=="genepop")  
      colnamesv <- m1[,1]
    m1 <- m1[,2:ncol(m1)]
  }
  if (!rows.are.individuals) 
        m1 <- t(m1)
  if (format.in=="genepop"){
    n.individuals <- nrow(m1)
    n.variables <- ncol(m1)
  }
  if (format.in=="structure"){
    n.individuals <- nrow(m1)/2
    n.variables <- ncol(m1)
  }
#  str(m1)
#  print(aletters)
  first.allele <- second.allele <- outmatrix <- matrix("",ncol=n.variables,
                                          nrow=n.individuals)
  double.out <- matrix("",ncol=n.variables, nrow=2*n.individuals)
  if (format.in=="genepop"){
    for (i in 1:n.individuals)
      for (j in 1:n.variables){
        first.allele[i,j] <- substring(m1[i,j],1,alength)
        second.allele[i,j] <- substring(m1[i,j],alength+1,2*alength)
      }
  }
  if (format.in=="structure"){
    for (i in 1:n.individuals)
      for (j in 1:n.variables){
        first.allele[i,j] <- m1[2*i-1,j]
        second.allele[i,j] <- m1[2*i,j]
      }
  }     
  double.allele <- rbind(first.allele,second.allele)
#  print(double.allele[c(1:10,300:319),])
  maxk <- 0
  for (j in 1:n.variables){
    jlevels <- levels(as.factor(double.allele[,j]))
    jlevels <- jlevels[jlevels!=orig.nachar]
    if (length(jlevels)>maxk) maxk <- length(jlevels)
#    print(jlevels)
    double.out[double.allele[,j]==orig.nachar,j] <- new.nachar
    if (format.out=="prabclus"){
      for (k in 1:length(jlevels))
        double.out[double.allele[,j]==jlevels[k],j] <- aletters[k]
      outmatrix[,j] <- mapply(paste,double.out[1:n.individuals,j],
                            double.out[(n.individuals+1):(2*n.individuals),j],
                            sep="")
      out <- structure(outmatrix,alevels=aletters[1:maxk])
    }
    if (format.out=="structurama"){
      for (k in 1:length(jlevels))
        double.out[double.allele[,j]==jlevels[k],j] <- jlevels[k]
      outmatrix[,j] <- mapply(paste,"(",double.out[1:n.individuals,j],",",
                            double.out[(n.individuals+1):(2*n.individuals),j],
                              ")",
                            sep="")
      out <- structure(outmatrix,alevels=jlevels)
    }
      
#    print(outmatrix[,j])
  }
  if(!is.null(outfile)){
    if (firstcolname)
      outmatrix <- cbind(colnamesv,outmatrix)
    write.table(outmatrix,outfile,quote=FALSE,col.names=FALSE)
  }     
  out
}

# alleleobject$amatrix rows are individuals
allele2zeroone <- function(alleleobject){
  lnumbers <- cumlnumbers <- numeric(0)
  cumlnumbers[1] <- 0
  for (i in 1:alleleobject$n.variables){
    lnumbers[i] <- sum(alleleobject$leveldist[,i]>0)
    cumlnumbers[i+1] <- cumlnumbers[i]+lnumbers[i]
  }
  nzeroone <- cumlnumbers[alleleobject$n.variables+1]
  imatrix <- matrix(0,nrow=alleleobject$n.individuals,ncol=nzeroone)
  for (i in 1:alleleobject$n.variables)
    for (j in 1:lnumbers[i])
      for (k in 1:alleleobject$n.individuals)
        imatrix[k,cumlnumbers[i]+j] <-
          as.integer(alleleobject$alevels[j] %in%
                     alleleobject$charmatrix[((k-1)*2+1):(k*2),i])
  imatrix
}

# See allelepop.nb for namode
# This assumes diploid organisms
# If nachar is only one of a pair, both will be set to nachar
# Note that the first index for allelelist is variable, the second one is
# individual. This is different from amatrix.
# Note! Taking allelematrix from alleleconvert is needed to guarantee a
# level coding of alleles starting from the same levels (usually "A")
# for all loci, which prevents allele2zeroone from producing empty columns.
# Example:
# alleleconvert(strmatrix=alleleobject$amatrix,alength=1,orig.nachar="-")
# distance can be "alleledist", "chord", "ney" and "none".
alleleinit <- function (file = NULL, allelematrix=NULL,
                        rows.are.individuals = TRUE, 
    neighborhood = "none", distance = "alleledist", namode="variables",
                       nachar="-", distcount=FALSE) 
{
    ieq <- function(level,fvector)
      sum(fvector==level,na.rm=TRUE)

    if(is.null(allelematrix))
      m1 <- read.table(file,stringsAsFactors=FALSE)
    else
      m1 <- as.data.frame(allelematrix,stringsAsFactors=FALSE)
    if (!rows.are.individuals) 
        m1 <- t(m1)    
#    m1 <- as.matrix(m1[, regperspec > 0])
#    regperspec <- regperspec[regperspec > 0]
    n.individuals <- nrow(m1)
    n.variables <- ncol(m1)
    nado <- !(namode=="none")
    napair <- paste(nachar,nachar,sep="")
    naprob=NULL
    
    allelelist <- list()
    for (i in 1:n.variables){
      allelelist[[i]] <- list()
      for (j in 1:n.individuals){
        allelelist[[i]][[j]] <- c(substring(m1[j,i],1,1),substring(m1[j,i],2,2))
        if (nado & nachar %in% allelelist[[i]][[j]]){
          allelelist[[i]][[j]] <- NA
          m1[j,i] <- napair
        }
      }
    }
    
    if (nado){
      nam1 <- nastats(m1,nastr=napair)
      nasum <- sum(nam1$narow)
      nallele <- n.individuals*n.variables
      if (namode=="single")
        naprob <- nasum/nallele
      if (namode=="individuals")
        naprob <- nam1$narow/n.variables
      if (namode=="variables")
        naprob <- nam1$nacol/n.individuals
    }
    
    prab <- t(as.matrix(as.data.frame(lapply(lapply(m1,as.factor),as.integer))))
    regperspec <- apply(prab, 2, sum)
    specperreg <- apply(prab, 1, sum)

    if (is.null(attr(allelematrix,"alevels"))){
      alevels <- levels(as.factor(unlist(allelelist)))
    }
    else
      alevels <- attr(allelematrix,"alevels")
    n.levels <- length(alevels)
    leveldist <- matrix(0,ncol=n.variables,nrow=n.levels)
    for (i in 1:n.variables){
      leveldist[,i] <- sapply(alevels,ieq,unlist(allelelist[[i]]))
      if (sum(leveldist[,i])==0) warning("Only NAs in variable ",i)
    }
                
    nb <- list()
    if (is.list(neighborhood)) 
        nb <- neighborhood
    else {
        if (identical(neighborhood,"none")) 
            for (i in 1:n.individuals) nb[[i]] <- numeric(0)
        else for (i in 1:n.individuals)
          nb <- c(nb, list(scan(file = neighborhood, 
            skip = i - 1, nlines = 1, quiet = TRUE)))
    }
    nbtest(nb, n.individuals)
    
    distmat <- switch(distance, alleledist = alleledist(allelelist,
                                  n.individuals,n.variables, count=distcount),
#                      chord = chorddist(allelelist,
#                                  n.individuals,n.variables, count=distcount),
#                      nei = neidist(allelelist,
#                                  n.individuals,n.variables, count=distcount),
                      none = NULL)

    ext.nblist <- build.ext.nblist(nb, n.individuals)
    charmatrix <- build.charmatrix(allelelist, n.individuals, n.variables)
    
    out <- list(distmat = distmat, amatrix = m1,
#                allelelist=allelelist,
                charmatrix=charmatrix,
                nb = nb, ext.nblist=ext.nblist,
                n.variables = n.variables, n.individuals = n.individuals,
                n.levels=n.levels, 
                n.species=n.individuals,
                alevels=alevels, leveldist=leveldist,
                prab=prab,
                regperspec=regperspec, specperreg=specperreg,
        distance = distance, namode=namode, naprob=naprob, nasum=nasum,
                nachar=nachar,
                spatial = (!identical(neighborhood, 
            "none")))
    class(out) <- "alleleobject"
    out
  }


print.alleleobject <- function (x, ...) 
{
    cat("Diploid allele/locus-data object\n")
    cat("Number of individuals: ",x$n.individuals,"\n")
    cat("Number of loci: ",x$n.variables,"\n")
    cat("Alleles (all loci): ",x$alevels,"\n")
    if (x$distance!="none")
      cat("Object contains between-individuals dissimilarity matrix of type ",
          x$distance,".\n")
    if (x$spatial)
      cat("Object contains neighborhood list nb between individuals.\n")
    cat("Object contains ",x$nasum," missing loci data.\n")
    cat("Further potentially informative components:\n")
    cat(" amatrix (individuals*loci matrix),\n")
    cat(" leveldist (distribution of alleles per locus), \n")
    cat(" naprob (probabilities for missing values, governed by namode=",
        x$namode,"), \n")
    cat("components charmatrix, ext.nblist, prab, regperspec, specperreg don't\n")
    cat("contain additional information and are only needed for efficient data handling.\n")
}



stressvals <-  function(x,mdsdim=1:12,trace=FALSE){
  oregions <- order(x$specperreg)
  prabo1 <- x$prab[oregions, ]
  ospecies <- do.call("order", as.data.frame(t(prabo1)))
  dm <- x$distmat[ospecies, ospecies]
  mindm <- min(dm[dm > 0])/10
  for (i in 1:(x$n.species - 1))
    for (j in (i + 1):x$n.species)
    if (dm[i, j] < mindm) 
            dm[i, j] <- dm[j, i] <- mindm
  mdsout <- list()
  for (mdim in mdsdim)
    mdsout[[mdim]] <- isoMDS(dm, k = mdim,trace=FALSE)
  MDSstress <- numeric(0)
  for (i in mdsdim)
    MDSstress[i] <- mdsout[[i]]$stress
  out <- list(MDSstress=MDSstress,mdsout=mdsout)
  out
}


# file.format could be "degminsec", "decimal4", "decimal2"
coord2dist <- function(file=NULL, coordmatrix=NULL, cut=NULL,
                       file.format="degminsec",
                       output.dist=FALSE, radius=6378.137,
                       fp=1/298.257223563, neighbors=FALSE){
  if (is.null(coordmatrix))
    m1 <- as.matrix(read.table(file))
  else
    m1 <- as.matrix(coordmatrix)
  ap <- radius
  n <- nrow(m1)
  p <- ncol(m1)
#  print(p)
  m2 <- matrix(0,ncol=2,nrow=n)
  distmatrix <- matrix(0,ncol=n,nrow=n)
  if (file.format=="decimal2"){
    if (!(p==2)) stop("decimal2 input needs two columns.")
    else
      m2 <- m1
  }
  if (file.format=="decimal4"){
    if (!(p==4)) stop("decimal4 input needs four columns.")
    else{
      m2[,1] <- sign(m1[,1])*(abs(m1[,1])+m1[,2]/100)
      m2[,2] <- sign(m1[,3])*(abs(m1[,3])+m1[,4]/100)
    }
  }
  if (file.format=="degminsec"){
    if (!(p==6)) stop("degminsec input needs six columns.")
    else{
      m2[,1] <- sign(m1[,1])*(abs(m1[,1]) + (m1[,2] + m1[,3] / 60) / 60)
      m2[,2] <- sign(m1[,4])*(abs(m1[,4]) + (m1[,5] + m1[,6] / 60) / 60)
    }
  }
  m2 <- m2*pi/180
  for (i in 1:(n-1))
    for(j in (i+1):n){
      Fp <- (m2[i,1]+m2[j,1])/2
      Gp <- (m2[i,1]-m2[j,1])/2
      lp <- (m2[i,2]-m2[j,2])/2
      Sp <- sin(Gp)^2*cos(lp)^2 +cos(Fp)^2*sin(lp)^2
      Cp <- cos(Gp)^2*cos(lp)^2 +sin(Fp)^2*sin(lp)^2
      wp <- atan(sqrt(Sp/Cp))
      Rp <- sqrt(Sp*Cp)/wp
      Dp <- 2*wp*ap
      H1 <- (3*Rp-1)/(2*Cp)
      H2 <- (3*Rp+1)/(2*Sp)
      if (Gp==0 & lp==0)
        distmatrix[i,j] <- distmatrix[j,i] <- 0
      else
        distmatrix[i,j] <- distmatrix[j,i] <-
          Dp*(1 + fp*H1*sin(Fp)^2*cos(Gp)^2 - fp*H2*cos(Fp)^2*sin(Gp)^2)
    }
  if (output.dist) distmatrix <- as.dist(distmatrix)
  if (neighbors){
    nblist <- geo2neighbor(distmatrix,cut)
    out <- list(distmatrix=distmatrix,nblist=nblist)
  }
  else
    out <- distmatrix
  out
}


  
      
