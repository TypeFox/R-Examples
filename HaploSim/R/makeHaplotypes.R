hPed2Ped <- function(hPed)
  {
    hPed[!hPed[,2]%in%hPed[,1],2] <- NA
    hPed[!hPed[,3]%in%hPed[,1],3] <- NA
    hPed[apply(hPed[,-1],1,function(x)any(is.na(x))),-1] <- rep(NA,2)
                                        # only one missing haplotype is not possible
    hPed[,-1] <- t(apply(hPed[,-1],1,function(x)x[order(x)]))
                                        # ordering parental haplotypes avoids finding inverted identical combinations
    hPairs <- apply(hPed[,2:3],1,paste,collapse = "/")
    uPairs <- !duplicated(hPairs)&hPairs!="NA/NA"
    nPairs <- sum(uPairs)
    uHaplos <- sapply(hPed[,1],function(x)!(x%in%c(hPed[,2],hPed[,3])))
    nnPair <- sum(uHaplos)
    nInd <-   nPairs + nnPair
    
    ped <- data.frame(matrix(NA,ncol = 5,nrow = nInd))
    names(ped) <- c("ID","p0","p1","H0","H1")

    ped[1:nPairs,4] <- hPed[uPairs,2]
    ped[1:nPairs,5] <- hPed[uPairs,3]

    ped[nPairs + seq(1,nnPair),4] <- hPed[uHaplos,1]

    pPed <- apply(ped[,4:5],1,paste,collapse = "/")
    mm <- match(ped[,4],hPed[,1])
    p0 <- match(hPairs[mm],pPed)
    mm <- match(ped[,5],hPed[,1])
    p1 <- match(hPairs[mm],pPed) 
    ped[,2] <- p0
    ped[,3] <- p1
    ped[,1] <- 1:nInd
    return(ped)
  }

buildhPedigree <- function(hPedigree = NULL,hList)
  {
    if(!validhaploListObject(hList))stop("hList is not a valid object of class haploList")
    hped <- data.frame(t(sapply(hList,function(x)
                                c(x@hID,x@phID0,x@phID1)))) 
    names(hped) <- c("hID","phID0","phID1")
    hPedigree <- rbind(hPedigree,hped)
    return(hPedigree)
  }

haploList <- function(list = NULL,hList = NULL,nDec,genDist,nChrom = 1)
  {
    if(!is.null(hList))
      {
        if(!validhaploListObject(hList))
          stop("hList shouhld be and object of class haploList")
        hList <- new("haploList",nDec = hList@nDec,genDist = hList@genDist,nChrom = hList@nChrom)
      }
    else
      hList <- new("haploList",nDec = as(nDec,'integer'),genDist = genDist,nChrom = as(nChrom,"integer"))
    if(!is.null(list))
      {
        if(any(sapply(list,function(x)class(x)!="haplotype")))
          stop("List can only contain objects of class haplotype")
        hList[1:length(list)] <- list
      }
    hList
  }

## small function to match calls within a function to an underlying function
match.Call <- function(Call = match.call(),fn)
  {
    Call <- as.list(Call)
    mm <- pmatch(names(formals(fn)),names(Call))
    Call[mm[!is.na(mm)]]
  }
    
SampleHaplotypes <- function(orig = NULL,nHaplotypes = 10,nMeioses = 2,gg = NULL,...)
  {
    hID <- 0
    call <- match.call()
    Call <- match.Call(Call = call,haploList)
    Call$hList <- orig
    hList <- do.call(haploList,Call)
    if(is.null(orig))
      {
        Call <- match.Call(Call = call,SampleBaseHaplotype)
        for(i in 1:nHaplotypes)
          {
            hList[[length(hList) + 1]] <- do.call(SampleBaseHaplotype,Call)
            hID <- hID + 1
            hList[[length(hList)]]@hID <- hID
          }
      }
    else 
      {
        nTraits <- getnTrait(hList) ## if the data is wrong, stops
        hID <- max(sapply(orig,function(x)x@hID))
        if(is.null(gg))
          gg <- sample(1:length(orig),size = length(orig))
        ii <- seq(1,length(gg),by = 2)
        Call <-   match.Call(Call = call,SampleHaplotype)
        Call$nDec <- orig@nDec
        Call$genDist <- orig@genDist
        Call$nChrom <- orig@nChrom
        Call$checkValidity <- FALSE ## otherwise checks everytime
        for(i in ii)
          {
            Call$H0 <- orig[[gg[i]]]
            Call$H1 <- orig[[gg[i + 1]]]
            for(n in 1:nMeioses)
              {
                hList[[length(hList) + 1]] <- do.call(SampleHaplotype,Call)
                hID <- hID + 1
                hList[[length(hList)]]@hID <- hID
              }
          }
      }
    return(hList)
  }

calcRecPos <- function(pos,rec) ## nonpublic function, only for use in  SampleHaplotype
{
  rr <- table(unlist(sapply(rec,function(x)
                            match(TRUE,x<pos))))
  return(as.integer(names(rr))[rr%%2 != 0])
  ##  dubbel recombinant leidt tot geen recombinant
}

getOrig <- function(orig,rec,nLoc) ## nonpublic function, only for use in SampleHaplotype
  {
    orig <- rep(orig,nLoc)
    nRec <- length(rec)
    end <- numeric()
    begin <- rec[seq(1,nRec,by = 2)]
    if(nRec>1)end <- rec[seq(2,nRec,by = 2)] - 1
    if(length(end)==length(begin)-1) end <- c(end,nLoc)
    for(.r in 1:length(begin))
      {
        oo <- orig[begin[.r]]
        orig[begin[.r]:end[.r]] <- 3 - oo
      }
    return(matrix(c(orig,1:nLoc),ncol = 2))
  }

SampleBaseHaplotype <- function(genDist,nDec,nLoc,pSnp = seq(0,1,length.out = nLoc))
  {
    hh <- new('haplotype')
    nSNPS <- genDist * 10^nDec
    if(is.null(pSnp)) ## optie om met pSnp Snp posities te geven
      pSnp <- runif(nLoc,min = 0,max = nSNPS) ## aantal loci en afronding bepaalt werkelijke aantal loci
    else if(max(pSnp)>genDist)stop("SNP positions are only allowed between 0 and genome size")
    else
      {
        pSnp <- pSnp[sample(1:length(pSnp),size = 0.5*length(pSnp),replace = FALSE)] ## makes frequency = 0.5 of each marker locus
        pSnp <- 10^nDec*pSnp
      }
    pSnp <- unique(as.integer(round(pSnp)))
    pSnp <- pSnp[order(pSnp)]
    hh@snp <- pSnp
    return(hh)
  }

## nonpublic function, useful for function SampleHaplotype
validhaplotypePair <- function(H0 = NULL,H1 = NULL,nSNPS){
  if(any(c(max(H0@snp),max(H1@snp))>nSNPS))
    return(FALSE)
  nTrait <- rep(NA,2)
  if(length(H0@qtl)>0)nTrait[1] <- length(H0@qtl[[1]])
  if(length(H1@qtl)>0)nTrait[2] <- length(H1@qtl[[1]])
  nTrait <- nTrait[!is.na(nTrait)]
  if(length(unique(nTrait))==2)
    return(FALSE)
  return(TRUE)
}

SampleHaplotype <- function(H0 = NULL,H1 = NULL,genDist,nDec,nChrom = 1,prMut = 1E-5,QTL = F,checkValidity = TRUE)
  {
    nSNPS <- genDist * 10^nDec
    {
      if(is.null(H0)&is.null(H1))sampleNew <- TRUE
      else if(is.null(H0)|is.null(H1))
        stop("Impossible to sample a meiosis whith only one haplotype object")
      else sampleNew <- FALSE
    }
    if(checkValidity)
      if(!validhaplotypePair(H0,H1,nSNPS)){
        cat("H0 and H1 are unvalid haplotype pair \n")
        cat("in combination with genDist and nDec.\n")
        stop()
      }
    hh <- new('haplotype')
    nRec <- rpois(1,genDist) ## on the Morgan scale
    ## haldane mapping function:
    ## crossover value is function of genetic distance. Probability of a crossover in a segment of length genDist is Poisson(genDist)
   
    rec <- runif(n = nRec,min = 0,max = nSNPS)
    if(nChrom>1){
      chrBreaks <- seq(1,nSNPS,length.out = nChrom + 1)[-c(1,nChrom + 1)]
      rr <- rbinom(1,size = length(chrBreaks),prob = 0.5)
      rec <- c(rec,chrBreaks[sample(1:length(chrBreaks),size = rr)])
    }
    rec <- rec[order(rec)]
    nRecSNP <- nRecQTL <- 0
    orig <- sample(1:2,1)
    ## snp part
    if(nRec>0)
      {
        snpPos <- unique(c(H0@snp,H1@snp))
        snpPos <- snpPos[order(snpPos)]
        rSNP <- calcRecPos(pos = snpPos,rec = rec) #removes double recombinant within segments
        nRecSNP <- length(rSNP)
      }
    if(nRecSNP>0)
      {
        oo <- getOrig(orig,rSNP,length(snpPos))
        snp <- matrix(0,ncol = length(snpPos),nrow = 2)
        snp[1,match(H0@snp,snpPos)] <- 1
        snp[2,match(H1@snp,snpPos)] <- 1
        snp <- snpPos[snp[oo] == 1]
      }
    else ## 
      {
        if(orig == 1)
          snp <- H0@snp
        else
          snp <- H1@snp
      }
    nMut <- rpois(1,prMut*nSNPS)
    if(nMut>0)
      {
        msnp <- unique(as.integer(round(runif(nMut,min = 0,max = nSNPS))))
        bsnp <- msnp[msnp%in%snp] ## snp that were present
        snp <- c(snp,msnp)
        snp <- snp[!snp%in%bsnp] ## mutations on existing snp that mutate back to nonnsp
        snp <- snp[!duplicated(snp)] ## mutation at snp position removes allele (only 1 alleles stored)
      }
    hh@snp <- snp[order(snp)]
    ## qtl part, if required
    if(QTL)
      {
        if(any(c(length(H0@qtl),length(H1@qtl))>0)){
          if(nRec>0)
            {
              qtlPos0 <- as.integer(names(H0@qtl))
              qtlPos1 <- as.integer(names(H1@qtl))
              qtlPos <- unique(c(qtlPos0,qtlPos1))
              qtlPos <- qtlPos[order(qtlPos)]
              rQTL <- calcRecPos(pos = qtlPos,rec = rec) ## also accounts for different chromosomes
              nRecQTL <- length(rQTL)
            }
          if(nRecQTL>0)
            {
              oo <- getOrig(orig,rQTL,length(qtlPos))
              oo[,2] <- qtlPos
              o0 <- oo[oo[,1] == 1,2]
              o1 <- oo[oo[,1] == 2,2]
              qtl <- H0@qtl[match(o0,qtlPos0)]
              qtl <- c(qtl,H1@qtl[match(o1,qtlPos1)])
              qtl <- qtl[!sapply(qtl,is.null)]
              ns <- as.numeric(names(qtl))
              qtl <- qtl[order(ns)]
            }
          else
            {
              if(orig == 1)
                qtl <- H0@qtl
              else
                qtl <- H1@qtl
            }
          hh@qtl <- qtl
        }
      }
    hh@phID0 <- H0@hID
    hh@phID1 <- H1@hID
    return(hh)
  }

RemoveHomozygotes <- function(hList)
  {
    if(!validhaploListObject(hList))stop("hList is not a valid object of class haploList")
    nHaplos <- length(hList)
    lTable <- table(unlist(lapply(hList,function(x)x@snp)))
    hSNPs <- names(lTable)[lTable == nHaplos]
    for(l in 1:length(hList))
      {
        hh <- hList[[l]]
        hh@snp <- hh@snp[!hh@snp%in%hSNPs]
        hList[[l]] <- hh
      }
    return(hList)
  }

AssignSingleQTL <- function (hList, nqtl = NA, frqtl = NA, sigma2qtl = NULL, shQTL = 1, 
                       scQTL = 1, MAF = 0.1, rmCausSNP = TRUE){
  if (!validhaploListObject(hList)) 
    stop("hList is not a valid object of class haploList")
  lTable <- table(unlist(lapply(hList, function(x) x@snp)))/length(hList)
  lTable <- lTable[lTable > MAF & lTable < (1 - MAF)]
  if (is.na(nqtl)) 
    nqtl = round(length(lTable) * frqtl)
  pqtl <- round(seq(1, length(lTable), length.out = nqtl))
  a <- numeric(nqtl)
  if (is.null(sigma2qtl)) 
    sigma2qtl <- rgamma(n = nqtl, shape = shQTL, scale = scQTL)
  else sigma2qtl <- rep(sigma2qtl, length.out = nqtl)
  for (e in 1:nqtl) a[e] <- sqrt(sigma2qtl[e]/(2 * lTable[pqtl[e]] * 
                                               (1 - lTable[pqtl[e]])))
  pqtl <- as.integer(names(lTable[pqtl]))
  for (l in 1:length(hList)) {
    x <- hList[[l]]
    pos <- pqtl[pqtl %in% x@snp]
    qtl <- as.list(a[pqtl %in% x@snp])
    names(qtl) <- as.character(pos)
    x@qtl <- qtl
    if (rmCausSNP) 
      x@snp <- x@snp[!x@snp %in% pqtl]
        hList[[l]] <- x
  }
  return(hList)
}

AssignQTL <- function(hList,nqtl = NA,frqtl = NA,sigma2qtl = NULL,shQTL = 1,scQTL = 1,
                      nTraits = 1,overlap = 0,MAF = 0.1,rmCausSNP = TRUE)
  {
    if(!validhaploListObject(hList))stop("hList is not a valid object of class haploList")
    lTable <- table(unlist(lapply(hList,function(x)x@snp)))/length(hList)
    lTable <- lTable[lTable>MAF&lTable<(1-MAF)]
    if(is.na(nqtl))
      nqtl = round(length(lTable) * frqtl)
    if(is.null(sigma2qtl))
      sigma2qtl <- lapply(1:nTraits,function(x){
        rgamma(n = nqtl,shape = shQTL,scale = scQTL)
      })
    if(is.numeric(sigma2qtl))
      sigma2qtl <- lapply(1:nTraits,function(x)
                          rep(sigma2qtl,length.out = nqtl))
    else if(is.list(sigma2qtl)&length(sigma2qtl)== nTraits)
      sigma2qtl <- lapply(1:nTraits,function(x)
                          rep(sigma2qtl[[x]],length.out = nqtl))
    else{
      cat("sigma2qtl was not specified correctly.\n")
      cat("If specified as a list, length should be equal to number of traits.\n")
      stop()
    }
    if(overlap<0|overlap>1)
      stop("overlap should be between 0 and 1.\n")
    pqtl <- as.integer(seq(1,length(lTable),length.out = 2 + nTraits*nqtl*(1-overlap))) ## qtl positions along the genome
    pqtl <- pqtl[-c(1,length(pqtl))] ## removes the two outer qtl
    qtl <- list()
    a <- numeric(nqtl) ## equal number of QTL for each trait
    for(t in 1:nTraits)
      {
        if(t>1){
          pqtl <- pqtl[!pqtl%in%pos]
          pos <- c(pqtl[seq(1,length(pqtl),length.out = nqtl*(1-overlap))],pos[seq(1,length(pos),length.out = overlap)])
          pos <- pos[order(pos)]
        }
        else
          pos <- pqtl[seq(1,length(pqtl),length.out = nqtl)]
        for(e in 1:nqtl)
          a[e] <- sqrt(sigma2qtl[[t]][e]) / (2*lTable[pos[e]]*(1-lTable[pos[e]]))
        qtl[[t]] <- list(pqtl = as.numeric(names(lTable[pos])),a = a)
      }
    pqtl <- unlist(lapply(qtl,function(x)x$pqtl))
    a <- unlist(lapply(qtl,function(x)x$a))
    trait <- unlist(lapply(1:nTraits,function(x)rep(x,nqtl)))
    qtlT <- tapply(a,list(pqtl,trait),sum)
    qtl <- lapply(1:nrow(qtlT),function(x)qtlT[x,])
    names(qtl) <- rownames(qtlT)
    qtls <- as.numeric(names(qtl))
    for(l in 1:length(hList))
      {
        x <- hList[[l]]
        x@qtl <- qtl[qtls%in%x@snp]
        if (rmCausSNP) 
          x@snp <- x@snp[!x@snp %in% qtls]
        hList[[l]] <- x
      }
    return(hList)
  }

ListQTL <- function(hList,nqtl = NA,frqtl = NA,sigma2qtl = NULL,shQTL = 1,scQTL = 1,
                      nTraits = 1,overlap = 0,MAF = 0.1)
  ## function to get a data.frame of the QTL positions and their effects
  ## does not put them in the object of class haploList
  ## to avoid redundant data
  {
    if(!validhaploListObject(hList))stop("hList is not a valid object of class haploList")
    lTable <- table(unlist(lapply(hList,function(x)x@snp)))/length(hList)
    lTable <- lTable[lTable>MAF&lTable<(1-MAF)]
    if(is.na(nqtl))
      nqtl = round(length(lTable) * frqtl)
    if(is.null(sigma2qtl))
      sigma2qtl <- lapply(1:nTraits,function(x){
        rgamma(n = nqtl,shape = shQTL,scale = scQTL)
      })
    if(is.numeric(sigma2qtl))
      sigma2qtl <- lapply(1:nTraits,function(x)
                          rep(sigma2qtl,length.out = nqtl))
    else if(is.list(sigma2qtl)&length(sigma2qtl)== nTraits)
      sigma2qtl <- lapply(1:nTraits,function(x)
                          rep(sigma2qtl[[x]],length.out = nqtl))
    else{
      cat("sigma2qtl was not specified correctly.\n")
      cat("If specified as a list, length should be equal to number of traits.\n")
      stop()
    }
    if(overlap<0|overlap>1)
      stop("overlap should be between 0 and 1.\n")
    pqtl <- as.integer(seq(1,length(lTable),length.out = 2 + nTraits*nqtl*(1-overlap))) ## qtl positions along the genome
    pqtl <- pqtl[-c(1,length(pqtl))] ## removes the two outer qtl
    qtl <- list()
    a <- numeric(nqtl) ## equal number of QTL for each trait
    for(t in 1:nTraits)
      {
        if(t>1){
          pqtl <- pqtl[!pqtl%in%pos]
          pos <- c(pqtl[seq(1,length(pqtl),length.out = nqtl*(1-overlap))],pos[seq(1,length(pos),length.out = overlap)])
          pos <- pos[order(pos)]
        }
        else
          pos <- pqtl[seq(1,length(pqtl),length.out = nqtl)]
        for(e in 1:nqtl)
          a[e] <- sqrt(sigma2qtl[[t]][e]) / (2*lTable[pos[e]]*(1-lTable[pos[e]]))
        qtl[[t]] <- list(pqtl = as.numeric(names(lTable[pos])),a = a)
      }
    pqtl <- unlist(lapply(qtl,function(x)x$pqtl))
    a <- unlist(lapply(qtl,function(x)x$a))
    trait <- unlist(lapply(1:nTraits,function(x)rep(x,nqtl)))
    qtlT <- tapply(a,list(pqtl,trait),sum)
    qtl <- lapply(1:nrow(qtlT),function(x)qtlT[x,])
    names(qtl) <- rownames(qtlT)
    return(qtl)
  }

        
