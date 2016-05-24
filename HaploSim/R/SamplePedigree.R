is.ordered.ped <- function(ped)
  {
    if(ncol(ped)<3)
      stop("pedigree should have at least 3 columns")
    colnames(ped)[1:3] <- c("ID","P0","P1")
    ped$r <- 1:nrow(ped)
    ped$m0 <- match(ped$P0,ped$ID,nomatch = 0)
    ped$m1 <- match(ped$P1,ped$ID,nomatch = 0)
    return(with(ped,all(all(m0 <= r)|all(m1 <= r))))
}
  
SamplePedigree <- function(orig,ped,...)
  {
    if(!is.data.frame(ped))stop("ped should be data.frame")
    if(ncol(ped)<3)stop("ped should have at least tree columns")
    if(!validhaploListObject(orig))stop("orig is not a valid object of class haploList")
    if(length(orig)<2)
      stop("provide at least two haplotypes for the base population")
    hID <- max(sapply(orig,function(x)x@hID))
    ## check if pedigree is ordered
    if(!is.ordered.ped(ped))
      stop("Pedigree is not ordered, first order it.")
    ped$hID1 <- ped$hID0 <- 0
    hList <- haploList(hList = orig)
    Call <- match.Call(Call = match.call(),SampleHaplotype)
    if(is.null(Call$nDec))
      Call$nDec <- orig@nDec
    if(is.null(Call$genDist))
      Call$genDist <- orig@genDist
    for(i in 1:nrow(ped))
      for(p in 1:2)
        {
          pID <- match(ped[i,1+p],ped[,1])
          if(!is.na(pID))
            {
              Call$H0 <- hList[[as.character(ped$hID0[pID])]]
              Call$H1 <- hList[[as.character(ped$hID1[pID])]]
            }
          else
            {
              hh <- sample(1:length(orig),size = 2)
              Call$H0 <- orig[[hh[1]]]
              Call$H1 <- orig[[hh[2]]]
            }
          haplo <- do.call(SampleHaplotype,Call)
          hID <- hID + 1
          haplo@hID <- hID
          if(p == 1)
            ped$hID0[i] <- haplo@hID
          else
            ped$hID1[i] <- haplo@hID
          hList[[as.character(haplo@hID)]] <- haplo
        }
    return(list(ped = ped,hList = hList))
  }
