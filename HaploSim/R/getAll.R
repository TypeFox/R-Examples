## nonpublic function, useful for other functions

getnTrait <- function(hList)
  {
    qtl <- sapply(hList,function(x)length(x@qtl)>0)
    if(!any(qtl))
      return(0)
    nTrait <- sapply(hList[qtl],function(x)length(x@qtl[[1]]))
    nTrait <- unique(nTrait)
    if(length(nTrait)>1)
      stop("Haplotypes with a different number of traits in data")
    else
      return(nTrait)
  }

getAll <- function(hList,what = c("snp","qtl"),removeHomozygotes = FALSE,translatePos = TRUE)
  {
    if(!validhaploListObject(hList))stop("hList is not a valid object of class haploList")
    what <- match.arg(what)
    if(what == "snp")
      {
        loci <- unique(unlist(lapply(hList,function(x)x@snp)))
        loci <- loci[order(loci)]
        all <- matrix(0,nrow = length(hList),ncol = length(loci))
        if(translatePos)
          colnames(all) <- round(loci/(10^hList@nDec),hList@nDec)
        else
          colnames(all) <- loci
        for(i in  1:length(hList))
          all[i,match(hList[[i]]@snp,loci)] <- 1
        if(removeHomozygotes)
          all <- all[,apply(all,2,function(x)length(unique(x))==2)]
      }
    else
      {
        nTrait <- getnTrait(hList)
        if(nTrait==0)
          return(NULL)
        loci <- as.numeric(unique(unlist(lapply(hList,function(x)names(x@qtl)))))
        loci <- loci[order(loci)]
        all <- array(0,dim = c(length(hList),length(loci),nTrait))
        if(translatePos)
          colnames(all) <- round(loci/(10^hList@nDec),hList@nDec)
        else
          colnames(all) <- loci
        for(i in 1:length(hList)){
          if(length(hList[[i]]@qtl)>0)
            all[i,match(names(hList[[i]]@qtl),loci),] <- as.numeric(matrix(unlist(hList[[i]]@qtl),ncol = nTrait,byrow = T))
        }
      }

    return(all)
  }

## getAll <- function(hList,what = c("snp","qtl"),removeHomozygotes = TRUE)
##   {
##     if(!validhaploListObject(hList))stop("hList is not a valid object of class haploList")
##     what <- match.arg(what)
##     if(what == "snp")
##       {
##         loci <- unique(unlist(lapply(hList,function(x)x@snp)))
##         loci <- loci[order(loci)]
##         all <- matrix(0,nrow = length(hList),ncol = length(loci))
##         colnames(all) <- round(loci/(10^hList@nDec),hList@nDec)        
##         for(i in  1:length(hList))
##           all[i,match(hList[[i]]@snp,loci)] <- 1
##         if(removeHomozygotes)
##           all <- all[,apply(all,2,function(x)length(unique(x))==2)]
##       }
##     else
##       {
##         nTrait <- getnTrait(hList)
##         if(nTrait==0)
##           return(NULL)
##         loci <- as.numeric(unique(unlist(lapply(hList,function(x)names(x@qtl)))))
##         loci <- loci[order(loci)]
##         all <- array(0,dim = c(length(hList),length(loci),nTrait))
##         colnames(all) <- round(loci/(10^hList@nDec),hList@nDec)
##         for(i in 1:length(hList)){
##           if(length(hList[[i]]@qtl)>0)
##             all[i,match(names(hList[[i]]@qtl),loci),] <- as.numeric(matrix(unlist(hList[[i]]@qtl),ncol = nTrait,byrow = T))
##         }
##       }

##     return(all)
##   }


    
