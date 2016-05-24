#Function: ghap.profile
#License: GPLv3 or later
#Modification date: 2 Feb 2016
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute individual profiles based on haplotype allele scores

ghap.profile <- function(
  score,
  haplo,
  only.active.samples=TRUE,
  ncores=1
){
  
  
  #Check if haplo is a GHap.haplo object
  if(class(haplo) != "GHap.haplo"){
    stop("Argument haplo must be a GHap.haplo object.")
  }
  
  #Check if inactive alleles and samples should be reactived
  if(only.active.samples == FALSE){
    haplo$id.in <- rep(TRUE,times=haplo$nsamples)
    haplo$nsamples.in<-length(which(haplo$id.in))
  }
  
  #Prepare score vectors
  score$HAPLOtmp <- paste(score$BLOCK,score$CHR,score$BP1,score$BP2,score$ALLELE,sep="_")
  haps <- NULL
  haps$HAPLOtmp <- paste(haplo$block,haplo$chr,haplo$bp1,haplo$bp2,haplo$allele,sep="_")
  haps$IDX <- 1:length(haps$HAPLOtmp)
  haps <- data.frame(haps,stringsAsFactors = FALSE)
  haps <- merge(x = haps, y = score, by.x = "HAPLOtmp", by.y = "HAPLOtmp", all.x=TRUE,sort=FALSE)
  haps <- haps[is.na(haps$SCORE) == FALSE,]
  
  #Log message
  if(nrow(haps) != nrow(score)){
    stop("From ",nrow(score)," haplotype alleles declared ",nrow(haps)," were found.\n")
  }
  
  #score iterate function
  score.FUN <- function(j){
    x <- haplo$genotypes[haps$IDX,j]
    return(sum(x*haps$SCORE))
  }
  
  #Compute haplotype regression statistics
  a <- mclapply(FUN=score.FUN,X=which(haplo$id.in),mc.cores = ncores)
  out <- NULL
  out$POP <- haplo$pop[haplo$id.in]
  out$ID <- haplo$id[haplo$id.in]
  out$PROFILE <- unlist(a)
  out <- data.frame(out,stringsAsFactors = FALSE)
  
  #Return object
  return(out)
  
}