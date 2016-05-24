pair_linkdisequ <- function(bial1,bial2,populations){


# Create Combinations of sites between the matrices 
# Look for an internal function for that 

combinations    <- NULL
vek             <- 1:dim(bial1)[2]
vek2            <- 1:dim(bial2)[2]

for(xx in vek){

val <- rep(xx,dim(bial2)[2])
mat <- rbind(val,vek2)

combinations <- cbind(combinations,mat) 

}

npops    <- length(populations)
reslist  <- vector("list",npops)

for(yy in 1:npops){

bial1_sub <- bial1[populations[[yy]],,drop=FALSE]
bial2_sub <- bial2[populations[[yy]],,drop=FALSE]

RETURN <- rep(NA,5)

res <- apply(combinations,2,function(x){
	
     site1 <- bial1_sub[,x[1]]
     site2 <- bial2_sub[,x[2]]

     site1x  <- site1[!is.na(site1)]
     site2x  <- site2[!is.na(site2)]
     stAa    <- unique(site1x)
     stBb    <- unique(site2x)
     if(length(stAa)<2|length(stBb)<2){return(RETURN)} # one of the site is not polymorphic
     
     
     # first site 
     #f1 <- stAa[1]
     #s1 <- stAa[2]     
     # site 1 ------------------ Frequenzen 0/1
     #first   <- sum(site1==f1,na.rm=TRUE)
     #second  <- sum(site1==s1,na.rm=TRUE)
     
     ones    <- sum(site1==1,na.rm=TRUE)
     zeros   <- sum(site1==0,na.rm=TRUE)
     
     alle     <- ones + zeros
     
     if(ones>=zeros){freqsite1 <- ones/alle;id <- 1}else{freqsite1 <- zeros/alle;id <- 0} # what is the mutation
     
     # second site
     #f2 <- stBb[1]
     #s2 <- stBb[2]
     # site 2  -----------------
     #first   <- sum(site2==f2,na.rm=TRUE)
     #second  <- sum(site2==s2,na.rm=TRUE)
     
     ones    <- sum(site2==1,na.rm=TRUE)
     zeros   <- sum(site2==0,na.rm=TRUE)
     
     alle     <- ones + zeros
     if(ones>=zeros){freqsite2 <- ones/alle;id2 <- 1}else{freqsite2 <- zeros/alle;id2 <- 0} # what is the mutation

     aa     <- site1==id
     bb     <- site2==id2
     cc     <- which(aa==TRUE & bb==TRUE)
     x      <- length(cc)
     d_raw  <- x/alle - freqsite1*freqsite2 #### d_raw
     # -------------------------------------------
    
     
     freqsite1_low <- 1-freqsite1
     freqsite2_low <- 1-freqsite2
     r2 <- (d_raw*d_raw)/(freqsite1*freqsite1_low*freqsite2*freqsite2_low)
     r <- sqrt(r2)
     if(d_raw<0){x <- min(freqsite1*freqsite2,freqsite1_low*freqsite2_low);r <- -1*r}else{x <- min(freqsite1*freqsite2_low,freqsite1_low*freqsite2)}
     d_prime <- d_raw/x   #### d_prime

     # --------------------------------------- #
    # M      <- matrix(0,2,2)
    # id     <- which(site1==stAa[1] & site2==stBb[1])
    # count  <- length(id)
    # M[1,1] <- count
    # id     <- which(site1==stAa[1] & site2==stBb[2])
    # count  <- length(id)
    # M[1,2] <- count
    # id     <- which(site1==stAa[2] & site2==stBb[1])
    # count  <- length(id)
    # M[2,1] <- count
    # id     <- which(site1==stAa[2] & site2==stBb[2])
    # count  <- length(id)
    # M[2,2] <- count
     # ----------------------------------------- #

     # p=fisherextest(M[1,1],M[1,2],M[2,1],M[2,2])

     # d_dist ----------------------------------
     
     # if(length(colnames(matrix_pol))>1){
     #   d_dist <- abs(matrix_pol_names[X[1]]- matrix_pol_names[X[2]])
     #}else{
     #   d_dist <- NaN
     #}
     
  ## RETURN

  RETURN[1] <- d_raw
  RETURN[2] <- d_prime
  RETURN[3] <- r2
  RETURN[4] <- r
  RETURN[5] <- NA
   
  return(RETURN)
     
})

 rownames(res)     <- c("d_raw","d_prime","r2","r","d_dist")
 colnames(res)     <- apply(combinations,2,function(x){return(paste(x[1],"/",x[2],sep=""))})
 reslist[[yy]]     <- res
 reslist           <- as.matrix(reslist)
 rownames(reslist) <- paste("pop",1:npops)
 colnames(reslist) <- "Linkdisequilibrium"

}# End of iteration over pops

return(reslist)

}
