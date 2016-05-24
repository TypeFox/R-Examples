linkdisequ_FAST <- function(matrix_pol,populations){

### NO NA in DATA

## NOTE: Gaps in the biallelic matrix should be deleted
npops        <- length(populations)
sitelength   <- dim(matrix_pol)[2]
if(sitelength<2){return(list(res=as.matrix(NaN),Zns=NaN,ZA=NaN,ZZ=NaN))}

if(length(colnames(matrix_pol))==0){
 colnames(matrix_pol) <- 1:sitelength
}



#numsitepairs     <- choose(sitelength,2)

# INIT # -----------------------------------#
init          <- vector(,npops)
names         <- paste("pop",1:npops)
# ------------------------------------------#
znsvek        <- init
names(znsvek) <- names

#reslist       <- vector("list",npops) 

 ZA            <- init
 names(ZA)     <- names
 ZZ            <- init
 names(ZZ)     <- names

#-------------------------------------------#

 segsites <- get_segsites_FAST(matrix_pol,populations) # positions of the segsites of each population

# RETURN <- rep(NA,5)

for(xx in 1:npops){

 popmatrix        <- matrix_pol[populations[[xx]],segsites[[xx]],drop=FALSE]
 #segsites_pop     <- .Call("get_segsites_C",popmatrix)
 #segsites.pos     <- which(segsites_pop)
 n.segsites.pop   <- length(segsites[[xx]])
 
 # print(n.segsites.pop)

 if(n.segsites.pop<=1){next}

 # matrix_pol_names <- as.numeric(colnames(matrix_pol[,segsites[[xx]]]))

  # sitepairs        <- combn(n.segsites.pop,2)
 
 #### NAMES ----------------------------------------
#pp <- combn(1:length(segsites[[xx]]),2)
#nn <- paste("s",pp[1,1],"/s",pp[2,1],sep="")
#if(dim(pp)[2]>1){ # more than 2 sites
# for(yy in 2:dim(pp)[2]){
#    m <- paste("s",pp[1,yy],"/s",pp[2,yy],sep="")
#    nn <- c(nn,m)
# }
#}#END if
#### ---------
 
  # site_length  <- dim(popmatrix)[1]
  EINSEN       <- colSums(popmatrix, na.rm = TRUE)
  NULLEN       <- colSums(popmatrix==0, na.rm = TRUE) # include.unknown=TRUE

  res          <- .Call("R2_C", popmatrix, EINSEN, NULLEN)  
 
  # NULLEN       <- site_length - EINSEN

  # BOOL         <- EINSEN>=NULLEN

  #NULLEN       <- site_length - EINSEN
  #ID           <- EINSEN>=site_length
  #ID2          <- NULLEN>=site_length
  #FREQSITE     <- EINSEN[ID]/site_length
  #FREQSITE_LOW <- NULLEN[!ID]/site_length
 
  #Apply
  # res <-  apply(sitepairs,2,function(X){
     
  #    site1   <- popmatrix[,X[1]]
  #    site2   <- popmatrix[,X[2]]
     
     #site1x  <- site1[!is.na(site1)]
     #site2x  <- site2[!is.na(site2)]
     #stAa    <- sum(site1)
     #stBb    <- sum(site2)
     #if(length(stAa)<2|length(stBb)<2){return(RETURN)} # one of the site is not polymorphic
     
     
     # first site 
     #f1 <- stAa[1]
     #s1 <- stAa[2]     
     # site 1 ------------------ Frequenzen 0/1
     #first   <- sum(site1==f1,na.rm=TRUE)
     #second  <- sum(site1==s1,na.rm=TRUE)
     
    # ones     <- EINSEN[X[1]] #sum(site1)
    # zeros    <- NULLEN[X[1]] #site_length - ones
      
     #  if(BOOL[X[1]]){freqsite1 <- ones/site_length;id <- 1}else{freqsite1 <- zeros/site_length;id <- 0}
    #  if(ones>=zeros){freqsite1 <- ones/site_length;id <- 1}else{freqsite1 <- zeros/site_length;id <- 0} # what is the mutation
     
     # second site
     #f2 <- stBb[1]
     #s2 <- stBb[2]
     # site 2  -----------------
     #first   <- sum(site2==f2,na.rm=TRUE)
     #second  <- sum(site2==s2,na.rm=TRUE)
     
    # ones    <- EINSEN[X[2]]  #sum(site2)
    # zeros   <- NULLEN[X[2]]  #site_length - ones
     
     # alle     <- ones + zeros

       #if(BOOL[X[2]]){freqsite2 <- ones/site_length;id2 <- 1}else{freqsite2 <- zeros/site_length;id2 <- 0}
    #  if(ones>=zeros){freqsite2 <- ones/site_length;id2 <- 1}else{freqsite2 <- zeros/site_length;id2 <- 0} # what is the mutation

    # aa     <- site1==id
    # bb     <- site2==id2
    # cc     <- which(aa&bb)
    # x      <- length(cc)
     
     #  x    <- sum(aa&bb)  
    # d_raw  <- x/site_length - freqsite1*freqsite2 #### d_raw
     # -------------------------------------------
    
     
    # freqsite1_low <- 1-freqsite1
    # freqsite2_low <- 1-freqsite2
    # r2 <- (d_raw*d_raw)/(freqsite1*freqsite1_low*freqsite2*freqsite2_low)
  
  #   r <- sqrt(r2)
  #   if(d_raw<0){x <- min(freqsite1*freqsite2,freqsite1_low*freqsite2_low);r <- -1*r}else{x <- min (freqsite1*freqsite2_low,freqsite1_low*freqsite2)}
#     d_prime <- d_raw/x   #### d_prime

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
  #RETURN[1] <- d_raw
  #RETURN[2] <- d_prime
  #RETURN[3] <- r2
  #RETURN[4] <- r
  #RETURN[5] <- d_dist
   
 # return(r2)
     
#})

#rownames(res) <- c("d_raw","d_prime","r2","r","d_dist")
#colnames(res) <- nn
#reslist[[xx]] <- res

# ZNS -------------------------------------------------------
#############################################################
 
 # Zns      <- sum(res[3,],na.rm=TRUE)/numsitepairs
 Zns        <- mean(res, na.rm = TRUE)
 znsvek[xx] <- Zns

# ZA/ZZ (Rozas) ---------------------------------------------
#############################################################

   # adjacent <- apply(sitepairs,2,function(x){return(x[2]==(x[1]+1))})

if(n.segsites.pop==2){

	adjacent <- 1

}else{
	adjacent <- rep(1,n.segsites.pop-1)
	fill     <- n.segsites.pop - 1
	for (yy in 2:length(adjacent)) {
 	adjacent[yy] <- adjacent[yy-1] + fill
 	fill <- fill - 1  
	}

}
   
   ZA[xx]   <- sum(res[adjacent], na.rm = TRUE)/(n.segsites.pop-1)
   ZZ[xx]   <- ZA[xx] - znsvek[xx]
 

}# End of iteration over pops

 #reslist           <- as.matrix(reslist)
 #rownames(reslist) <- paste("pop",1:npops)
 #colnames(reslist) <- "Linkdisequilibrium"
 
 return(list(Zns=znsvek,ZA=ZA,ZZ=ZZ))

}



