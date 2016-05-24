linkdisequ <- function(matrix_pol,populations){



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
reslist       <- vector("list",npops) 
ZA            <- init
names(ZA)     <- names
ZZ            <- init
names(ZZ)     <- names

#-------------------------------------------#

segsites <- get_segsites(matrix_pol,populations) # positions of the segsites of each population

RETURN <- rep(NA,5)

for(xx in 1:npops){


 if(length(segsites[[xx]])<=1){next}
 matrix_pol_names <- as.numeric(colnames(matrix_pol[,segsites[[xx]]]))
 sitepairs        <- combn(length(segsites[[xx]]),2)
 popmatrix        <- matrix_pol[populations[[xx]],segsites[[xx]],drop=FALSE]
 
 #### NAMES ----------------------------------------
pp <- combn(1:length(segsites[[xx]]),2)
nn <- paste("s",pp[1,1],"/s",pp[2,1],sep="")
if(dim(pp)[2]>1){ # more than 2 sites
 for(yy in 2:dim(pp)[2]){
    m <- paste("s",pp[1,yy],"/s",pp[2,yy],sep="")
    nn <- c(nn,m)
 }
}#END if
#### ---------
 
 
 #Apply
  res <-  apply(sitepairs,2,function(X){
     
     site1   <- popmatrix[,X[1]]
     site2   <- popmatrix[,X[2]]
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
     valid_comp <- length(which(!is.na(site1) & !is.na(site2))) # a new FIX for include.unknown=TRUE
     d_raw  <- x/valid_comp - freqsite1*freqsite2 #### d_raw
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
        d_dist <- abs(matrix_pol_names[X[1]]- matrix_pol_names[X[2]])
     #}else{
     #   d_dist <- NaN
     #}
     
  ## RETURN
  RETURN[1] <- d_raw
  RETURN[2] <- d_prime
  RETURN[3] <- r2
  RETURN[4] <- r
  RETURN[5] <- d_dist
   
  return(RETURN)
     
})

rownames(res) <- c("d_raw","d_prime","r2","r","d_dist")
colnames(res) <- nn
reslist[[xx]] <- res

# ZNS -------------------------------------------------------
#############################################################
 
 # Zns        <- sum(res[3,],na.rm=TRUE)/numsitepairs
 Zns        <- mean(res[3,],na.rm=TRUE)
 znsvek[xx] <- Zns

# ZA/ZZ (Rozas) ---------------------------------------------
#############################################################

adjacent <- apply(sitepairs,2,function(x){return(x[2]==(x[1]+1))})
ZA[xx]   <- sum(res[3,adjacent],na.rm=TRUE)/(length(segsites[[xx]])-1)
ZZ[xx]   <- ZA[xx] - znsvek[xx]
 

}# End of iteration over pops

 reslist           <- as.matrix(reslist)
 rownames(reslist) <- paste("pop",1:npops)
 colnames(reslist) <- "Linkdisequilibrium"
 
 return(list(res=reslist,Zns=znsvek,ZA=ZA,ZZ=ZZ))

}



