calc_hwhafsth_FAST <- function(matrix_pol,populations,outgroup=FALSE,simulation=FALSE){



# Imprtant for Coalescent Simulation
if(simulation){
rownames(matrix_pol) <- 1:dim(matrix_pol)[1]
}

# Only one polymorphic site
if(is.vector(matrix_pol)){
   matrix_pol <- as.matrix(matrix_pol)
   warning("#---------------> Only one polymorphic site <-------------------#")
}


if(length(rownames(matrix_pol))==0){
 rownames(matrix_pol) <- 1:dim(matrix_pol)[1]
}


popsize <- length(populations)
npops   <- popsize

init    <- numeric(npops)
init2   <- matrix(0,npops,npops)


matrix_hap <- matrix_pol
# matrix_hap <- as.numeric(matrix_hap)
# matrix_hap <- matrix(matrix_hap,dim(matrix_pol)[1] , dim(matrix_pol)[2])


# If gaps in every site # matrix_hap is empty
if(dim(matrix_hap)[2]==0){

 return(list(hapw=0, fsthALL=NaN, PIW_nei=0,sfreqh= as.matrix(NaN) ))

}
#-------------------------------------------

#### INIT hapwvek ######
hapwvek   <- init # define hpwvek
nam       <- paste("pop", 1:npops)
#######################
nh        <- vector(,npops)
names(nh) <- nam


## Calculate Frequencies ###################################
############################################################

rownames(matrix_hap) <- rownames(matrix_pol)
matrix_hap_sub       <- matrix_hap[unique(unlist(populations)),,drop=FALSE]#### Wegen sfreq !!!!


# uniquematrix         <- unique(matrix_hap_sub)
 duplids               <- .Call("my_unique_C", matrix_hap_sub)
 uniquematrix          <- matrix_hap_sub[!duplids,,drop=FALSE]

nhgesamt             <- dim(uniquematrix)[1]
sfreqh               <- matrix(0,npops,nhgesamt)
 
rownames(sfreqh)     <- nam
colnames(sfreqh)     <- rownames(uniquematrix)
rownames(matrix_hap) <- NULL



 for(xx in 1:npops){

  #sfreqh[xx,] <- get_sfreqh(uniquematrix,matrix_hap[populations[[xx]],,drop=FALSE])  
  sfreqh[xx,]  <- .Call("C_get_sfreqh_C",uniquematrix,matrix_hap[populations[[xx]],,drop=FALSE]) 
  nh[xx]       <- length(which(sfreqh[xx,]!=0)) # number of haplotypes for each population
 
 }


if(nhgesamt==1){

return(list(hapw=0,fsthALL=NaN,PIW_nei=0,sfreqh=sfreqh))

}



# Create happairs ! perhaps there is a better internal R - Function !!!

#if(npops > 1){

# happairsbetween <- NULL
# vek             <- 1:dim(uniquematrix)[1]

# for(xx in 1:dim(uniquematrix)[1]){

# com <- vek[-xx]
# val <- rep(xx,length(com))
# mat <- rbind(val,com)

# happairsbetween <- cbind(happairsbetween,mat) 
# }

#}

# happairswithin      <- combn(dim(uniquematrix)[1],2)

#-------------------------------------------------------------------------


for(xx in 1:npops){
  
  #   div     <- apply(happairswithin,2,function(hap){
  #              freq       <- sfreqh[xx,]
  #              return(freq[hap[1]]*freq[hap[2]])
  #   }) 
                     
  # hapwvek[xx] <- sum(div,na.rm=TRUE)

  hapwvek[xx] <- .Call("combnsum_C", sfreqh[xx,,drop=FALSE]) 
  p_size      <- length(populations[[xx]])
  

  if(p_size>1){   
     hapwvek[xx]     <- hapwvek[xx]/((p_size*(p_size-1))/2)
 }
}


# Do the same for population pairs

if(npops > 1){

 pairs <- combn(npops,2)
 #--------------

 hapamatrix   <- matrix(NaN,npops,npops)

 # Apply
 hapavek     <- apply(pairs,2,function(x){
     
      hapa <- NaN
     
      m1_size <- length(populations[[ x[1] ]]) # size of population 1
      m2_size <- length(populations[[ x[2] ]]) # size of population 2     
     
     freqpop1 <- sfreqh[x[1],,drop=FALSE]
     freqpop2 <- sfreqh[x[2],,drop=FALSE]
     
          #  div     <-  apply(happairsbetween,2,function(hap){
          #              return(freqpop1[hap[1]]*freqpop2[hap[2]])
          #  })
                       
     #hapa  <- sum(div)
     hapa  <- .Call("combnsum2_C",freqpop1,freqpop2)
         
     hapa  <- hapa/(m1_size*m2_size)
     hapamatrix[x[1],x[2]]  <<- hapa
  
     return(hapa)

}) 

hapavek           <- as.matrix(hapavek)

}else{hapavek <- as.matrix(NaN);hapamatrix <- NaN} # End of if (popsize > 1)


# FST HAPLOTYPE # ---------------------------------------------------------- 

if(npops > 1){

poppairs   <- combn(npops,2)   

   shapw  <- sum(hapwvek[1:npops],na.rm=TRUE)
   ncw    <- npops
   
   
   sh   <- apply(poppairs,2,function(x){
         pop1 <- x[1]
         pop2 <- x[2]
         return(hapamatrix[pop1,pop2])      
   })
   
   shapa  <- sum(sh,na.rm=TRUE)
   nca    <- length(sh)
   
   if(shapa){
      fsthALL <- 1 - (shapw/ncw)/(shapa/nca)
   }else{fsthALL <- NaN}

}else{fsthALL <- NaN} # End if(only_populations >1)

###################### END FSTHALL #########################################
  

## calcPi (within) #########################################################

PIW_nei              <- rep(NaN,npops)

# start comment -----
# to slow at the moment for whole genome FIXME

# rownames(matrix_hap) <- rownames(matrix_pol)

#for(xx in 1:npops){
   
#   freqq    <- sfreqh[xx, ]/length(populations[[xx]]) 
#   freqq    <- freqq[which(freqq!=0)]

   
#   if(nh[xx]>1){vergl <- combn(nh[xx],2)}else{PIW_nei[xx]<-0;next;}   
   
#   res      <- apply(vergl,2,function(x){
#   hap1     <- matrix_hap[names(freqq[x[1]]),]
#   hap2     <- matrix_hap[names(freqq[x[2]]),]
#   div      <- hap1!=hap2 # ohne substituitionsmodell
#   div      <- sum(div)
#   return(2*freqq[x[1]]*freqq[x[2]]*div)               
#   })

#n            <- length(populations[[xx]])
#PIW_nei[xx]  <- (n/(n-1))*sum(res,na.rm=TRUE)

#}
# end comment -----


return(list(hapw=hapwvek,fsthALL=fsthALL,sfreqh=sfreqh,PIW_nei=PIW_nei))

} # End of Function hwhafsth_FAST 




