calc_diversities <- function(matrix_pol,populations,pi=FALSE){

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
nucwvek   <- init
PIW_nei   <- init
nam       <- paste("pop", 1:npops)
#######################
nh        <- vector(,npops)
names(nh) <- nam


## Calculate Frequencies ###################################
############################################################

rownames(matrix_hap) <- rownames(matrix_pol)
matrix_hap_sub       <- matrix_hap[unique(unlist(populations)),,drop=FALSE]#### Wegen sfreq !!!!


# uniquematrix       <- unique(matrix_hap_sub)

duplids             <- .Call("my_unique_C", matrix_hap_sub)
uniquematrix        <- matrix_hap_sub[!duplids,,drop=FALSE]

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

return(list(hapwvek=hapwvek,nucw=nucwvek,PIW_nei=PIW_nei,sfreqh=sfreqh))

}

# happairswithin      <- combn(dim(uniquematrix)[1],2)

#-------------------------------------------------------------------------
# calculate haplotype diversity

for(xx in 1:npops){
  
#     div     <- apply(happairswithin,2,function(hap){
#                freq       <- sfreqh[xx,]
#                return(freq[hap[1]]*freq[hap[2]])
#     }) 
                     
#  hapwvek[xx] <- sum(div,na.rm=TRUE)
  hapwvek[xx] <- .Call("combnsum_C", sfreqh[xx,,drop=FALSE])  
  p_size      <- length(populations[[xx]])
  

  if(p_size>1){   
     hapwvek[xx]     <- hapwvek[xx]/((p_size*(p_size-1))/2)
 }
}

#---------------------------------------------------------
# calculate nucleotide diversity

for(xx in 1:npops){

popmat      <- matrix_pol[populations[[xx]],,drop=FALSE]
nucwvek[xx] <- calc_nuc_diversity_within(popmat)

}


## calcPi (within) #########################################################


# start comment -----
# to slow at the moment for whole genome FIXME

if(pi){

 rownames(matrix_hap) <- rownames(matrix_pol)

for(xx in 1:npops){
   
   freqq    <- sfreqh[xx, ]/length(populations[[xx]]) 
   freqq    <- freqq[which(freqq!=0)]

   
   if(nh[xx]>1){vergl <- combn(nh[xx],2)}else{PIW_nei[xx]<-0;next;}   
   
   res      <- apply(vergl,2,function(x){
   hap1     <- matrix_hap[names(freqq[x[1]]),]
   hap2     <- matrix_hap[names(freqq[x[2]]),]
   div      <- hap1!=hap2 # ohne substituitionsmodell
   div      <- sum(div, na.rm=TRUE)
   return(2*freqq[x[1]]*freqq[x[2]]*div)               
   })
n            <- length(populations[[xx]])
PIW_nei[xx]  <- (n/(n-1))*sum(res,na.rm=TRUE)

}

}# if pi
# end comment -----


return(list(hapw=hapwvek,nucw=nucwvek,sfreqh=sfreqh,PIW_nei=PIW_nei))

} # End of Function calc_diversities

# SUBFUNCTIONS
calc_nuc_diversity_within <- function(matrix){
  
  n.individuals   <- dim(matrix)[1]

  if(n.individuals==1){return(0)}

  #n.vergleiche    <- choose(dim(matrix)[1],2)
  
  erg <- apply(matrix,2,function(x){
             
         einsen       <- sum(x==1, na.rm=TRUE)
         nullen       <- sum(x==0, na.rm=TRUE)
	 all          <- einsen + nullen
         n.vergleiche <- (all*(all-1))/2  
         if(n.vergleiche==0){return(0)}                                   
         div          <- (einsen * nullen)/n.vergleiche 
         return(div)
         })

  erg <- sum(erg)

return(erg)
}

