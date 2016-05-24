calc_phi_st <- function(bial, populations){

# Phi_ST

res       <- calc_pi(bial,populations)
pi        <- res$pi
whole.pi  <- res$whole.pi

if(length(populations)==1){
 return(list(phi_ST=NaN, whole.pi=whole.pi)) 
}

if(is.na(whole.pi)){return(list(phi_ST=1,whole.pi=0))} # only one haplotype !

phi_ST    <- (whole.pi - mean(pi))/whole.pi

return(list(phi_ST=phi_ST, whole.pi=whole.pi))

}
# End of function calc_pi 

# --------------------------------------------------------------------------

# calc_Pi

calc_pi    <- function(bial, populations){


npops      <- length(populations)
wholepop   <- unique(unlist(populations))
popbial    <- bial[wholepop,,drop=FALSE]
uniquebial <- unique(popbial)

# Make Distance Matrix
n.haps     <-  dim(uniquebial)[1]

if(n.haps<=1){return(list(pi=NaN,whole.pi=NaN))}

n.snps     <- dim(uniquebial)[2]
GLxx       <- new.env()
GLxx$DistM <- matrix(,n.haps,n.haps)
happairs   <-  combn(n.haps,2)
apply(happairs,2,function(x){
     hap1                  <-  uniquebial[x[1],]
     hap2                  <-  uniquebial[x[2],]
     GLxx$DistM[x[1],x[2]] <-  sum(hap1!=hap2, na.rm=TRUE)
     })

# Calc Frequencies of Haplotypes in each population
hap.freqs <- matrix(0,npops,n.haps)
counts    <- matrix(0,npops,n.haps)

 for (xx in 1:npops){
  counts[xx,]    <- .Call("C_get_sfreqh_C",uniquebial,bial[populations[[xx]],,drop=FALSE]) 
  hap.freqs[xx,] <- counts[xx,]/length(populations[[xx]]) 
 }
 
 whole.hap.freqs <- colSums(counts)/length(wholepop) # the whole population


if(npops>1){

# Finally calc pi for each pop
 temp <- apply(happairs,2,function(x){
          freq.hap.1 <- hap.freqs[,x[1]]
          freq.hap.2 <- hap.freqs[,x[2]]
          erg        <- 2*freq.hap.1*freq.hap.2*GLxx$DistM[x[1],x[2]]
          return(erg)
         })

 pi       <- rowSums(temp)

}else{ pi <- NaN }

# And calc pi for the whole pop
temp <- apply(happairs,2,function(x){
         freq.hap.1 <- whole.hap.freqs[x[1]]
         freq.hap.2 <- whole.hap.freqs[x[2]]
         erg        <- 2*freq.hap.1*freq.hap.2*GLxx$DistM[x[1],x[2]]
         return(erg)
        })


whole.pi <- sum(temp)


return(list(pi=pi,whole.pi=whole.pi))

}# End of function calc_Pi
