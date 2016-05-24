compute_intern_achaz <- function(bial,populations,outgroup){

npops        <- length(populations)
Achaz.value  <- rep(NaN,npops)

if(outgroup[1]==FALSE){

for (xx in 1:npops ){

  pop.matrix <- bial[populations[[xx]],,drop=FALSE]
  FREQ_ERG   <- .Call("compute_FREQ_C",pop.matrix)
  FREQ_ERG   <- rowSums(FREQ_ERG)
   
  init    <- vector(,floor(length(populations[[xx]])/2))
  w1      <- init
  w2      <- init
  sfreqn  <- init  
  popsize <- length(populations[[xx]])
  sfreq   <- FREQ_ERG[2:popsize] # monomorph positions are not relevant # ACHTUNG bei sehr kleinen Populationen !!!
  
  for(x in 1:floor(popsize/2)) {
				if(x == 1) {w1[x] <- 0
				}else{
					if(x == popsize-x){w1[x] <- popsize/2.0
					}else{w1[x] <- popsize}
				}
				if(x == 1){w2[x] <- 0
				}else{
					if(x == popsize-x){ w2[x] <- popsize/(x*(popsize - x)*2.0)
					}else{ w2[x] <- popsize/(x*(popsize - x))}
				}
				if(x == popsize-x){ sfreqn[x] <- as.integer(((sfreq[x]+ sfreq[popsize-x])/2))
				}else{sfreqn[x] <- sfreq[x] + sfreq[popsize-x]
				}
  }

Achaz.value[xx] <- freqtestn_achaz(popsize,sfreqn,0,w1,w2)
# print(Achaz.value)

}# End over all pops

return(Achaz.value)

}# End no Outgroup

#------------------------------------- WITH OUTGROUP !!!

if(outgroup[1]!=FALSE){

   
  for (xx in 1:npops ){

    popsize    <- length(populations[[xx]]) 
    pop.matrix <- bial[populations[[xx]],,drop=FALSE]
    out.matrix <- bial[outgroup,,drop=FALSE]
    
    # Get ANCESTRAL  
    ancestral  <- .Call("verify_ancestral_C",out.matrix)
    #------------ 
    all.matrix <- rbind(pop.matrix,ancestral)
    FREQ_ERG   <- .Call("compute_FREQOUT_C",all.matrix)
    FREQ_ERG   <- rowSums(FREQ_ERG)

    w1 <- c(0, seq(length(populations[[xx]])-2,1)) # ?? 0 3 2 1 # Population size 5
    w2 <- c(0,1/(3:length(populations[[xx]])-1))   # ?? 0 1/2 1/3 1/4 # Population size 5
    # Yach[xx]   <- freqtesto_achaz(length(populations[[xx]]),FREQ_ERG[[xx]],THETA_ERG["S",xx],w1,w2)
    
    Achaz.value[xx] <- freqtesto_achaz(popsize,FREQ_ERG,0,w1,w2)

  }

return(Achaz.value)

}# End of outgroup included






}
