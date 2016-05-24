#Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
#Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127
#modified

complike <- function(bial,populations, freq.table){

# freq.table is a list 

npops    <- length(populations)
CL       <- rep(NaN,npops)
CLR      <- rep(NaN,npops)

freqs    <- jointfreqdist(bial,populations,outgroup=FALSE)[[1]]

for(xx in 1:npops){
 
  # Calculate logCL
  # CL[xx]    <- sum(log(freqs[xx,][freqs[xx,]>0]))
  
  Counts1     <- table(freqs[xx,])
  p1          <- as.numeric(names(Counts1))
  val         <- (Counts1/sum(Counts1))^Counts1 # p^Counts  # Counts^p 
  CL[xx]      <- prod(val[val>0]) #sum(log(val[val>0])) 

	# CLR 
	if(is.list(freq.table)){
	 p2        <- as.numeric(names(freq.table[[xx]]))
         pids      <- match(p1,p2)
         Counts2   <- freq.table[[xx]][pids]
         val       <- (Counts2/sum(freq.table[[xx]]))^Counts1
         val       <- prod(val[val>0])
 	 CLR[xx]   <- 2*( log(CL[xx]) - log(val) )   #sum(log(val[val>0])))    	
	}

 
}

return(list(CL=CL,CLR=CLR))

}

complike_FST <- function(bial, populations, freq.table){

# freq.table is a list 
freq.table <- list(freq.table)
#npops     <- length(populations)
npops    <- 1 # because FST is only one value 
CL       <- rep(NaN,npops)
CLR      <- rep(NaN,npops)

freqs    <- site_FST(bial,populations)

for(xx in 1:npops){
 
  # Calculate logCL
  # CL[xx]    <- sum(log(freqs[xx,][freqs[xx,]>0]))
  
  Counts1     <- table(freqs)
  p1          <- as.numeric(names(Counts1))
  val         <- (Counts1/sum(Counts1))^Counts1 # p^Counts  # Counts^p 
  CL[xx]      <- prod(val[val>0]) #sum(log(val[val>0])) 

	# CLR 
	if(is.list(freq.table)){
	 p2        <- as.numeric(names(freq.table[[xx]]))
         pids      <- match(p1,p2)
         Counts2   <- freq.table[[xx]][pids]
         val       <- (Counts2/sum(freq.table[[xx]]))^Counts1
         val       <- prod(val[val>0])
 	 CLR[xx]   <- 2*( log(CL[xx]) - log(val) )   #sum(log(val[val>0])))    	
	}

 
}

return(list(CL=CL,CLR=CLR))

}
