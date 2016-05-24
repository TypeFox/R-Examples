
site_FST <- function(bial, populations){

pop.within.div <- matrix(0,length(populations),dim(bial)[2])

for ( xx in 1:length(populations) ){

	popbial             <- bial[populations[[xx]],,drop=FALSE]
	pop.within.div[xx,] <- calc_nuc_diversity_within_per_site(popbial)

}

#print(pop.within.div)

# calc mean 
mean.div         <- apply(pop.within.div,2,mean)
mean.div.between <- calc_average_nuc_diversity_between_per_site(bial, populations)

FST <- 1 - (mean.div/mean.div.between)
FST[is.na(FST)] <- 0

return(FST)

}

## SUBFUNCTIONS
calc_nuc_diversity_within_per_site <- function(matrix){
  
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

  #erg <- sum(erg)

return(erg)
}

calc_average_nuc_diversity_between_per_site <- function(matrix, populations){
  
  n.individuals   <- dim(matrix)[1]
  if(n.individuals==1){return(0)}
  if(length(populations)<2){return(NaN)}
  
  npops       <- length(populations)
  pop.pairs   <- combn(npops,2)
  n.pop.pairs <- choose(npops,2)  

  erg <- apply(matrix,2,function(x){
	
	res <- apply(pop.pairs,2,function(y){
	      pop1vec  <- x[populations[[y[1]]]]
              pop2vec  <- x[populations[[y[2]]]]
              pop1zero <- sum(pop1vec==0, na.rm=TRUE)
              pop2zero <- sum(pop2vec==0, na.rm=TRUE)
              pop1one  <- sum(pop1vec==1, na.rm=TRUE)
              pop2one  <- sum(pop2vec==1, na.rm=TRUE)
	      pop1size <- pop1zero + pop1one
              pop2size <- pop2zero + pop2one
	      vergl    <- pop1size * pop2size
              divers   <- (pop1zero*pop2one + pop1one*pop2zero)/vergl
	      return(divers)	  	
	})
  return(mean(res))
  })  
return(erg)
}



