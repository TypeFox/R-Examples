eqtlDir.P <- function(genoGroups,gex,mc=mc,nper){

  output <- c()
    innerFunction <- function(i)
    {
      if(sum(is.na(genoGroups[,i]))<(nrow(genoGroups))){
	PI <- length(table(genoGroups[,i]))
	# All individuals have the same genotype, do nothing
	if(PI==1){
	  output[i] <- 1.25
	
	# Then the 2 groups comparison
	} else if (PI==2){
	  output[i] <- gmw(gex[!is.na(genoGroups[,i])],genoGroups[!is.na(genoGroups[,i]),i],test="mw",type="permut",alternative="two.sided",nper=nper)

	# And the three group comparison
	} else if (PI==3){
	  p1 <- gmw(gex[!is.na(genoGroups[,i])],genoGroups[!is.na(genoGroups[,i]),i],test="triple",type="permutation",alternative="greater",alg="Csub",nper=2000)$p.values
	  p2 <- gmw(gex[!is.na(genoGroups[,i])],createGroups(genoGroups[!is.na(genoGroups[,i]),i],c(2,1,0)),test="triple",type="permutation",alternative="greater",alg="Csub",nper=2000)$p.values
	  output[i] <- min(2*min(p1,p2),1)
	  
	}
      } else {
	output[i] <- 1.5
      }
   }
  output <- unlist(mclapply(1:ncol(genoGroups),innerFunction,mc.cores=mc))
  output
}
