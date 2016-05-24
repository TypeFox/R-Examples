CR_gen_geno <-
function(map, map.function, ped, seed, transpos, transval)
{
  genome <- combine_chr(map, map.function)
  genome2 <- genome
  if (transval>0) 
	  genome2[transpos] <- 0

  num.loci <- length(genome)+1

  # select row numbers of founders in pedigree
  indx <- which(ped[,2]==0 & ped[,3]==0)
  n.founders <- length(indx)

  genC <- .C("gengeno", as.double(genome), as.double(genome2), as.integer(ped[,1]), as.integer(ped[,2]), as.integer(ped[,3]), as.integer(num.loci), as.integer(n.founders), as.integer(nrow(ped)), as.integer(seed), as.integer(transpos[1]-1), as.integer(transval-1), out=integer(length=nrow(ped)*(num.loci*2)), PACKAGE="mpMap") 
  
  # convert genC results into sim.data
  sim.data <- matrix(data=genC$out, nrow=nrow(ped), ncol=2*num.loci, byrow=TRUE)
  colnames(sim.data) <- c(paste("ibd1_", names(unlist(map)), sep=""), 
	paste("ibd2_", names(unlist(map)), sep=""))

  geno <- list()
  # make sure that sim.data has the correct elements - founders, finals, etc.
  geno$founders <- sim.data[1:n.founders,1:ncol(sim.data)]
  geno$finals <- sim.data[ped[,4]==1, 1:ncol(sim.data)]
  geno$id <- ped[ped[,4]==1,1] 
  geno$fid <- ped[indx,1]
  geno$pedigree <- ped
  geno$map <- map

  return(geno)
}

