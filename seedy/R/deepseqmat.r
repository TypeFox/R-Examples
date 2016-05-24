deepseqmat <- function(X) {
      if (!"obs.strain"%in%names(X)) {
        stop("Input must be complete genomic sampling - set 'full' to TRUE in simulation functions")
      }
      observed.strains <- unique(unlist(X$obs.strain)) # unique observed strains
      polymorphic.loci <- NULL
      for (i in 1:length(observed.strains)) {
        polymorphic.loci <- c(polymorphic.loci, X$libr[[which(X$librstrains==observed.strains[i])]])
      }
      # polymorphisms observed over all samples
      polymorphic.loci <- unique(polymorphic.loci) 
      if (sum(is.na(polymorphic.loci)>0)) {
      	polymorphic.loci <- polymorphic.loci[-is.na(polymorphic.loci)]
      }
      
      n <- length(X$obs.freq)
      host.polys <- list() # polymorphisms by host
      host.poly.freq <- list() # frequencies for each polymorphism
      for (i in 1:n) {
        hoststrains <- unlist(X$obs.strain[[i]])
        hostfreq <- unlist(X$obs.freq[[i]])
        hp <- NULL
        hf <- NULL
        for (j in 1:length(hoststrains)) {
          hp <- c(hp, X$libr[[which(X$librstrains==hoststrains[j])]])
          hf <- c(hf, rep(hostfreq[j],length(X$libr[[which(X$librstrains==hoststrains[j])]])))
        }
        if (sum(is.na(hp)>0)) {
          hp <- hp[-which(is.na(hp))]
        }
        host.polys[[i]] <- unique(hp)
        hf1 <- NULL
        for (j in 1:length(host.polys[[i]])) {
          hf1 <- c(hf1, sum(hf[which(hp==host.polys[[i]][j])]))
        }
        host.poly.freq[[i]] <- hf1/sum(X$obs.freq[[i]])
      }
      
      polytable <- matrix(0, length(polymorphic.loci), n)
      
      for (i in 1:n) {
        for (j in 1:length(host.polys[[i]])) {
          polytable[which(polymorphic.loci==host.polys[[i]][j]),i] <- host.poly.freq[[i]][j]
        }
      }
      rownames(polytable) <- paste("locus", polymorphic.loci)
      if ("sampledata"%in%names(X)) {
         colnames(polytable) <- paste("patient", X$sampledata[,1])     	
      } else {
         colnames(polytable) <- paste("sample", 1:n)      		
	  }
      return(polytable)
}


