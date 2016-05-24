

HWExact <- function(GenotypeCounts) 
{
      
        nSNP <- nrow(GenotypeCounts)
       
		
        obs.hom1 <- GenotypeCounts$nAA * (GenotypeCounts$nAA >= 
            GenotypeCounts$naa) + GenotypeCounts$naa * (GenotypeCounts$nAA < 
            GenotypeCounts$naa)
			
        obs.hom2 <- GenotypeCounts$nAA * (GenotypeCounts$nAA < 
            GenotypeCounts$naa) + GenotypeCounts$naa * (GenotypeCounts$nAA >= 
            GenotypeCounts$naa)
			
			
			
			
			
        obs.hets <- GenotypeCounts$nAa
        n1 <- 2 * obs.hom1 + obs.hets
        sel.bad <- n1 == 0
        nGood <- nSNP - sum(sel.bad)
        hwep <- rep(NA, times = nSNP)
        hwep.good <- .C("SNPHWE", 
					    PACKAGE = "GWASExactHW",
					    as.integer(nGood),
						as.integer(obs.hets[!sel.bad]), 
						as.integer(obs.hom1[!sel.bad]), 
						as.integer(obs.hom2[!sel.bad]), 
						hwe = double(nGood))$hwe
        hwep[!sel.bad] <- hwep.good
        return(hwep)
}
    
    
    
	
	

    
