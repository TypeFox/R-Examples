# TODO: TV naming of the function (BH while all sorts of tests)
IsoTestBH <- function (rp, FDR, type =  c("BH", "BY"),
          stat = c("E2", "Williams", "Marcus", "M", "ModifM")){

  type <- match.arg(type)
  stat <- match.arg(stat)
        
  Probe.ID <- rp[,1]
    
  rpraw <- switch(stat,
      E2 = rp[,2],
      Williams = rp[,3],
      Marcus = rp[,4],
      M = rp[,5],
      ModifM = rp[,6])
  
 ## procs <- c("BH", "BY") # only "BH" and "BY" are needed
 ## res <- mt.rawp2adjp(rpraw, procs) # from multtest
 ## adjp <- res$adjp[order(res$index), ]
  
    adjp <- cbind(rpraw, p.adjust(rpraw, "BH"),p.adjust(rpraw,"BY"))
  

  # TODO: TV use names Bonferroni   Holm Hochberg     SidakSS     SidakSD          BH
  place.keep33 <- if (type == "BH"){
    which(adjp[,2] <= FDR)
  } else { # type == "BY"  
    which(adjp[,3] <= FDR)
  }

  sign.Probe.ID <- Probe.ID[place.keep33] 
  
  if (type == "BH")  {
    sign.genes <- data.frame(sign.Probe.ID,
                            place.keep33,
                            adjp[adjp[,2] <= FDR,1],
                            adjp[adjp[,2] <= FDR,2])
  } else {
    sign.genes <- data.frame(sign.Probe.ID,
                            place.keep33,
                            adjp[adjp[,3] <= FDR,1],
                            adjp[adjp[,3] <= FDR,3])     
  }    

  names(sign.genes) <- c("Probe.ID", "row.name", "raw p-values",
                         paste(type, "adjusted p values", sep = " "))
  return(sign.genes)
}
