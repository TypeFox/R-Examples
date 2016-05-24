emp.calc <- function(tab.pop.pop, value="rxy", ref.pop="NA")
{
         
    number.loci <- (ncol(tab.pop.pop)-2)/2
    
    empirical.share.ls <- vector("list",number.loci)
    
    # Calculation of value for each locus in population tab.pop.pop
    for (i in 1:number.loci)
    {
      message(paste("---","Calculations for empirical values are performed for Locus",i,"----",Sys.time()),"\n")
      
      # Empirical share calculated for each locus
      empirical.share.ls[[i]] <- allele.sharing(tab.pop.pop,tab.pop.pop,(i*2)+1,FALSE, value, ref.pop)
      names(empirical.share.ls)[i] <- paste(names(tab.pop.pop)[(i*2)+1],names(tab.pop.pop)[(i*2)+2],sep="-")
      
    }
    
    # Mean empirical
    empirical.list <- do.call("cbind",lapply(empirical.share.ls,array))	
    empirical.list <- matrix(rowMeans(empirical.list,na.rm=TRUE),nrow(empirical.share.ls[[1]]))
    row.names(empirical.list) <- row.names(empirical.share.ls[[1]])
    colnames(empirical.list) <- colnames(empirical.share.ls[[1]])
			
    return(empirical.list)
}
