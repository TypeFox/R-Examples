GSEforGPL <-
function (con, GPLid) 
{
  connect <- con$connect

  ###########Get idExperiment for the GSE
  if(!is.character(GSEinDB(con)))
    {
      results <- GSEinDB(con)
      i <- which(results$Chip %in% GPLid)
      k <- which(!GPLid %in% results$Chip)
		
      if(length(i)==0){stop(paste("GSE records corresponding to platform record",GPLid,"have not yet been loaded in the compendium database",sep=" "))}
      if(length(k)!=0){cat(paste("GSE records corresponding to platform record",GPLid[k],"have not yet been loaded in the compendium database",sep=" "),"\n")}

      results <- results[i,]
      rownames(results) <- c(1:nrow(results))
      results
    }else{
      stop("No experimental data has been loaded in the compendium database yet")
    }
}

