GDSforGSE <-
function (con, GSEid) 
{
  if(length(GSEid)>1){stop("Please make sure to enter single GSE ID at a time")}  

  connect <- con$connect

  ###########Get idExperiment for the GSE
  if(!is.character(GSEinDB(con)))
    {
      results <- GSEinDB(con)
      i <- which(results$Experiment %in% GSEid)
				
      if(length(i)==0){stop(paste("Series record",GSEid,"has not been loaded in the compendium database yet",sep=" "))}

      results <- results[i,]
      j <- which(is.na(results[,"GDS"]))
      noGDS <- unique(results[j,"Experiment"]) ### GSE record without a GDS
      k <- which(!GSEid %in% results$Experiment) ### If GSE record is not yet in the compendium database

      results <- results[which(!is.na(results[,"GDS"])),]			
      if(nrow(results)==0){stop(paste("Series record",GSEid,"does not have a GDS",sep=" "))}
      else if(length(noGDS)!=0){stop(paste("Series record",noGDS,"does not have a GDS",sep=" "))}
      if(length(k)!=0){cat(paste("Series record",GSEid[k],"has not been loaded in the compendium database yet",sep=" "),"\n")}
      rownames(results)=c(1:nrow(results))
      results
    }else{
      stop("No experimental data has been loaded in the compendium database yet")
    }	
}

