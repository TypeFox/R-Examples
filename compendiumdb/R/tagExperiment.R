tagExperiment <-
function (con, GSEid, tag)
{
  con <- con$connect
  
  query_GSE <- paste("SELECT expname FROM experiment WHERE expname='",GSEid,"'",sep="")
  rs <- dbSendQuery(con, query_GSE)
  gse <- fetch (rs, n= -1)
  dbClearResult(rs)
  if(nrow(gse)==0){stop(paste("Series record",GSEid,"has not been loaded in the compendium yet",sep=" "))}
	
  query_tag <- paste("UPDATE experiment SET tag='",tag,"' WHERE expname='",GSEid,"'",sep="")

  rs <- dbSendQuery(con, query_tag)
  cat("GSE record",GSEid,"has been tagged \n")
}
