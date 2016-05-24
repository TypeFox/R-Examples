removeGSE <-
  function(con, GSEid){

   user <- con$user
   password <- con$password
   host <- con$host
   dbname <- con$dbname
   port <- con$port


    query_GSE <- paste("SELECT expname FROM experiment WHERE expname='",GSEid,"'",sep="")
    rs <- dbSendQuery(con$connect, query_GSE)
    gse <- fetch (rs, n= -1)
    dbClearResult(rs)
    if(nrow(gse)==0){stop(paste("Series record",GSEid,"has not been loaded in the compendium database yet",sep=" "))}

    GPLid <- GSEinDB(con, GSEid)$Chip
      
    dir <- path.package("compendiumdb")

    plFile <- paste(dir,"/scripts/Perl/deleteAllforGSE.pl",sep="")
    plFile <- gsub("^","\"",plFile)
    plFile <- gsub("$","\"",plFile) 

    system(paste("perl",plFile,GSEid,user,password,host,port,dbname))

    ## remove the corresponding GPLs if there are no experiments with this GPL ID left
    plFile <- paste(dir,"/scripts/Perl/deleteAllforGPL.pl",sep="")
    plFile <- gsub("^","\"",plFile)
    plFile <- gsub("$","\"",plFile)

    newerror<- tryCatch(
    		GSEinDB(con),
       	error = function(e) e
    )

    for (id in GPLid){
      if (is.null(dim(newerror)) || !(id %in% newerror$Chip)){
       system(paste("perl",plFile,id,user,password,host,port,dbname))
     }
    }
}


