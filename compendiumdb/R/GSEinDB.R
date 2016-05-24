GSEinDB <-
  function (con, GSEid = NULL) 
  {
    Sys.setenv(CYGWIN="nodosfilewarning")
    subGSEinDB <- function(con,GSEid=NULL)
    {      
      connect <- con$connect
    
      query_GSM_inDB <- paste("SELECT eh.idExperiment as id_Compendium,count(h.hybid) as Samples, chip.db_platform_id as Chip, h.expdesign as experimentDesign FROM experiment_has_hyb eh 
                              INNER JOIN hyb h ON (eh.hybid =h.hybid) 
                              INNER JOIN chip ON h.idchip=chip.idchip 
                              GROUP BY  eh.idExperiment, h.expdesign, h.idchip ORDER BY eh.idExperiment",sep="")
      
      query_GSEquery_inDB <- paste("SELECT distinct e.idExperiment as id_Compendium, e.expname as Experiment, e.tag as Tag, chip.db_platform_id as Chip, org.ncbiorgid as OrganismNCBIid, org.officialname as OrganismName, e.date_loaded, hyb.expdesign as experimentDesign  
                                   FROM experiment e 
                                   JOIN experiment_has_hyb eh ON e.idExperiment=eh.idExperiment 
                                   JOIN hyb ON eh.hybid=hyb.hybid 
                                   JOIN chip ON hyb.idchip=chip.idchip 
                                   JOIN organism org ON chip.idorganism=org.idorganism ORDER BY e.idExperiment,hyb.expdesign",sep="")
      
      query_GDS_inDB <- paste("SELECT distinct e.idExperiment as id_Compendium, e.expname as Experiment, e.tag as Tag, chip.db_platform_id as Chip, gds.GDS 
                              FROM experiment e 
                              JOIN experiment_has_hyb eh ON e.idExperiment=eh.idExperiment 
                              JOIN hyb ON eh.hybid=hyb.hybid
                              JOIN gds on hyb.idGDS=gds.idGDS 
                              JOIN chip ON hyb.idchip=chip.idchip order by e.idExperiment",sep="")
      
      
      rs <- dbSendQuery(connect, query_GSM_inDB)
      samples_inDB <- fetch (rs, n= -1)
      dbClearResult(rs)
      if(!nrow(samples_inDB)){
        stop("The compendium database is empty: no GSE record has been loaded yet")
      }
      
      rs <- dbSendQuery(connect, query_GSEquery_inDB)
      idExperiment_inDB <- fetch (rs, n= -1)
      #Date <- idExperiment_inDB["date_loaded"]
      #idExperiment_inDB <- idExperiment_inDB[-grep("date_loaded",colnames(idExperiment_inDB))]
      dbClearResult(rs)
      
      rs <- dbSendQuery(connect, query_GDS_inDB)
      gds_inDB <- fetch (rs, n= -1)
      dbClearResult(rs)
      
      results <- merge(idExperiment_inDB,samples_inDB,by=c("id_Compendium","Chip"),all.x=TRUE)
      
      idExperiment_inDB_na<-is.na(idExperiment_inDB)      
      reshead <- c("id_Compendium","Experiment","experimentDesign","Chip","Samples","Tag","OrganismNCBIid","OrganismName","DateLoaded","GDS")
      if(nrow(gds_inDB)!=0){
        y <- merge(results,gds_inDB,by=c("Experiment","Chip"),all.x=TRUE)
        y <- y[order(y$id_Compendium.x),c("id_Compendium.x","Experiment","experimentDesign.x","Chip","Samples","Tag.x","OrganismNCBIid","OrganismName","date_loaded","GDS")]
        y <- unique(y) 
	 rownames(y) <- c(1:nrow(y))
        results <- y
      }else{
        results <- results[,c("id_Compendium","Experiment","experimentDesign.x","Chip","Samples","Tag","OrganismNCBIid","OrganismName","date_loaded")]
        results <- cbind(results,NA)
      }
      
      colnames(results) <- reshead
       
      temp <- array(NA,dim=c(1,ncol(results)))
      colnames(temp) <- colnames(results)	
      
      noGSE=character()
      if(length(GSEid)!=0){
        for(id in GSEid){
          if(id %in% results[,"Experiment"]){ # For exact matching 
            temp <- rbind(temp,results[results$Experiment==id,])
          }else{ 
            noGSE <- c(noGSE,id)
          }
        }
        finalResult <- temp[-1,]
        if(length(noGSE)!=0){stop(c(noGSE," has not been loaded in the compendium database yet"))}

        if(nrow(finalResult)){
          rownames(finalResult) <- c(1:nrow(finalResult))
          return(finalResult)
        }        
      }else{
        if((is.null(idExperiment_inDB_na))){
          stop("The compendium database is empty: no GSE record has been loaded yet")
        }else{ return(results)} 
      }
    }
    
    possibleError <- tryCatch(
      subGSEinDB(con,GSEid=GSEid),
      error = function(e) e
    )
    
    if(length(grep("Table .* doesn't exist",as.character(possibleError)))){
      stop("Please check if the database schema has already been loaded. You can load the schema using the function loadDatabaseSchema")
    }else if(length(grep("compendium database is empty",as.character(possibleError)))){
      stop("The compendium database does not have any GSE records loaded. Load the data using the function loadDataToCompendium")
    }else if(length(grep("not been loaded in the compendium database",as.character(possibleError)))){
      stop(paste(possibleError,collapse=" "))
    }else {
      return(possibleError)
    }
  }
