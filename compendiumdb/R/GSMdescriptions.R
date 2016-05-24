GSMdescriptions <-
function (con, GSEid, GPLid = "")
{

  if (length(GSEid)!=1) stop("Please make sure to enter a single GSE ID")
  connect <-con$connect

  ###########Get idExperiment for the GSE
  idGSE <- GSEinDB(con,GSEid)

  idGDS <- idGSE$GDS
  idExperiment <- unique(idGSE$id_Compendium)
  idExperiment_inDB_null <- is.null(idExperiment)

  ## Experimental design as determined by expDesign.pl
  if(GPLid!=""){
    if(length(GPLid)>1){stop("Please make sure to enter a single GPL corresponding to a single GSE")}
    expDesign <- idGSE[idGSE$Chip==GPLid,"experimentDesign"]
    if(length(expDesign)==0){stop(paste(GPLid,"does not correspond to",GSEid))}
  }else{
    expDesign <- unique(idGSE$experimentDesign)
    if(length(expDesign)>1){
      expDesign <- "SC"
    }
  }

  if (!idExperiment_inDB_null)
    {
      query1 <- "SELECT hyb.barcode, hd.hybid, hd.hyb_type, hd.hyb_description from hyb_has_description hd 
				INNER JOIN hyb ON hd.hybid=hyb.hybid 
				INNER JOIN experiment_has_hyb eh ON hyb.hybid=eh.hybid 
				INNER JOIN chip ON hyb.idchip=chip.idchip
				WHERE eh.idExperiment="
      if(GPLid!=""){
        query_hyb <- paste(query1,idExperiment," and chip.db_platform_id='",GPLid,"'",sep="")
      }else{
        query_hyb <- paste(query1,idExperiment,sep="")
      }

      rs <- dbSendQuery(connect, query_hyb)
      data <- fetch (rs, n= -1) 
      dbClearResult(rs)
	
      GSM <- unique(data[,1])
      hybid <- unique(data[,2])
      type <- unique(data[,3])

      flag <- 0		
      if(expDesign=='CR' || expDesign=='DC' || expDesign=='DS'){ ##  two-channel design
        queryS <- "SELECT hyb.barcode,hyb.hybid, sample.sampletitle, sample.samplelabel, sample.samplecharacteristics, sample.samplesource, sample.samplenumber FROM sample
				INNER JOIN hyb_has_sample hs ON sample.idsample=hs.idsample 
				INNER JOIN hyb on hs.hybid=hyb.hybid 
				INNER JOIN chip ON hyb.idchip=chip.idchip
				INNER JOIN experiment_has_hyb eh on hyb.hybid=eh.hybid
				INNER JOIN experiment e on eh.idExperiment= e.idExperiment WHERE e.idExperiment="
			
        query_sample <- paste(queryS,idExperiment,sep="")
        rs <- dbSendQuery(connect, query_sample)
        sampleData <- fetch (rs, n= -1) 
        dbClearResult(rs)

        i <- match(sampleData$barcode,GSM)
        sampleData <- sampleData[which(!is.na(i)),]
				
        labels <- unique(sampleData$samplelabel)
        sampleDesign <- array(NA,dim=c(length(GSM)+1,length(labels)+3))
				
        sampleDesign[1,] <- c(labels,"barcode","samplesource_ch1","samplesource_ch2")
        rownames(sampleDesign) <- c("header",GSM)
				
        sampleD <- apply(sampleDesign,2,function(x){
          y=sampleData[sampleData$samplelabel==x[1],c("barcode","samplecharacteristics")]
          y[match(rownames(sampleDesign),y$barcode),"samplecharacteristics"]
        })

        colnames(sampleD) <- sampleDesign[1,]
        rownames(sampleD) <- rownames(sampleDesign)
        sampleD <- sampleD[-1,]
        sampleD[,"barcode"] <- rownames(sampleD)
        sampleD <- data.frame(sampleD)		

        source1 <- sampleData[sampleData$samplenumber==1,c("barcode","samplesource")]
        source2 <- sampleData[sampleData$samplenumber==2,c("barcode","samplesource")]

        i <- match(sampleD$barcode,source1$barcode)
        sampleD$samplesource_ch1 <- source1$samplesource[i]
        i <- match(sampleD$barcode,source2$barcode)
        sampleD$samplesource_ch2 <- source2$samplesource[i]
        flag <- 1
      }
							
      ###########Create character matrix with descriptions
      first_out <- array(NA,dim=c((length(GSM)+1), (length(type)+3) ) )
      first_out[1,(length(type)+3)]<-"GPL"
      for(i in 1:(length(type)))
        {
          type_one <- as.character(type[i])
          first_out[1,i+2] <- type_one
        }	

      for(j in 1:length(GSM))
        {							
          first_out[j+1,1] <- hybid[j]

          GSMID <- as.character(GSM[j])
          first_out[j+1,2] <- GSMID
          ###########Get platform for the GSM
          platform_query_part1 <- "select distinct chip.db_platform_id from chip inner join hyb on hyb.idchip=chip.idchip where hyb.barcode = '"
          platform_query_part2 <- "'"
          query_platform <- paste(platform_query_part1,GSM[j],platform_query_part2,sep="")
          rs <- dbSendQuery(connect, query_platform)
          platform <- fetch (rs, n= -1)
          platform <- as.character(platform)
          first_out[j+1,(length(type)+3)] <- platform 
        }

      for (k in 1:nrow(data))
        {
          hybID <- data[k,2]
          type_one <- data[k,3]
          description_one <- data[k,4]
          first_out[first_out[,1]==hybID, first_out[1,]== type_one] <- description_one
        }

      first_out[1,1:2] <- c("sampleID","barcode")
      colnames(first_out) <- first_out[1,]
      first_out <- first_out[-1,]
      first_out <- first_out[,-grep("sampleID",colnames(first_out))]
      rownames(first_out) <- first_out[,"barcode"]

      if(flag){ ##Two-channel experiment
        table_for_choice <- merge(sampleD,first_out,by=c("barcode"),all.x=TRUE)
        if(!is.na(match("samplechar",colnames(table_for_choice)))){
          if(unique(!is.na(idGDS))){
            table_for_choice <- table_for_choice[,-grep("samplechar",colnames(table_for_choice))]
          }else{
		indx <- which(!is.na(match(colnames(table_for_choice),c("samplechar","samplesource"))))
		table_for_choice <- table_for_choice[,-indx]
          }
        }
      }else{
        table_for_choice <- first_out
      }				
      rownames(table_for_choice) <- table_for_choice[,"barcode"]
      table_for_choice <- table_for_choice[,-grep("barcode",colnames(table_for_choice))]

      if(GPLid!=""){
        table_for_choice <- table_for_choice[table_for_choice[,"GPL"]==GPLid,]
      }
      as.matrix(table_for_choice)
    }
}

