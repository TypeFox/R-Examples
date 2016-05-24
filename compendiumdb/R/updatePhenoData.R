updatePhenoData<-
function(con, GSEid, data)
{
  if (mode(data)!="character") stop("Argument 'data' must be a character matrix")
  if (length(grep("GSM",rownames(data)))!=nrow(data)) stop("All row names of the argument 'data' should be valid GSM IDs")

  # dummies to trick check()
  GPL <- NULL
  Row.names <- NULL 

  pdata <- GSMdescriptions(con,GSEid=GSEid)
 
  commonRows <- intersect(rownames(pdata),rownames(data))
  
  if (length(commonRows)!=length(rownames(pdata))) {
      stop("Row names of the argument 'data' should be identical to those retrieved via the function GSMdescription for the specified GSE record")
  }else {
      response <- readline(paste0("You are going to update the phenotypic data of ",GSEid,". The update will overwrite the phenotypic data currently stored in the compendium database.\nType 'y' if you want to continue, type 'n' if you don't want to update the phenotypic data: "))
      if(response=='n')warning(paste("Phenotypic data of",GSEid,"has not been updated"),call.=FALSE)
  }

  #if (length(commonCol)!=length(colnames(pdata))) stop("Column names of the argument 'data' should be identical to those retrieved via the function GSMdescription for the specified GSE record")

  connect <- con$connect
  barcode <- paste(rownames(data),collapse="','")
  barcode <- sub("^","'",barcode)
  barcode <- sub("$","'",barcode)

  db_query_hybid <- paste("SELECT hybid,barcode FROM hyb WHERE barcode IN(",barcode,")",sep="")
  rs <- dbSendQuery(connect,db_query_hybid)
  hybdata <- fetch(rs,n=-1)
  hybid <- hybdata[,1]
  if (nrow(hybdata)!=nrow(data)) stop("All row names of the argument 'data' should be valid GSM IDs of a GSE record loaded in the compendium database")

  db_delete_description <- paste("DELETE FROM hyb_has_description WHERE hybid IN(",paste(hybid,collapse=","),")",sep="")
  rs <- dbSendQuery(connect,db_delete_description)

  rownames(hybdata) <- hybdata$barcode
  newdata <- as.matrix(merge(hybdata,data.frame(data),by="row.names"))
  newdata <- subset(newdata,select=-c(Row.names,barcode))
  colnames(newdata) <- c("hybid",colnames(data))

  if("GPL" %in% colnames(newdata))
    {
      if (ncol(newdata)==1) stop("At least one of the column names should not be 'GPL'")
      newdata <- subset(newdata,select=-c(GPL))
    }

  for(i in 1:nrow(newdata))
    {
      hybid <- newdata[i,1]
      for(j in 2:ncol(newdata))
        {
          hyb_type <- colnames(newdata)[j]
          hyb_description <- newdata[i,j]
          
          query_insert <- paste("INSERT INTO hyb_has_description (hybid,hyb_type,hyb_description) VALUES(",hybid,",'",hyb_type,"','",hyb_description,"')",sep="")
          dbSendQuery(connect,query_insert)
        }
    }
    
    cat("Phenotypic data of",GSEid,"has been updated \n")
}
