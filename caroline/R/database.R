dbWriteTable2 <- function(con, table.name, df, fill.null = TRUE, add.id=TRUE, row.names=FALSE, pg.update.seq=FALSE, ...){
  fields <- dbListFields(con, table.name)
  fields <- fields[!grepl('\\.\\.pg\\.dropped',fields)]
  
  ## add id column if missing
  if(add.id){
    last.id.list <- dbGetQuery(con, paste("SELECT id FROM", table.name,"ORDER BY id DESC LIMIT 1"))
    if(length(last.id.list)==0)
      n <- 0
    else
      n <- last.id.list[[1]]
    df$id <- 1:nrow(df) + n
  }
  ## look for unloadable columns in the df
  names(df) <- tolower(names(df))
  names(df) <- gsub("\\.",'_',names(df))
  
  clmn.match <- match(names(df), fields)
  if(any(is.na(clmn.match)))
    warning(paste("Found '",names(df)[is.na(clmn.match)], "' not in fields of '", table.name,"' table. Omiting.\n", sep=''))
  
  ## add missing fields to df
  field.match <- match(fields, names(df))
  if(sum(is.na(field.match))>0 & fill.null == TRUE){
    message("creating NAs/NULLs for for fields of table that are missing in your df")
    nl <- as.list(rep(NA, sum(is.na(field.match))))
    df.nms.orgnl <- names(df)
    df <- cbind(df, nl)
    names(df) <- c(df.nms.orgnl, fields[is.na(field.match)])
  } 	
  
  ## reorder df columns as per field order
  reordered.names <- names(df)[match(fields, names(df))]
  if(any(is.na(reordered.names)))
    stop('Too many unmatched columns to database column list. Stopping')
  df <- df[ ,reordered.names]

  
  ## BEGIN ERROR CHECKING
  r <- dbSendQuery(con, paste("SELECT * FROM", table.name,"ORDER BY id DESC LIMIT 1"))
  db.col.info <- dbColumnInfo(r); rownames(db.col.info) <- db.col.info$name

  ## check for na's which might prevent a load
  null.OK <- nv(db.col.info,'nullOK')
  reqd.fields <- names(null.OK[!null.OK])
  na.cols <- sapply(df, function(x) any(is.na(x)) )
  req.miss <- na.cols[reqd.fields]
  if(any(req.miss))
    stop(paste("Didn't load df because required field(s)", paste(names(req.miss)[req.miss],collapse=', '),"contained missing values"))
  
  ## check for length mismatches    
  db.precisions <- nv(db.col.info, 'precision')
  df.nchars <- sapply(df, function(c) max(nchar(c)))
  prec.reqd <- db.precisions > 0
  too.long <- db.precisions[prec.reqd] < df.nchars[prec.reqd]
  if(any(too.long))
    stop(paste("Didn't load df because fields", paste(names(df.nchars)[prec.reqd][too.long],collapse=', '),'were too long'))

  ## check for type mismatches
  db.sclasses <- nv(db.col.info,'Sclass')
  df.classes <- sapply(df, class)
  type.mismatches <- names(df.classes)[db.sclasses != df.classes & !na.cols]
  #if(length(type.mismatches)>0)
  #  warning(paste('The dataframe columns:',paste(type.mismatches, collapse=','),'may have type mismatches from their sclass mappings to the database table fields.'))
    
  dbClearResult(r)
  
  ## check unique constrains
  #r <- dbGetQuery(con, paste("SELECT constraint_name FROM information_schema.table_constraints WHERE table_name = '",table.name,"'", sep=''))
  #if(nrow(df) != length(unique(apply(df[,c('day','file')],1, paste, collapse='.')))
  #stop
  
  ## load table
  print(paste("loading", table.name, "table to database"))
  db.write <- dbWriteTable(con, table.name, df, row.names=row.names, ...)
  
  #updating postgresql sequence
  if(pg.update.seq){
    if(class(con)=='PostgreSQLConnection'){
      r <- dbSendQuery(con, paste("SELECT pg_catalog.setval(pg_get_serial_sequence('",table.name,"', 'id'), (SELECT MAX(id) FROM ",table.name,")+1);",sep=''))
      dbClearResult(r)
    }else{
      stop('pg.update.seq=TRUE flag not compatable with database connection type')
    }
  }
  
  if(db.write & add.id)
    invisible(df$id)
  else
    return(db.write)
}
