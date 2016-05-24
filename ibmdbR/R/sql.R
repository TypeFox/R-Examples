# 
# Copyright (c) 2010, 2014, IBM Corp. All rights reserved. 
# 		
# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version. 
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details. 
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>. 
# 
# 

utils::globalVariables(c("idaRGlobal"))

################ show tables ############################
idaShowTables <- function(showAll=FALSE, matchStr=NULL) {
  idaCheckConnection();
  if (!is.null(matchStr) && !is.character(matchStr))
    stop("This function can only be applied to character values.")
    
  wherePart <- "WHERE (OWNERTYPE = 'U')";
  if(!showAll) {
    currSchema <- idaGetCurrentSchema();	
    wherePart <- paste(wherePart," AND (TABSCHEMA= '",currSchema,"') ",sep='');
  }
  
  if (!is.null(matchStr))
    wherePart <- paste(wherePart," AND TABNAME LIKE '%", paste(matchStr,collapse='%',sep=''), "%'", sep='')
  
  idaQuery('SELECT distinct TABSCHEMA as "Schema", TABNAME as "Name", OWNER as "Owner", TYPE as "Type" from SYSCAT.TABLES ', wherePart ,'ORDER BY "Schema","Name"')
}

################ Internal schema and catalog utilities ############################
idaGetCurrentSchema <- function() {
  return(idaScalarQuery("select trim(current_schema) from sysibm.sysdummy1"));
}

prepareCatalogLookup <- function(tableName) {
  x <- parseTableName(tableName);
  return(paste("TABNAME='",x$table,"'",ifelse(!is.null(x$schema),paste(" AND TABSCHEMA = '",x$schema,"'",sep=''),''),sep=''));
}

parseTableName <- function(tableName,removeQuotes=T,disallowQuotes=F) {
  
  inTabNameQuotes <- F;
  sepTableHits <- c();
  
  for(i in 1:nchar(tableName)) {
    
    c <- substr(tableName,i,i);
    if(c=="\"")
      inTabNameQuotes <- !inTabNameQuotes;
    
    if((c==".")&&(!inTabNameQuotes))
      sepTableHits <- c(sepTableHits,i);	
  }
  
  if(inTabNameQuotes) {
    stop("Table name not formatted correctly");
  }
  
  #Only the table name, if DB is schema enabled, we need to 
  # read the current schema, as this could change after the
  # nz.data.frame was created
  if(length(sepTableHits)==0) {
    table <- tableName;
    schema <- idaGetCurrentSchema();
    #Table and schema name provided, if db is not schema enabled, ignore schema	
  } else if(length(sepTableHits)==1) {
    schema <- substr(tableName, 1, sepTableHits[1]-1);		
    table <- substr(tableName, sepTableHits[1]+1, nchar(tableName));
  } else {
    stop("Table name ", table, " is not well-formed");
  }
  
  #Remove quotes or set case
  if (nchar(table) > 2 && substr(table, 1, 1) == '"' && substr(table, nchar(table), nchar(table)) == '"') {
    
    if(disallowQuotes)
      stop("Quoted schema and table names are currently not supported for this function.")
    
    if(removeQuotes)
      table <- substr(table,2,nchar(table)-1)
  } else {
    table <- toupper(table);
  }
  
  #Remove quotes or set case
  if(!is.null(schema)) {
    if (nchar(schema) > 2 && substr(schema, 1, 1) == '"' && substr(schema, nchar(schema), nchar(schema)) == '"') {
      if(disallowQuotes)
        stop("Quoted schema and table names are currently not supported for this function.")
      
      if(removeQuotes)	
        schema <- substr(schema,2,nchar(schema)-1)
    } else { 
      schema <- toupper(schema);
    }
  }
  
  return(list(table=table,schema=schema));
}

################ Basic operations on tables ############################
idaExistTable <- function (tableName) {
  idaCheckConnection();
  cmp <- prepareCatalogLookup(tableName)
  tableExists <- as.logical(as.integer(idaScalarQuery("SELECT CAST(COUNT(*) AS INTEGER) FROM SYSCAT.TABLES WHERE ", cmp))>0);
  return(tableExists);
}

idaIsView <- function (tableName) {
  idaCheckConnection();
  cmp <- paste(prepareCatalogLookup(tableName), " AND TYPE = 'V'",sep=''); 
  viewExists <- as.logical(as.integer(idaScalarQuery("SELECT CAST(COUNT(*) AS INTEGER) FROM SYSCAT.TABLES WHERE ", cmp))>0);
  return(viewExists);
}

idaListTableColumns <- function (tableName) {
  idaCheckConnection();
  cmp <- prepareCatalogLookup(tableName)
  return(as.vector(idaQuery("SELECT COLNAME FROM SYSCAT.COLUMNS WHERE ",cmp," ORDER BY COLNO ASC")[[1]]))
}

idaDeleteViewOrTable <- function(tableName) {
  if(idaIsView(tableName)) {
    idaDropView(tableName)
  } else {
    idaDeleteTable(tableName);
  }
}

idaDeleteTable <- function(table) {
  if(inherits(table,"ida.data.frame")) {
    idaDeleteViewOrTable(table@table)
  } else {
    if (idaExistTable(table))
      idaQuery("DROP TABLE ", table)
    invisible()
  }
}

idaDropView <- function(v) {
  idaCheckConnection();
  try({idaQuery("DROP VIEW ", v)});
}

idaDeleteTempTables <- function() {
  idaCheckConnection();
  tempTables <- idaShowTables()
  tempTables <- tempTables[grep('DATA_FRAME_[0-9]+', tempTables$Name),]
  if (nrow(tempTables)==0) {
    cat("No temporary tables or views are found to be deleted.")
    return(invisible(NULL))
  }
  row.names(tempTables) <- 1:nrow(tempTables)
  cat("The temporary tables and views to be deleted are following:\n")
  print(tempTables)
  cat("\nHow do you like to delete the temporary tables and views?\n")
  cat("1) Delete all temporary tables and views without further confirmation.\n")
  cat("2) Ask conformation for deleting each temporary table or view.\n")
  cat("3) Cancel deletion and return.\n\n")
  while(TRUE){
    cat("Type a key to select your choice:\n")
    cat(" 1) 'a'       2) 'b'       3) 'c'\n")
    choice <- readline()
    if (choice=='a' || choice=='b' || choice=='c' )
      break
  }
  if (choice=='a') {
    for (i in 1:nrow(tempTables)) {
      tempTab = tempTables[i,]
      idaDeleteViewOrTable(paste('"', sub('\\s+$', '', tempTab$Schema), '"."', tempTab$Name, '"', sep=''))
    } 
  }
  else if (choice=='b') {
    for (i in 1:nrow(tempTables)) {
      tempTab <- tempTables[i,]
      schema <- sub('\\s+$', '', tempTab$Schema)
      cat("Delete ", ifelse(tempTab$Type=='T', 'table', 'view'), ": ", 
          schema, '.', tempTab$Name, "? yes[y], no[n] or exit[e] :", sep='')
      choice <- readline()
      if (choice=='e') {
        break
      }
      else if (choice=='y') {
        idaDeleteViewOrTable(paste('"', schema, '"."', tempTab$Name, '"', sep=''))
      }
    } 
  }
}

idaMaterialize <- function(idf,tableName) {
  
  v <- idaCreateView(idf)
   
  tryCatch({
        
        odbcSetAutoCommit(get("p_idaConnection",envir=idaRGlobal), autoCommit = FALSE)
        idaQuery("CREATE TABLE ", tableName, " LIKE ",v," NOT LOGGED INITIALLY");
        idaQuery("INSERT INTO ", tableName, " SELECT * FROM ",v);
        
        odbcEndTran(get("p_idaConnection",envir=idaRGlobal), commit = TRUE);
        
      },
      error=function(e){print(e);odbcEndTran(get("p_idaConnection",envir=idaRGlobal), commit = FALSE)},
      finally={odbcSetAutoCommit(get("p_idaConnection",envir=idaRGlobal), autoCommit = TRUE);idaDropView(v)}
  );
}

idaQuery <- function (..., as.is = TRUE, na.strings = "NA")  {
  
  idaCheckConnection();
  
  query <- paste(..., sep = "", collapse = "")
  
  if(exists("p_debug",envir=idaRGlobal)&&(get("p_debug",envir=idaRGlobal))) { 
    print(query)
  }
  
  result <- sqlQuery(get("p_idaConnection", envir=idaRGlobal), query, believeNRows = FALSE, 
      stringsAsFactors = FALSE, as.is = as.is, na.strings = na.strings)
  
  #we got an error message
  if (is.character(result) && length(result) > 0) {
    stop(paste(result, collapse = "\n"))
  }
  
  return(result)
}

idaScalarQuery <- function (..., as.is=TRUE) {
  return(idaQuery(..., as.is=as.is)[[1]][1])
}

idaGetValidTableName <- function (prefix = "DATA_FRAME_") {
  idaCheckConnection();
  prefix <-toupper(prefix); 
  while (TRUE) {
    name = paste(prefix, floor(runif(1,0,100000)), sep="")
    if (!idaExistTable(name))
      return(name)
  }
}

callSP <- function (spname, ...) {
  idaCheckConnection()
  
  args = list(...)
  views <- c()
  tmp <- c()
  for (name in names(args)) {
    value <- args[[name]]
    if (is.null(value))   next
    if (is.ida.data.frame(value)) {
      view <- idaCreateView(value)
      tmp[length(tmp) + 1] = paste(name, "=", view, sep = "")
      views <- append(views, view)
    }
    else if (is.character(value) && length(grep(" ", value))) 
      tmp[length(tmp) + 1] = paste(name, "=\"", value, "\"", sep = "")
    else if (length(value) > 1) 
      tmp[length(tmp) + 1] = paste(name, "=\"", paste(value, collapse = " "), "\"", sep = "")
    else 
      tmp[length(tmp) + 1] = paste(name, "=", value, sep = "")
  }
  res <- try(idaQuery("CALL ", spname, "('", paste(tmp, collapse = ","), "')"), silent=T)
  for (view in views) idaDropView(view)
  if(inherits(res, "try-error")) {
    fullError <- idaScalarQuery("values idax.last_message")
    if(nchar(fullError)==0)
      fullError <- res;
    
    stop(fullError)
  }
  return(invisible(res))
}

################ Table def ############################
idaTableDef <- function(bdf,  collapse=TRUE) {
  if (!is.ida.data.frame(bdf))
    stop("this method can be applied only on ida.data.frame")
  
  viewName <- NULL;	
  if(length(idadf.defined.columns(bdf)>0)) {
    viewName <- idaCreateView(bdf);
    bdf <- ida.data.frame(viewName);
  }	
  
  attrs <- idaQuery("SELECT COLNAME, TYPENAME FROM SYSCAT.COLUMNS WHERE ",prepareCatalogLookup(bdf@table), " ORDER BY COLNO")
  
  if(!is.null(viewName)) {
    idaDropView(viewName);
  }
  
  idx <- tolower(attrs[,1]) %in% tolower(bdf@cols)
  
  if (as.logical(collapse))
    return(paste(attrs[idx,1], attrs[idx,2], collapse=","))
  else {
    res <- data.frame(name=attrs[idx,1], type=attrs[idx,2])
    res$valType <- ifelse(res$type %in% c('VARCHAR','CHARACTER','VARGRAPHIC','GRAPHIC','CLOB'),'CATEGORICAL',ifelse(res$type %in% c('SMALLINT', 'INTEGER','BIGINT','REAL','DOUBLE','FLOAT','DECIMAL','NUMERIC'),'NUMERIC','NONE'))
    return(res)
  }
}

idaMaterialize <- function(idf,tableName) {
  
  v <- idaCreateView(idf)
  
  
  tryCatch({
        
        odbcSetAutoCommit(get("p_idaConnection",envir=idaRGlobal), autoCommit = FALSE)
        idaQuery("CREATE TABLE ", tableName, " LIKE ",v," ORGANIZE BY ROW NOT LOGGED INITIALLY");
        idaQuery("INSERT INTO ", tableName, " SELECT * FROM ",v);
        
        odbcEndTran(get("p_idaConnection",envir=idaRGlobal), commit = TRUE);
        
      },
      error=function(e){print(e);odbcEndTran(get("p_idaConnection",envir=idaRGlobal), commit = FALSE)},
      finally={odbcSetAutoCommit(get("p_idaConnection",envir=idaRGlobal), autoCommit = TRUE);idaDropView(v)}
  );
}

idaAppend <- function(df, table) {
  sqlSave(get("p_idaConnection",envir=idaRGlobal), dat=df, tablename = table, rownames=F,append=T);
}

##############################  Utilities ############################

removeQuotes <- function(x) {
  if (!is.character(x) || length(x)>1)
    stop("This function can only be applied on a character value")
  if (nchar(x) > 2 && substr(x, 1, 1) == '"' && substr(x, nchar(x), nchar(x)) == '"') {
    x <- substr(x, 2, nchar(x)-1)
    return(x)
  }
  else
  	return (x)
}

colName <- function(colDef) {
  if (!inherits(colDef, 'ida.col.def'))
    stop("This function can only be applied on ida.col.def object")
  
  term = removeQuotes(colDef@term)
  tab <- colDef@table
  if (term %in% tab@cols)
    return(term)
    
  defCols <- tab@colDefs
  if (is.null(defCols) || length(defCols)==0)
    return(NULL)
  for (colName in names(defCols)) {
    if (term == defCols[[colName]])
      return(colName)
  }
}
