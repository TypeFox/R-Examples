`RODM_apply_model` <- function( 
#
# Applies an Oracle Data Mining model
#
   database,                    # Database ODBC channel identifier
   data_table_name,             # Database table/view containing the training dataset
   model_name,                  # ODM model name
   supplemental_cols,           # Columns to retain in the output
   sql.log.file = NULL)         # File where to append the log of all the SQL calls made by this function
{
   if (!is.null(sql.log.file)) write(paste("--- SQL calls by ODM function: RODM_apply_model ", 
              date(), "---"), file = sql.log.file, append = TRUE, ncolumns = 1000)

   # Determine the mining function   
   query.string <- paste(
           "SELECT MINING_FUNCTION FROM USER_MINING_MODELS ",
           "WHERE MODEL_NAME = '", toupper(model_name), "'", sep="");
   func <- sqlQuery(database, query = query.string)

   slistcols <- ""
   for (i in 1:length(supplemental_cols)) {
     slistcols <- paste(slistcols,'"',toupper(supplemental_cols[i]),'", ',sep="")
   }

   if (func == 'CLASSIFICATION') {
     setop <- "PREDICTION_SET"
     predcol <- "PREDICTION"
     valcol <- "PROBABILITY"
   } else if (func == 'CLUSTERING') {
     setop <- "CLUSTER_SET"
     predcol <- "CLUSTER_ID"
     valcol <- "PROBABILITY"
   }

   # Apply ODM model
   if (func == 'REGRESSION') {
     query.string <- paste(
           "SELECT ", slistcols,
           "   prediction(", model_name, " using *) prediction ",
           "FROM ", data_table_name, " T", sep="");
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     apply.results <- sqlQuery(database, query = query.string)
   } else {
     query.string <- paste(
           "SELECT ", predcol, " FROM TABLE(", setop, "(",
           model_name, " USING))", sep="");
     classlist <- sqlQuery(database, query = query.string)
     pivotlist <- ""
     pivotlistcols <- ""
     for (i in 1:length(classlist[,1])) {
       if (i > 1) {pivotlist <- paste(pivotlist,", ",sep="")}
       pivotlist <- paste(pivotlist,"'",classlist[i,1],"'",sep="")
       pivotlistcols <- paste(pivotlistcols,'"',"'",classlist[i,1],"'",'", ',sep="")
     }
     query.string <- paste(
           "SELECT ", pivotlistcols, slistcols, predcol, " FROM (",
           " SELECT rnum, ", slistcols, 
                    "S.", predcol, " rodmp1, S.", valcol, " rodmp2, T.",
                    predcol, " FROM (",
           "  SELECT ",
                predcol, "(", model_name, " using *) ", predcol, ", ",
                 slistcols,
           "     row_number() over (order by 1) rnum, ",
                 setop, "(", model_name, " using *) resset ",
           "  FROM ", data_table_name, ") T, table(T.resset) S) ",
           " PIVOT (max(rodmp2) for rodmp1 in (", pivotlist, ")) ", sep="");
     if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
     apply.results <- sqlQuery(database, query = query.string)
  }

  model.apply.list <- list("model.apply.results" = apply.results)

  return(model.apply.list)
} # End of RODM_apply_model
