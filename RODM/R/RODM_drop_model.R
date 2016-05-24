`RODM_drop_model` <- function(
#
# Drop an Oracle Data Mining model
#
  database,
  model_name,
  sql.log.file = NULL) 
{
  query.string <- paste("BEGIN dbms_data_mining.drop_model('", model_name, 
                        "'); END;", sep="", collapse="")
  if (!is.null(sql.log.file)) write(query.string, file = sql.log.file, append = TRUE, ncolumns = 1000)
  res <- sqlQuery(database, query = query.string)
  if (identical(res, character(0))) invisible(TRUE) else stop(res)
} # end of RODM_drop_model
