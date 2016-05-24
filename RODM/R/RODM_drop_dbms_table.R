`RODM_drop_dbms_table` <-  function(
#
# Drop an Oracle Database table
#
database,                             # Name of the database ODBC channel
data_table_name)                      # Name of the table in the database
{
  res <- sqlQuery(database, 
                  paste("drop table ", data_table_name, " purge", 
                        sep="", collapse=""))
  if (identical(res, character(0))) invisible(TRUE) else res
} # end of RODM_drop_dbms_table
