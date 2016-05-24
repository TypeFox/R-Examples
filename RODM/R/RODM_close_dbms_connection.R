`RODM_close_dbms_connection` <- function(
#
# Close an ODBC connection to the Oracle Database
#
  database
)
{
  return(odbcClose(database))
} # end of RODM_close_dbms_connection
