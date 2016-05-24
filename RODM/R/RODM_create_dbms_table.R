`RODM_create_dbms_table` <- function(  
#
# Moves data (frame) table to an Oracle Database table
#
database,                             # Name of the database ODBC channel
data_table_name)                      # Name of the table in the database
{
  if (sqlQuery(database, 
        paste("select count(*) from user_tables where table_name = '",
              toupper(data_table_name), "'", sep="", collapse="")) > 0) {
    RODM_drop_dbms_table(database, data_table_name)
  }
  return(sqlSave(database, dat=eval(parse(text=data_table_name)), 
    tablename = data_table_name, rownames = FALSE, safer = FALSE, test = FALSE))
}  # end of RODM_create_dbms_table
