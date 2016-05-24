`RODM_list_dbms_models` <- function(database) {
#
# List Oracle Data Mining models present in the user's schema in the database.
#
  return(sqlQuery(database, "select * from user_mining_models"))
} # end of RODM_list_dbms_models
