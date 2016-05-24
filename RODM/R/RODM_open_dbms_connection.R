`RODM_open_dbms_connection` <- function(
#
# Open an ODBC connection to the Oracle Database
#
dsn,
uid = "",
pwd = "")
{
  library(RODBC)

  channel <- odbcConnect(dsn = dsn, uid = uid, pwd = pwd, case = "toupper", rows_at_time=1)
  driver <- odbcGetInfo(channel)[[1]]
  setSqlTypeInfo(driver, list(double="double precision", integer="integer", character="varchar(255)", logical="varchar(255)"))

  # Check that we are running against Oracle version 11 or later
  query.string <- "select to_number(substr(version,1,instr(version,'.')-1)) 
      from product_component_version where product like 'Oracle Database%'"
  dbvers <- sqlQuery(channel, query = query.string)
  if (dbvers < 11) {
    RODM_close_dbms_connection(channel)
    stop("RODM requires Oracle Database version 11 or higher")
  }

  # Check that ODM is installed
  query.string <- "select value from v$option where parameter = 'Data Mining'"
  odminst <- sqlQuery(channel, query = query.string)
  if (odminst != "TRUE") {
    RODM_close_dbms_connection(channel)
    stop("Oracle Data Mining option is not installed on the Oracle Database")
  }

  # Check for database privileges necessary to perform mining
  queryprefix.string <-
    "select count(*) from 
     (select privilege from user_sys_privs
      union all
      select privilege from role_sys_privs)
     where privilege in ";
  query.string <- paste(queryprefix.string, "('CREATE MINING MODEL', 'CREATE ANY MINING MODEL')")
  modelpriv <- sqlQuery(channel, query = query.string)
  query.string <- paste(queryprefix.string, "('CREATE TABLE', 'CREATE ANY TABLE')")
  tablepriv <- sqlQuery(channel, query = query.string)
  query.string <- paste(queryprefix.string, "('CREATE VIEW', 'CREATE ANY VIEW')")
  viewpriv <- sqlQuery(channel, query = query.string)
  if ((modelpriv == 0) | (tablepriv == 0) | (viewpriv == 0)) {
    RODM_close_dbms_connection(channel)
    if (modelpriv == 0) stop("Database user lacks CREATE MINING MODEL privilege, which is necessary for mining")
    if (tablepriv == 0) stop("Database user lacks CREATE TABLE privilege, which is necessary for mining")
    if (viewpriv == 0) stop("Database user lacks CREATE VIEW privilege, which is necessary for mining")
  }

  # Create a settings table (global temporary) if one does not exist yet
  if (sqlQuery(channel, 
        "select count(*) from user_tables where table_name = 'RODM_SETTINGS_TABLE'") 
      == 0) {
    query.string <- "CREATE GLOBAL TEMPORARY TABLE RODM_SETTINGS_TABLE (SETTING_NAME VARCHAR2(30), SETTING_VALUE VARCHAR2(4000)) ON COMMIT PRESERVE ROWS"
    sqlQuery(channel, query = query.string)
  }

  # Create tables for bias settings if they do not exist yet
  if (sqlQuery(channel, 
        "select count(*) from user_tables where table_name = 'RODM_NUM_PRIORS'") 
      == 0) {
    query.string <- "CREATE GLOBAL TEMPORARY TABLE RODM_NUM_PRIORS (TARGET_VALUE NUMBER, PRIOR_PROBABILITY NUMBER) ON COMMIT PRESERVE ROWS"
    sqlQuery(channel, query = query.string)
  }
  if (sqlQuery(channel, 
        "select count(*) from user_tables where table_name = 'RODM_CAT_PRIORS'") 
      == 0) {
    query.string <- "CREATE GLOBAL TEMPORARY TABLE RODM_CAT_PRIORS (TARGET_VALUE VARCHAR2(4000), PRIOR_PROBABILITY NUMBER) ON COMMIT PRESERVE ROWS"
    sqlQuery(channel, query = query.string)
  }
  if (sqlQuery(channel, 
        "select count(*) from user_tables where table_name = 'RODM_NUM_COSTS'") 
      == 0) {
    query.string <- "CREATE GLOBAL TEMPORARY TABLE RODM_NUM_COSTS (ACTUAL_TARGET_VALUE NUMBER, PREDICTED_TARGET_VALUE NUMBER, COST NUMBER) ON COMMIT PRESERVE ROWS"
    sqlQuery(channel, query = query.string)
  }
  if (sqlQuery(channel, 
        "select count(*) from user_tables where table_name = 'RODM_CAT_COSTS'") 
      == 0) {
    query.string <- "CREATE GLOBAL TEMPORARY TABLE RODM_CAT_COSTS (ACTUAL_TARGET_VALUE VARCHAR2(4000), PREDICTED_TARGET_VALUE VARCHAR2(4000), COST NUMBER) ON COMMIT PRESERVE ROWS"
    sqlQuery(channel, query = query.string)
  }

  return(channel)
} # end of RODM_open_dbms_connection
