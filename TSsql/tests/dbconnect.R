#  The real tests are in the database specific packages

require("TSdbi")

if(require("RMySQL") ) {
  # Before starting R you need to set user/passwd/host in ~/.my.cnf
  m <- dbDriver("MySQL")
  con <- try(dbConnect(m, dbname="test")) # pass user/passwd/host in ~/.my.cnf
  if(! inherits(con, "try-error")) {
     dbListTables(con) 
     cat("**************        disconnecting RMySQL test\n")
     dbDisconnect(con)
     } else  warning("RMySQL server error. Skipping tests.")
  dbUnloadDriver(m)
  } 

if(require("RSQLite") ) {

  m <- dbDriver("SQLite")
  f <- tempfile()
  con <- dbConnect(m, dbname=f) # no user/passwd/host
  dbListTables(con) 

  cat("**************        disconnecting SQLite test\n")
  dbDisconnect(con)
  unlink(f)
  #  dbUnloadDriver(m) complains about open connections.
  } 
