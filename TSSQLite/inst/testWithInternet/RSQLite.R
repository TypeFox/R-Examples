# sourcing HistQuote.TSsql requires the Internet

cat("************** RSQLite  Examples ******************************\n")

# no user/passwd/host
dbfile <- tempfile()
setup <- RSQLite::dbConnect(RSQLite::SQLite(), dbname=dbfile) 

RSQLite::dbListTables(setup) 

TSsql::removeTSdbTables(setup, yesIknowWhatIamDoing=TRUE)
TSsql::createTSdbTables(setup, index=FALSE)

DBI::dbListTables(setup) 
DBI::dbDisconnect(setup)

require("TSSQLite")

con <- try(TSconnect("SQLite", dbname=dbfile) )
if(inherits(con, "try-error")) stop("CreateTables did not work.")

source(system.file("TSsql/Populate.TSsql", package = "TSsql"))
source(system.file("TSsql/TSdbi.TSsql", package = "TSsql"))
source(system.file("TSsql/dbGetQuery.TSsql", package = "TSsql"))
source(system.file("TSsql/HistQuote.TSsql", package = "TSsql"))

cat("**************        remove test tables\n")
TSsql::removeTSdbTables(con, yesIknowWhatIamDoing=TRUE)

cat("**************        disconnecting test\n")
dbDisconnect(con)
unlink(dbfile)

