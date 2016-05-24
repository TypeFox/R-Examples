# INSERT MISSING SPECIES into TAXONOMY TABLE
# package: megaptera
# called by: stepC
# author: Christoph Heibl
# last update: 2014-09-23

dbAutoFillTaxonomy <- function(x){
  
  ## DEFINITIONS
  ## -----------
  gene <- sql.conform(x$gene)
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  
  ## open database connection
  ## ------------------------
  conn <- dbConnect(PostgreSQL(), host = x$db[["host"]], 
                    port = x$db[["port"]], dbname = x$db[["dbname"]], 
                    user = x$db[["user"]], password = x$db[["password"]])
  
  SQL <- paste("SELECT taxon from", acc.tab, 
               "WHERE status != 'excluded' OR status IS NULL",
               "EXCEPT SELECT spec from taxonomy")
  spec <- dbGetQuery(conn, SQL)$taxon
  
  
  autoFill <- function(z, conn){
    zz <- paste("SELECT * FROM taxonomy WHERE", 
                sql.wrap(strip.spec(z), term = "gen"))
    zz <- dbGetQuery(conn, zz)
    if ( nrow(zz) == 0 ) {
      warning("genus '", strip.spec(z), "' not in taxonomy")
      SQL <- paste("UPDATE", acc.tab, 
                   "SET", sql.wrap("excluded", term = "status"),
                   "WHERE", sql.wrap(z, term = "taxon", operator = "~", BOOL = NULL))
      dbSendQuery(conn, SQL)
    } else {
      zz <- zz[, -which(colnames(zz) == "spec")]
      zz <- unique(zz)[1, ]
      zz <- data.frame(zz, spec = z, stringsAsFactors = FALSE)
      SQL <- paste("INSERT INTO taxonomy",
                   paste("(", paste(names(zz), collapse = ", "), ")", sep = ""),
                   "VALUES",
                   paste("(", paste(paste("'", zz, "'", sep = ""), collapse = ", "), ")", sep = ""))
      dbSendQuery(conn, SQL)
    }
  }
  lapply(spec, autoFill, conn = conn)
  dbDisconnect(conn)
}