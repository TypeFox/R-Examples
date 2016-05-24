# ALIGN CONSPECIFIC SEQUENCES
# package: megaptera
# called by: USER
# author: Christoph Heibl
# last update: 2014-11-01

stepC <- function(x){
  
  start <- Sys.time()
  pipeline <- TRUE
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x@db)
  
  ## DEFINITIONS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  tax <- dbReadTaxonomy(conn)
  align.exe <- x@align.exe
  max.bp <- x@params@max.bp
  
  ## iniate logfile
  ## --------------
  logfile <- paste(gene, "stepC.log", sep = "-")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""), 
       "\n\nSTEP C: alignment of conspecific sequences\n",
       file = logfile)
  
  ## check
  if ( !dbExistsTable(conn, acc.tab) ) {
    dbDisconnect(conn)
    stop("stepB has not been called yet")
  } 
  
  ## get table of species available in database
  ## ------------------------------------------
  SQL <- paste("SELECT taxon AS spec, count(taxon)",
               "FROM", acc.tab, 
               "WHERE npos <=", max.bp, "AND status ~ 'raw|aligned'",
               "GROUP BY taxon") 
  tax <- dbGetQuery(conn, SQL)
  if ( nrow(tax) == 0 ) {
    dbDisconnect(conn)
    slog("no sequences - try to rerun stepB", file = logfile)
#     x$evaluate <- FALSE
    slog("\n\nSTEP C finished", file = logfile)
    td <- Sys.time() - start
    slog(" after", round(td, 2), attr(td, "units"), file = logfile)
    return(x)
  }
  
  ## mark single-sequence species in <status>
  ## ----------------------------------------
  SQL <- paste("UPDATE", acc.tab, 
               "SET status = 'single'",
               "WHERE", sql.wrap(tax$spec[tax$count == 1], term = "taxon", BOOL = NULL))
  lapply(SQL, dbSendQuery, conn = conn)
  
  ## mark sequences too long to align in <status>
  ## --------------------------------------------
  SQL <- paste("UPDATE", acc.tab, 
               "SET status = 'excluded (too long)'",
               "WHERE npos >", max.bp)
  dbSendQuery(conn, SQL)
  
  ## select species for which we have more than 1 acc.
  ## -------------------------------------------------
  spec <- tax$spec[tax$count > 1] 
  slog("\n", length(spec), "species need to be aligned", 
       file = logfile)
  
  ## aligning -- either sequential or parallel
  ## -----------------------------------------
  if ( length(spec) > 0 ) {
    cpus <- 12
    if ( length(spec) < cpus | x@parallel == FALSE ){
      lapply(spec, funX, megProj = x, acc.tab = acc.tab, max.bp = max.bp, 
               align.exe = align.exe, logfile = logfile)
    } else {
      sfInit(parallel = TRUE, cpus = cpus, type = "SOCK")
      sfLibrary("megaptera", character.only = TRUE)
      megProj <- x
      sfExport("spec", "megProj", "acc.tab", 
               "max.bp", "align.exe", "logfile")
      sfLapply(x = spec, fun = funX, megProj = megProj, acc.tab = acc.tab, max.bp = max.bp, 
               align.exe = align.exe, logfile = logfile)
      sfStop()
    }
  } 
  
  dbDisconnect(conn)
  
  slog("\n\nSTEP C finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), 
       "\n", file = logfile)
  invisible(x)
}
