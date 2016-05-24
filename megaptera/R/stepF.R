# construct species consensus sequences
# PACKAGE: megaptera
# AUTHOR: Christoph Heibl
# LAST UPDATE: 2014-10-30
# TO DO: line 82-94 update must be based on spec table

stepF <- function(x){
  
#   if ( !x$evaluate ) return(x)
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  
  ## PARAMETERS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  spec.tab <- paste("spec", gsub("^_", "", gene), sep = "_")
  align.exe <- x@align.exe
  max.bp <- x@params@max.bp
  max.dist <- x@params@max.dist
  
  ## iniate logfile
  ## --------------
  logfile <- paste(gene, "stepF.log", sep = "-")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""),
       "\nSTEP F: construct species consensus sequences\n", file = logfile)
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x@db)
  
  ## create table (if it does not exist)
  ## -----------------------------------
  if ( !dbExistsTable(conn, spec.tab) ) {
    slog("\ncreating TABLE", spec.tab, file = logfile)
    SQL <- paste(spec.tab, "_pk", sep = "")
    SQL <- paste("CREATE TABLE", spec.tab, 
                 "(spec text NOT NULL,",
                 "block text,",
                 "npos integer NOT NULL,", 
                 "dna text NOT NULL,",
                 "CONSTRAINT", SQL, "PRIMARY KEY ( spec ))")
    dbSendQuery(conn, SQL)
  }
  
  ## clear *_sel column in taxonomy table from previous calls to stepF
  ## -----------------------------------------------------------------
  #   SQL <- paste("UPDATE taxonomy SET ", gene, "_sel = NULL", sep = "")
  #   dbSendQuery(conn = conn, SQL)
  
  ## vector of species names in table <acc.tab>
  ## ------------------------------------------
  SQL <- paste("SELECT taxon AS spec, count(taxon)",
               "FROM", acc.tab, 
               "WHERE distreference <=", max.dist, 
                    "AND status !~ 'excluded|too'",
               "GROUP BY taxon")
  taxa <- dbGetQuery(conn, SQL)
  slog(paste("\n ", nrow(taxa), " species found in table '", acc.tab, "'", sep = ""), file = logfile)
  
  ## update <status> for accessions *NOT* selected
  ## ---------------------------------------------
  SQL <- c(paste("UPDATE", acc.tab, "SET status='too long' WHERE npos >", max.bp),
           paste("UPDATE", acc.tab, "SET status='too distant (from reference)' WHERE distreference >", max.dist))
  lapply(SQL, dbSendQuery, conn = conn)
  
  ## intersect with species names in query
  ## --------------------------------------
  #   taxa <- intersect(taxa, x$tax$spec)
  #   slog("\n", length(taxa), "species of which belong to current query", file = logfile)
  
  if ( nrow(taxa) == 0 ) {
    dbDisconnect(conn)
    stop("stepE has not been called yet")
  }
#   slog("\n  selecting: ", file = logfile)
#   slog("\n  - sequences shorter than", max.bp, file = logfile)
#   slog("\n  - sequences more similar to reference than", max.dist, file = logfile)
  
  ## check if sequences have been selected
#   fn <- paste(gene, "fas", sep = ".") # also needed for update == TRUE
#   if ( !update ){
#     if ( file.exists(fn) ){
#       seqs <- as.list(read.fas(fn))
#       already <- intersect(taxa, names(seqs))
#       taxa <-  setdiff(taxa, names(seqs))
#     } else {
#       already <- vector()
#     }
#   } else {
    already <- vector()
#   }
  
  if ( nrow(taxa) == 0 ){
    slog("\n\n  all suitable sequences already selected", file = logfile)
    dbDisconnect(conn)
  } else {
    
    #single.seq <- taxa$spec[taxa$count == 1]
    
    ## species consensus -- either sequential or parallel
    ## --------------------------------------------------
    if ( nrow(taxa) > 0 ) {
      cpus <- 12
      if ( nrow(taxa) < cpus | !x@parallel ){
        lapply(taxa$spec, speciesConsensus, megProj = x)
      } else {
        sfInit(parallel = TRUE, cpus = cpus, type = "SOCK")
        sfLibrary("megaptera", character.only = TRUE) 
        megProj <- x
        sfExport("megProj", "taxa")
        sfLapply(taxa$spec, speciesConsensus, megProj = megProj)
        sfStop()
      }
    }

  }
  dbDisconnect(conn)
  slog("\n\nSTEP F finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
  invisible(x)
}
