# Calculation of each accessions distance to reference sequence
# PACKAGE: megaptera
# AUTHOR: Christoph Heibl
# LAST UPDATE: 2014-11-03
# TO DO: 
# - sequences that do not overlap with reference sequence receive dist=1
#   although they might be 'good'

stepE <- function(x){
  
#   if ( !x$evaluate ) return(x)
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  align.exe <- x@align.exe
  max.bp <- x@params@max.bp
  max.dist <- x@params@max.dist
  parallel <- x@parallel
  
  ## iniate logfile
  ## --------------
  logfile <- paste(gene, "stepE.log", sep = "-")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""),
       "\nSTEP E: calculate genetic distances from reference", 
       file = logfile)
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x@db)
  
  ## set taxonomic rank for 
  ## reference sequence calculation
  ## ------------------------------
  rr <- x@taxon@reference.rank
  if ( rr == "auto" ){
    tax <- dbReadTaxonomy(conn)
    if ( x@taxon@species.list ){
      rr <- apply(tax, 2, function(x) length(unique(x)))
      rr <- head(rr, -1) ## delete synonyms column
      rr <- names(rr)[max(which(rr == 1))]
    } else {
      rr <- apply(tax, 2, grep, pattern = paste("^", "$", 
                                                sep = x@taxon@ingroup))
      rr <- names(rr)[sapply(rr, length) > 0]
    }
  }
  
  ## read reference sequences
  ## ------------------------
  reference <- dbReadReference(conn, gene)
  
  ## clear reference column if all entries are to be updated
  ## -------------------------------------------------------
  if ( x@update ) dbSendQuery(conn, paste("UPDATE", acc.tab, "SET distreference = NULL"))
  
  ## table of species names + reference sequence
  ## -------------------------------------------
  slog("\n.. reading species names: ", file = logfile)
  tax <- paste("SELECT DISTINCT spec,", rr, "AS ref",
               "FROM", acc.tab, "JOIN taxonomy ON spec = taxon ",
               "WHERE status !~ 'excluded|too' AND distreference IS NULL")
  tax <- dbGetQuery(conn, tax)
  slog(nrow(tax), "found ...", file = logfile)
  
  if ( nrow(tax) == 0 ) {
    slog("\n.. database is up to date -- nothing to do", file = logfile)
    dbDisconnect(conn)
#     x$evaluate <- FALSE
    slog("\n\nSTEP E finished", file = logfile)
    td <- Sys.time() - start
    slog(" after", round(td, 2), attr(td, "units"), file = logfile)
    return(x)
  }
  
  ## select the 'best' sequences
  ## ---------------------------
  slog("\n.. calculating genetic distance from 'reference' for ...", 
       file = logfile)
  
  ## distance from reference -- either sequential or parallel
  ## --------------------------------------------------------
  if ( nrow(tax) > 0 ) {
    cpus <- 12
    if ( nrow(tax) < cpus | !x@parallel ){
      apply(tax, 1, calcDistToRef, 
            megProj = x, max.bp = max.bp, reference = reference, 
            align.exe = align.exe, acc.tab = acc.tab, logfile = logfile)
    } else {
      slog("\n", file = logfile)
      sfInit(parallel = TRUE, cpus = cpus, type = "SOCK")
      sfLibrary("megaptera", character.only = TRUE)
      megProj <- x
      sfExport("reference", "megProj", "logfile", "acc.tab", "align.exe", "max.bp")
      sfApply(tax, 1, calcDistToRef, megProj = megProj, max.bp = max.bp,
              reference = reference, align.exe = align.exe, acc.tab = acc.tab, logfile = logfile)
      sfStop()
    }
  }
  
  dbDisconnect(conn)
  slog("\n\nSTEP E finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
  invisible(x)
}
