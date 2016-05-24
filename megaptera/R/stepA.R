# SEARCH AND DOWNLOAD TAXONOMY in serial/parallel execution
# package: megaptera
# author: Christoph Heibl
# last change 2014-11-03

stepA <- function(x){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  
  ## PARAMETERS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  
  ## iniate logfile
  ## --------------
  logfile <- paste(acc.tab, "stepA.log", sep = "-")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""), 
       "\nSTEP A: searching and downloading taxonomy from GenBank\n",
       file = logfile)
  
  conn <- dbconnect(x@db)
  
  if ( !dbExistsTable(conn, "taxonomy") | x@update ){
    
    ## download NCBI taxonomy for ingroup
    ingroup <- ncbiTaxonomy(x@taxon@ingroup, 
                            species.list = x@taxon@species.list,
                            kingdom = x@taxon@kingdom)
    ingroup <- fixTaxonomy(ingroup, auto = TRUE)
    dbUpdateTaxonomy(conn, ingroup)
    
    ## download NCBI taxonomy for ingroup
    outgroup <- ncbiTaxonomy(x@taxon@outgroup, 
                             species.list = x@taxon@species.list,
                             kingdom = x@taxon@kingdom)
    outgroup <- fixTaxonomy(outgroup, auto = TRUE)
    dbUpdateTaxonomy(conn, outgroup)
    
  } else {
    slog("\ntaxonomy already downloadad", file = logfile)
  }


  dbDisconnect(conn)
  
  slog("\n\nSTEP A finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
}