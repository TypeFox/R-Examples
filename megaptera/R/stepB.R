# SEARCH AND DOWNLOAD SEQUENCES in sereial/parallel execution
# package: megaptera
# author: Christoph Heibl
# last change 2014-11-06

stepB <- function(x){
  
  start <- Sys.time()
  
  ## CHECKS
  ## ------
  if ( !inherits(x, "megapteraProj") )
    stop("'x' is not of class 'megapteraProj'")
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x@db)
  
  ## PARAMETERS
  ## -----------
  update.seqs <- "all"
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  
  ## iniate logfile
  ## --------------
  logfile <- paste(acc.tab, "stepB.log", sep = "-")
  if ( file.exists(logfile) ) unlink(logfile)
  slog(paste("\nmegaptera", 
             packageDescription("megaptera")$Version),
       paste("\n", Sys.time(), sep = ""), 
       "\nSTEP B: searching and downloading sequences from GenBank\n",
       file = logfile)
  
  ## delete table (if it exists and is not to be updated)
  ## ----------------------------------------------------
  if ( dbExistsTable(conn, acc.tab) & update.seqs == "all" ) {
    dbRemoveTable(conn, acc.tab)
    slog("\n.. removing TABLE", acc.tab, "..", file = logfile)
  }
  
  ## create table (if it does not exist)
  ## -----------------------------------
  if ( !dbExistsTable(conn, acc.tab) ) {
    slog("\n.. creating TABLE", acc.tab, "..", file = logfile)
    SQL <- paste(acc.tab, "_pk", sep = "")
    SQL <- paste("CREATE TABLE", acc.tab, 
                 "(gi text NOT NULL,",
                 "taxon text NOT NULL,",
                 "spec_ncbi text NOT NULL,",
                 "status text,",
                 "genom text,",
                 "npos integer NOT NULL,", 
                 "distreference real,",
                 "dna text NOT NULL,",
                 "CONSTRAINT", SQL, "PRIMARY KEY ( gi ))")
    dbSendQuery(conn, SQL)
  }
  
  ## list of species or higher taxa to be be searched for
  ## ----------------------------------------------------
  slog("\n.. assembling taxon search list ..", file = logfile)
  
  if ( x@taxon@species.list ){
    st <- "spec"
    tax <- dbReadTaxonomy(conn)
    search.tax <- unique(tax[, st])
    search.tax <- gsub("_", " ", search.tax)
  } else {
    search.tax <- c(x@taxon@ingroup, 
                    x@taxon@outgroup)
  } 
  
  ## search and download sequences
  ## -----------------------------
  lapply(search.tax, downloadSequences, x = x)
  
  #   if ( dbExistsTable(conn, "locus") ){
  #     present.markers <- names(dbGetQuery(conn, paste("SELECT * FROM locus")))
  #     if ( update.seqs != "all" & cols[1] %in% present.markers ) {
  #       if ( update.seqs == "not.yet.searched" )
  #         not.search <- paste("SELECT spec FROM locus WHERE ", 
  #                             gene, "_gb IS NOT NULL", sep = "")
  #       if ( update.seqs == "not.yet.found" )
  #         not.search <- paste("SELECT spec FROM locus WHERE ", 
  #                             gene, "_gb > 0", sep = "")
  #       not.search <- dbGetQuery(conn, not.search)$spec
  #       search.tax <- search.tax[!names(search.tax) %in% not.search]
  #     } 
  #   }
  
  ## handle subspecies and varieties as species
  ## ------------------------------------------
  if ( !x@taxon@species.list ){
    taxon <- paste("SELECT taxon FROM", acc.tab)
    taxon <- unique(dbGetQuery(conn, taxon)$taxon)
    taxon <- data.frame(gb = taxon,  spec = strip.infraspec(taxon), 
                        stringsAsFactors = FALSE)
    
    id <- c(grep("^'", taxon$gb), ## "'Asterotremella_humicola'"
            which(taxon$gb != taxon$spec)) ## intraspecific ranks
    if ( length(id) > 0 ){
      taxon <- taxon[id, ] 
      if ( !x@taxon@hybrids )
        exclude <- grep("_x_|^x_", taxon$gb)
      exclude <- c(exclude, grep("^'", taxon$gb), grep("sp[.]|var[.]|cf[.]|aff[.]|hybrid$|environmental$|.[[:upper:]]", taxon$spec))
      
      if ( length(exclude) == 0 ){
        rename <- taxon
      } else {
        rename <- taxon[-exclude, ]
        exclude <- taxon[exclude, ]
        
        ## declare excluded taxa
        ## ---------------------
        SQL <- paste("UPDATE", acc.tab, 
                     "SET status = 'excluded (indet)'",
                     "WHERE", sql.wrap(exclude$gb, term = "taxon", operator = "~", BOOL = NULL))
        lapply(SQL, dbSendQuery, conn = conn)
      }
      ## rename infraspecific taxa
      ## -------------------------
      SQL <- paste("UPDATE", acc.tab, 
                   "SET", sql.wrap(rename$spec, term = "taxon", BOOL = NULL),
                   "WHERE", sql.wrap(rename$gb, term = "taxon", BOOL = NULL, operator = "~"))
      lapply(SQL, dbSendQuery, conn = conn)
    }
  } # end of IF (l.99)
  
  ## select sequences if there are > max.gi.per.spec
  ## -----------------------------------------------
  taxon <- paste("SELECT gi, taxon, npos", 
                 "FROM", acc.tab,
                 "WHERE npos <=", x@params@max.bp)
  taxon <- dbGetQuery(conn, taxon)
  freqs <- table(taxon$taxon)
  if ( any(id <- freqs > x@params@max.gi.per.spec) ){
    id <- names(id)[id]
    slog("\n..", length(id), "species have >", 
         x@params@max.gi.per.spec, "sequenes:", 
         file = logfile)
    for ( i in id ){
      acc <- taxon[taxon$taxon == i, ]
      acc <- acc[order(acc$npos, decreasing = TRUE), ]
      acc <- tail(acc$gi, -x@params@max.gi.per.spec)
      slog("\n -", i, "(", length(acc), "sequences excluded )",
           file = logfile)
      acc <- paste("UPDATE", acc.tab, 
                   "SET status='excluded (max.gi)'",
                   "WHERE", sql.wrap(acc, term = "gi", BOOL = "OR"))
      dbSendQuery(conn, acc)
      
    }
  }
  
  ## create and update relation <locus>
  ## ----------------------------------
  dbUpdateLocus(x)
  
  # summary
  # -------
  check <- paste("SELECT taxon FROM", acc.tab)
  check <- table(dbGetQuery(conn, check)$taxon)
  dbDisconnect(conn)
  total <- length(check)
  one <- length(check[check == 1])
  multiple <- length(check[check > 1])
  slog(paste("\n\nNumber of species in database:", total),
       paste("\nNumber of species with one accession: ", one, 
             " (", round(one * 100 / total, 1), "%)", sep = ""),
       paste("\nNumber of species with multiple accessions: ", multiple, 
             " (", round(multiple * 100 / total, 1), "%)", sep = ""),
       paste("\nTotal number of species:", one + multiple),
       paste("\nTotal number of accessions:", sum(check)), 
       file = logfile)
  
  
  slog("\n\nSTEP B finished", file = logfile)
  td <- Sys.time() - start
  slog(" after", round(td, 2), attr(td, "units"), file = logfile)
}
