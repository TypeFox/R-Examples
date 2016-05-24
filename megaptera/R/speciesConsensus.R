# construct species consensus sequences
# PACKAGE: megaptera
# CALLED BY: stepE
# AUTHOR: Christoph Heibl
# LAST UPDATE: 2014-10-30

speciesConsensus <- function(megProj, spec){
  
  ## PARAMETERS
  ## -----------
  gene <- megProj@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  spec.tab <- paste("spec", gsub("^_", "", gene), sep = "_")
  align.exe <- megProj@align.exe
  max.bp <- megProj@params@max.bp
  max.dist <- megProj@params@max.dist
  logfile <- paste(gene, "stepE.log", sep = "-")
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(megProj@db)
  
  ## read species alignments
  ## -----------------------
  obj <- dbReadDNA(spec, conn = conn, tab.name = acc.tab, 
                   max.bp = max.bp, max.dist = max.dist, 
                   enforce.binomial = TRUE)
  
  ## check for aberrant sequences
  ## SHOULD BE REMOVED TO STEP D?!
  ## ------------------------------
  #   d <- lapply(obj, myDist)
  #   md <- sapply(d, max, na.rm = TRUE)
  #   md <- names(md[md > .1])
  #   if ( length(md) > 0 ){
  #     md <- sql.wrap(md, term = "taxon")
  #     md <- paste("(", md, ")", sep = "")
  #     md <- paste("SELECT * FROM", acc.tab, "WHERE npos <=", max.bp, "AND distreference <=", max.dist, 
  #                 "AND", md, "ORDER BY taxon")
  #     md <- dbGetQuery(conn = conn, md)
  #     s <- split(md["distreference"], f = md$taxon)
  #     s <- sapply(s, min) * 1.25
  #     gis <- vector()
  #     for ( i in seq_along(s) ){
  #       id <- md[md$taxon == names(s)[i] & md$distreference < s[i], "gi"]
  #       gis <- c(gis, id)
  #       id <- paste(id, collapse = "|")
  #       txn <- names(s)[i]
  #       obj[[txn]] <- obj[[txn]][grep(id, rownames(obj[[txn]])), ]
  #       obj[[txn]] <- deleteEmptyCells(obj[[txn]], quiet = TRUE)
  #     }
  #     gis <- setdiff(md$gi, gis)
  #   } else {
  #     gis <- NULL
  #   }
  
  ## update status field
  ## -------------------
#   SQL <- paste("UPDATE", acc.tab, 
#                "SET status=status || '-selected'", 
#                "WHERE status !~ 'selected' AND npos <=", max.bp,
#                "AND distreference <=", max.dist)
#   dbSendQuery(conn, SQL)
  #   if ( !is.null(gis) ){
  #     SQL <- c(SQL,
  #              paste("UPDATE " , acc.tab, " SET status='too distant (from conspecific)' WHERE gi='", gis, "'", sep = ""))
  #   }
  
  
  ## update locus table
  ## ---------------------
  nb.acc <- nrow(obj)
  SQL <- paste("UPDATE locus", 
               "SET", sql.wrap(nb.acc, term = paste(gene, "sel", sep = "_")),
               "WHERE", sql.wrap(spec))
  dbSendQuery(conn, SQL)
#   tax <- dbReadTable(conn, "locus")$spec
#  # tax <- setdiff(tax, union(spec, already))
#   tax <- setdiff(tax, spec)
#   SQL <- paste("UPDATE locus",
#                "SET", sql.wrap(0, term = paste(gene, "sel", sep = "_")),
#                "WHERE", sql.wrap(tax, BOOL = NULL))
#   lapply(SQL, dbSendQuery, conn = conn)
  
  ## species consensus sequences
  ## ---------------------------
  obj <- specCons(obj, log = logfile)
  obj <- list(obj)
  names(obj) <- spec
  class(obj) <- "DNAbin"
  
  write.dna.spectable(conn, spec.tab, obj)

  dbDisconnect(conn)
}