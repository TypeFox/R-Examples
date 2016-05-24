# CREATE and UPDATE the TAXONOMY TABLE
# package: megaptera
# author: Christoph Heibl
# last update: 2014-10-30
# TO DO: taxon

dbUpdateTaxonomy <- function(conn, taxonomy){
  
  if ( !inherits(conn, "PostgreSQLConnection") )
    stop("object 'conn' is not a valid PostgreSQL connection")
  
  ## write data to taxonomy
  if ( !missing(taxonomy) ){
    ## enforce synonyms column
    if ( !"synonyms" %in% names(taxonomy) ){
      taxonomy <- data.frame(taxonomy, synonyms = "-")
    }
    
    taxonomy <- sqlTaxonomyHeader(taxonomy)
    taxonomy$spec <- gsub(" ", "_", taxonomy$spec)
    
    if ( !dbExistsTable(conn, "taxonomy") ){
      ## create new taxonomy table
      ## -------------------------
      dbWriteTable(conn, "taxonomy", taxonomy, row.names = FALSE)
      sql <- "ALTER TABLE taxonomy ADD PRIMARY KEY (spec)"
      dbSendQuery(conn, sql)
      
    } else {
      ## update existing taxonomy table
      ## ------------------------------
      present <- dbReadTable(conn, "taxonomy")
      
      ## drop ranks higher than already existing in database
      ## ---------------------------------------------------
      id <- which(!is.na(match(colnames(taxonomy), colnames(present))))[1]
      taxonomy <- taxonomy[, id:ncol(taxonomy)]
      
      ## add columns missing in taxonomy
      ## -------------------------------
      id <- match(colnames(present), colnames(taxonomy))
      taxonomy <- t((t(taxonomy)[id, ]))
      id <- which(is.na(id))
      taxonomy[, id] <- "-"
      colnames(taxonomy)[id] <- colnames(present)[id]
      
      taxonomy <- data.frame(taxonomy)
      dbWriteTable(conn, "taxonomy", 
                   taxonomy[!taxonomy$spec %in% present$spec, ], 
                   row.names = FALSE, append = TRUE)
  } ## end of: write data to taxonony
  
  ## check consistency of taxonomy
  
  ## 1. mysterious wrong-genus-error
  ## -------------------------------
  tax <- dbReadTable(conn, "taxonomy")
  test <- tax$gen != strip.spec(tax$spec)
  if ( any(test) ){
    test <- paste("UPDATE taxonomy",
                  "SET", sql.wrap(strip.spec(tax$spec[test]), term = "gen", BOOL = NULL),
                  "WHERE", sql.wrap(tax$spec[test], term = "spec", BOOL = NULL))
    lapply(test, dbSendQuery, conn = conn)
   
  }
  
  

    
    ## update those tuples that have changed
    ## -------------------------------------
#     taxonomy <- taxonomy[taxonomy$spec %in% present$spec, ]
#     spec <- present$spec
#     syn <- which(names(present) %in% "synonyms")
#     for ( i in spec[11:100]){
#       id <- taxonomy[taxonomy$spec == i, -syn] == 
#         present[present$spec == i, -syn]
#       if ( !all(id) ){
#         print(i)
#         print("huu")
#       }
#     }
    
    
    ### ab hier beginnt der alte Code
    ### -----------------------------------------------
    
#     gene <- sql.conform(sproj$gene)
#     
#     ## 1. Detect and correct synonymous species names
#     ## ----------------------------------------------
#     spec <- dbGetQuery(conn, "SELECT spec FROM taxonomy")$spec
#     x <- lapply(spec, grep, x = sproj$tax$synonyms)
#     x <- data.frame(replace = spec[sapply(x, length) > 0],
#                     by = sproj$tax$spec[unlist(x[sapply(x, length) > 0])])
#     x <- x[!x$replace %in% x$by, ]
#     if ( nrow(x) > 0 ){
#       message(nrow(x), " synonymized species detected and corrected")
#       SQL <- paste("UPDATE taxonomy SET", sql.wrap(x$by, BOOL = NULL), 
#                    "WHERE", sql.wrap(x$replace, BOOL = NULL))
#       ## synonyms whose accepted name is already in database
#       id <- lapply(x$by, grep, x = spec)
#       id <- which(id > 0)
#       if ( length(id) > 0 ){
#         SQL[id] <- paste("DELETE FROM taxonomy WHERE", sql.wrap(x$replace[id], BOOL = NULL))
#       }
#       lapply(SQL, dbSendQuery, conn = conn)
#       ## correct also accessions table(s)
#       tabs <- dbListTables(conn)
#       tabs <- tabs[-grep("taxonomy|benchmark", tabs)]
#       SQL <- paste("SET", sql.wrap(x$by, term = "taxon", BOOL = NULL), 
#                    "WHERE", sql.wrap(x$replace, term = "taxon", BOOL = NULL))
#       SQL <- lapply(paste("UPDATE", tabs), function(x, y) paste(x, y), y = SQL)
#       lapply(unlist(SQL), dbSendQuery, conn = conn)
#     }
#     ## detect and correct 'lonesome' species in accession tables
#     ## ---------------------------------------------------------
#     spec <- dbGetQuery(conn, "SELECT spec FROM taxonomy")$spec
#     ssp <- vector()
#     for ( i in tabs ){
#       ti <- paste("SELECT DISTINCT taxon FROM", i)
#       si <- dbGetQuery(conn, ti)$taxon
#       #     print(ti); print(si[!si %in% spec])
#       ssp <- c(ssp, si[!si %in% spec])
#     }
#     ssp <- unique(ssp)
#     if ( length(ssp) > 0 ) warning(length(ssp), "species from accession tables not in taxonomy table")
#     
#     
#     ## 2. Correct higher rank taxonomy from database and from sproj
#     ## ------------------------------------------------------------
#     id <- sproj$tax$spec %in% tax$spec
#     message(length(which(!id)), " species not yet contained in database")
#     TAX <- sproj$tax[id, ]
#     TAX[sproj$synonyms] <- NULL
#     if ( !all(names(tax) == names(TAX)) )
#       warning("taxonomy tables contain different higher ranks")
#     TAX <- TAX[match(tax$spec, TAX$spec), ]
#     
#     ## update taxonomy in database if classification in query has changed
#     ## ------------------------------------------------------------------
#     ID <- which(!apply(TAX == tax, 1, all)) ## rownumber of species with conflicting higher ranks
#     if ( length(ID) > 0 ){
#       for ( i in ID ){
#         id <- tax[i, ] == TAX[i, ]
#         if ( !all(id) ){
#           cn <- colnames(id)[!id]
#           SQL <- as.vector(apply(TAX[i, cn, drop = FALSE], 1, as.character))
#           SQL <- paste(cn, "='", SQL, "'", sep = "")
#           SQL <- paste(SQL, collapse = ", ")
#           SQL <- paste("UPDATE taxonomy SET", SQL, 
#                        "WHERE", sql.wrap(tax$spec[i]))
#           cat("\n\t", SQL)
#           dbSendQuery(conn, SQL)
#         }
#       } # end of FOR loop over i
#     } else {
#       message("no taxonomic conflict found comparing ranks higher than species")
#     }
#     
#     ## names of taxa and their frequency
#     ## ---------------------------------
#     
#     tax <- paste("SELECT taxon, spec_ncbi FROM", gene)
#     tax <- dbGetQuery(conn, tax)
#     tax <- data.frame(table(tax$taxon, dnn = "taxon"))
#     
#     ## update <*_gb> columns in taxonomy
#     ## ---------------------------------
#     SQL <- paste("UPDATE taxonomy SET ", gene,  
#                  "_gb=", tax$Freq, " WHERE spec='", tax$taxon, "'",
#                  sep = "")
#     lapply(SQL, dbSendQuery, conn = conn)
#     
#     if ( !missing(set.null) ){
#       SQL <- paste("UPDATE taxonomy SET ", gene,  
#                    "_gb=NULL WHERE spec~'", set.null, "'",
#                    sep = "")
#       dbSendQuery(conn, SQL)
#     }
  }
}