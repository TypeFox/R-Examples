alignGenus <- function(genus, megProj){
  
  ## PARAMETERS
  ## -----------
  gene <- megProj@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  spec.tab <- paste("spec", gsub("^_", "", gene), sep = "_")
  logfile <- paste(gene, "stepF.log", sep = "-")
  align.exe <- megProj@align.exe
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(megProj@db)

  slog("\n-- ", genus, file = logfile)
  seqs <- dbReadDNA(conn, spec.tab, taxon = genus)
  seqs <- as.list(seqs)
  if ( length(seqs) > 1 ) {
    
    ## aligning with benchmark
#     b <- read.benchmark.db(conn, gene)
#     br <- unique(x$tax[x$tax$gen == i, x$benchmark.rank])
#     # benchmark is not available:
#     if ( is.null(b[[br]]) ){
#       br <- names(b)[1]
#       if ( is.null(b[[br]]) )
#         stop("no reference sequence available")
#     }
#     ## important: benchmark first due to --adjustdirection
#     gen[[i]] <- do.call(c, list(benchmark = as.list(b[[br]]), 
#                                 gen[[i]]))
    seqs <- mafft(seqs, method = "genafpair", maxiterate = 1000,
                      options = "--adjustdirectionaccurately", 
                      path = align.exe, quiet = TRUE)
    ## keep species names compatible with data base
    ## TO DO: add remark about reverse complement
#     rownames(gen[[i]]) <- gsub("_R_", "", rownames(gen[[i]]))
#     if ( is.list(gen[[i]]) ) gen[[i]] <- as.matrix(gen[[i]]) # e.g. Encalypta
#     if ( gen[[i]][1] == 0 ) gen[[i]] <- NULL
#     ## FILTER OUT WRONG SEQUENCES 
#     gen[[i]] <- filter.alignment(gen[[i]], x, logfile = logfile)
#     if ( all(rownames(gen[[i]]) %in% "benchmark") ){
#       
#     } else {
#       gen[[i]] <- deleteEmptyCells(gen[[i]][setdiff(rownames(gen[[i]]), "benchmark"),], quiet = TRUE)
#       d <- dist.dna(gen[[i]], model = "raw", pairwise.deletion = TRUE, 
#                     as.matrix = TRUE)
#       if ( !quiet ) slog(" --", nrow(gen[[i]]), 
#                          "species -- genetic distance",
#                          paste(round(range(d, na.rm = TRUE), 2), 
#                                collapse = "-"), 
#                          file = logfile)
#     }
  }
  else {
    seqs <- as.matrix(seqs)
    slog(" -- 1 species", file = logfile)  
  }
dbDisconnect(conn)
seqs
}