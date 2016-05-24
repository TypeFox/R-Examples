funX <- function(megProj, acc.tab, taxa, max.bp, align.exe, logfile){
  
  conn <- dbconnect(megProj@db)
  seqs <- dbReadDNA(conn, acc.tab, taxa, max.bp, 
                    enforce.binomial = TRUE)
  slog(paste("\n-- ", ifelse(is.list(seqs), length(seqs), nrow(seqs)), 
             " seqs. of ", taxa,  sep = ""), file = logfile) 
  seqs <- mafft(seqs, method = "auto", path = align.exe)
  ## Vitis vinifera + trnLF:
  ## bad sequences make alignment longer than max.bp
  if ( ncol(seqs) > max.bp ){
    d <- dist.dna(seqs, model = "N", pairwise.deletion = TRUE, 
                  as.matrix = TRUE)
    n.zero <- apply(d, 2, function(x) length(x[x == 0]))
    n.zero <- which.max(n.zero)
    exclude <- rownames(d)[d[, n.zero] > 10]
    seqs <- deleteEmptyCells(seqs[!rownames(seqs) %in% exclude, ])
    
    exclude <- lapply(exclude, splitGiTaxon)
    exclude <- do.call(rbind, exclude)[, "gi"]
    SQL <- paste("UPDATE", acc.tab, 
                 "SET status = 'excluded (too distant)'",
                 "WHERE", sql.wrap(exclude, term = "gi"))
    dbSendQuery(conn, SQL)
    slog("\n-- NOTE:", length(exclude), "seqs. of", taxa, "excluded", file = logfile)
  }
  dbWriteDNA(conn, acc.tab, seqs, enforce.binomial = FALSE, 
             status = "aligned")
  dbDisconnect(conn)
}