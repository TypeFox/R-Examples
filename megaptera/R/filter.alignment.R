filter.alignment <- function(DNAbin, megProj, step = "F", logfile){
  
  a <- DNAbin
  gene <- megProj@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  spec.tab <- paste("spec", gsub("^_", "", gene), sep = "_")
  threshold1 <- megProj@params@filter1
  threshold2 <- megProj@params@filter2
  threshold3 <- megProj@params@filter3
  threshold4 <- megProj@params@filter4
  
  dbc <- dbconnect(megProj@db)
  
  NNmean <- function(x, .mode = "min"){
    if ( all(is.na(x)) ){
      x <- NULL
    } else {
      x <- ifelse(.mode == "min",
                  min(x[!is.na(x)]),
                  median(x[!is.na(x)]))
    }
    x
  }
  
  ## FILTER 1: number of informative nucleotides per site
  ## -----------------------------------------------------
  m <- coverage(a, "nucleotides")
  aa <- a[, m > threshold1]
  #         write.nex(aa, paste(gene, "filter1.nex", sep = "."))
  
  ## FILTER 2: mean pairwise identity
  ## --------------------------------
  dd <- meanPairwiseIdentity(aa)
  #         write(dd, "meanPairwiseIdentity.txt")
  aaa <- deleteEmptyCells(aa[, dd > threshold2], quiet = TRUE)
  #         write.nex(aaa, paste(gene, "filter2.nex", sep = "."))
  s <- setdiff(rownames(aa), rownames(aaa))
  if ( length(s) > 0 ){
    slog("\n\n", length(s), "species with mean pairwise identity >", 
         threshold2, "excluded\n\t-",
         paste(s, collapse = "\n\t- "), file = logfile)
    SQL <- c(paste("UPDATE locus SET ", gene, "_blocks='excluded (step", 
                 step, "-filter2)' WHERE ", sql.wrap(s), sep = ""),
             paste("UPDATE ", spec.tab,
                 " SET block='excluded (step", step, "-filter2)'",
                 "WHERE ", sql.wrap(s), sep = ""))
    lapply(SQL, dbSendQuery, conn = dbc)
  }
  
  ## FILTER 3: mean number of oberlapping nucleotides
  ## ------------------------------------------------
  no <- nucleotideOverlap(aaa) # nn overlapping nucs
  rs <- sort(unlist(apply(no, 1, NNmean, .mode = "mean")))/ncol(aaa) # ranked species
  #         write(rs, "meanPairwiseOverlap.txt")
  id <- names(rs)[rs < threshold3]
  if ( length(id) > 0 ){
    slog("\n\n", length(id), "species with <", round(threshold3 * 100, 0), 
         "% mean overlapping nucleotides excluded\n\t-",
         paste(id, collapse = "\n\t- "), file = logfile)
    SQL <- c(paste("UPDATE locus SET ", gene, "_blocks='excluded (step", 
                   step, "-filter3)' WHERE ", sql.wrap(id), sep = ""),
             paste("UPDATE ", spec.tab,
                   " SET block='excluded (step", step, "-filter3)'",
                   "WHERE ", sql.wrap(id), sep = ""))
    lapply(SQL, dbSendQuery, conn = dbc)
  }
  aaa <- deleteEmptyCells(aaa[names(rs)[rs >= threshold3], ],
                            quiet = TRUE)
  #         write.nex(aaa, paste(gene, "filter3.nex", sep = "."))
  
  ## FILTER 4: mean distance
  ## ------------------------------------------------
  d <- dist.dna(aaa, model = "N", pairwise.deletion = TRUE, 
                as.matrix = TRUE)
  d <- d/ncol(aaa)
  diag(d) <- NA
  ee <- sort(unlist(apply(d, 1, NNmean, .mode = "min")))
  #         write(dd, "meanPairwiseDistance.txt")
  bad <- names(ee)[ee >= threshold4]
  good <- names(ee)[ee < threshold4]
  ## all sequences are bad:
  if ( "benchmark" %in% bad ) {
    ph <- bad; bad <- good; good <- ph
  }
  ## ignore threshold4 if exceeded by all species in alignment
  if ( length(good) == 0 ){
    slog("\n CAUTION: threshold4 ignored", file = logfile)
    ph <- bad; bad <- good; good <- ph
  }
  if ( length(bad) > 0 ){
    slog("\n\n", length(bad), "species isolated by a distance >", 
         round(threshold4, 5), "excluded\n\t-",
         paste(bad, collapse = "\n\t- "), file = logfile)
    bad <- gsub("_R_", "", bad)
    SQL <- c(paste("UPDATE locus SET ", gene, "_blocks='excluded (step", 
                   step, "-filter4)' WHERE ", sql.wrap(bad), sep = ""),
             paste("UPDATE ", spec.tab,
                   " SET block='excluded (step", step, "-filter4)'",
                   "WHERE ", sql.wrap(bad), sep = ""))
    lapply(SQL, dbSendQuery, conn = dbc)
  }
  dbDisconnect(dbc)
  DNAbin[good, ]
  #         write.nex(a, paste(gene, "filter4.nex", sep = "."))
}
