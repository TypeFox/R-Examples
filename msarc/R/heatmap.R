# Generate heat maps showing similarity of mass spec experiments

# Overall method:
#   1) load several mass spec experiments
#   2) Extract set of all observed Uniprot IDs
#   3) Assign score to each ID per each file
#      a) either just binary (present/absent)
#      b) or by Mascot score (-1 == missing)
#   4) Generate clustering/heat map of resulting table, somehow or other.

makeTable <- function(samples,useScore=F) {
  ids <- vector()
  for (sample in samples) {
    ids <- c(ids,sample$data$uniprot)
  }
  ids <- unique(ids)
  tbl <- data.frame(ids)
  rownames(tbl) <- ids
  for (sample in samples) {
    if (useScore) {
      map <- match(ids,sample$data$uniprot)
      tbl[[sample$filename]] <- sample$data$score[map]
      if (sum(is.na(map)) > 0) {
        tbl[[sample$filename]][is.na(tbl[[sample$filename]])] <- 0
      }
    } else {
      tbl[[sample$filename]] <- as.numeric(ids %in% sample$data$uniprot)
    }
  }
  tbl$ids <- NULL
  return(tbl)
}

msarc.plotHeatmap <- function(msalist,method="euclidean",useScore=T,...) {
  tbl <- makeTable(msalist,useScore=useScore)
  tbl <- data.matrix(tbl)
  d = dist(t(tbl),diag=T,upper=T,method=method)
  heatmap.2(as.matrix(d),margins=c(14,14),...)
  return(invisible(d))
}
