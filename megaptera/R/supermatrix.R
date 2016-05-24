# CONCATENATE LOCI
# PACKAGE: megaptera
# AUTHOR: Christoph Heibl (at gmx.net)
# LAST CHANGE: 2014-03-12
# TO DO: marker-wise deletion of species (line 38)

supermatrix <- function(conn, subset, outgroup, min.n.seq, exclude){
  
  ## join taxonomy and locus tables
  ## ------------------------------
  tab <- "SELECT * FROM taxonomy INNER JOIN locus USING (spec)"
  tab <- dbGetQuery(conn, tab)
  tab[is.na(tab)] <- "excluded"
  cols <- grep("_blocks", names(tab))
  tab <- tab[apply(tab[, cols], 1, 
          function(x) length(grep("excluded", x)) != length(x)), ]
  
  if ( !missing(exclude) ){
    tab <- tab[!tab$spec %in% exclude, ]
  }
    
  ## subset
  ## ------
  if ( !missing(subset) | !missing(outgroup) ){
    id.subset <- grep(paste(subset[-1], collapse = "|"), tab[, subset[1]])
    if ( length(id.subset) == 0 ) stop("subset empty")
    id.outgroup <- grep(paste(outgroup[-1], collapse = "|"), tab[, outgroup[1]])
    if ( length(id.outgroup) == 0 ) stop("outgroup empty")
    tab <- tab[ union(id.subset, id.outgroup), ]
  }
  
  ## determine tables that contain more than min.n.seq species
  ## ---------------------------------------------------------
  if ( missing(min.n.seq) ){
    tabnames <- names(tab[cols])
  } else {
    tabnames <- apply(tab[, cols], 2, function(x) length(grep("selected", x)))
    tabnames <- names(tabnames[tabnames >= min.n.seq])
  }
  tabnames <- paste("spec", gsub("_blocks", "", tabnames), sep = "_")
  tabnames <- gsub("__", "_", tabnames)
  
  ## select and read individual alignments
  ## -------------------------------------
  x <- lapply(tabnames, dbReadDNA, 
              conn = conn, taxon = ".+", regex = TRUE,
              ignore.excluded = TRUE)
  names(x) <- gsub("^spec_", "", tabnames)
  ngene <- length(x)
  
  ## subset
  ## ------
  subset.alignment <- function(a, s) deleteEmptyCells(a[rownames(a) %in% s, ])
  x <- lapply(x, subset.alignment, s = tab$spec)
  
# -->> add here markerwise deletion of species
  
  ## partitions file for RAxML
  ## -------------------------
  p <- cbind(rep(1, length(x)), sapply(x, ncol))
  for ( i in 2:nrow(p) ){
    p[i, 1] <- p[i - 1, 2] + 1
    p[i, 2] <- p[i, 1] + p[i, 2] -1
  }
  p <- paste("DNA, ", rownames(p), " = ", p[, 1], "-", p[, 2], sep = "")
  write(p, "partitions.txt")
  
  ## create SUPERMATRIX
  ## ------------------
  x <- do.call(cbind.DNAbin, c(x, fill.with.gaps = TRUE))
  
  ## write data as PHY and NEX
  ## -------------------------
  fn <- paste("supermatrix", nrow(x), ngene, ncol(x), sep = "_")
  write.nex(x, paste(fn, "nex", sep = "."))
  write.phy(x, paste(fn, "phy", sep = "."))
  files <- c(paste(fn, "phy", sep = "."), "partitions.txt")
  
  ## outgroup
  ## ---------
  if ( !missing(outgroup)){
    ogp <- tab[tab[outgroup[1]] == outgroup[2], "spec"]
    og <- paste(ogp, collapse = ",")
    clip <- pipe("pbcopy", "w")
    write(og, file = clip)
    close(clip)
    write(og, "outgroup.txt")
    files <- c(files, "outgroup.txt")
  } else {
    ogp <- NULL
  }
  
  ## zip and return
  ## --------------
  zip(zipfile = fn, files = files)
  list(zipname = fn,
       supermatrix = x,
       outgroup = ogp,
       partitions = p)
}

