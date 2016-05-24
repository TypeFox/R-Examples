# Clean up shared file into a traditional OTU table

# 1. remove label, numOTus
# 2. transpose
# 3. shorten sample names

clean_shared <- function(shared_file, trim.names){
  x <- read.table(shared_file, header = T, row.names = 2)
  x$label <- NULL
  x$numOtus <- NULL
  if(trim.names == T){
    split.names <- strsplit(rownames(x), split = "\\.")
    samples <- c()
    for(i in 1:length(split.names)){
      samples[i] <- split.names[[i]][1]
    }
    rownames(x) <- make.unique(samples)
  }
  x <- t(x)
  return(x)
}

# Keep only the OTU column of the taxonomy file, and only OTUs that survived the subsampling step in the OTU table
clean_taxonomy <- function(taxonomy_file, table, remove_bootstrap){
  y <- read.csv(taxonomy_file, header = T)
  keep <- match(rownames(table), y$seqID, nomatch = NA)
  y <- y[keep,]
  colnames(y) <- c("OTU", "Kingdom", "Phylum", "Class", "Order", "Lineage", "Clade", "Tribe")
  if(remove_bootstrap == T){
    y$Kingdom <- gsub("\\(\\d*\\)", "", y$Kingdom)
    y$Phylum <- gsub("\\(\\d*\\)", "", y$Phylum)
    y$Class <- gsub("\\(\\d*\\)", "", y$Class)
    y$Order <- gsub("\\(\\d*\\)", "", y$Order)
    y$Lineage <- gsub("\\(\\d*\\)", "", y$Lineage)
    y$Clade <- gsub("\\(\\d*\\)", "", y$Clade)
    y$Tribe <- gsub("\\(\\d*\\)", "", y$Tribe)
  } 
  return(y)
}

