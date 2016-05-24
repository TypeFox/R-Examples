# Data processing functions (post-mothur)

# Select data from a single bog lake and layer
bog_subset <- function(bog_id, table){
  sample_ids <- colnames(table)
  bog_names <- substr(sample_ids, start = 1, stop = 3)
  subset <- table[, grep(bog_id, bog_names, ignore.case = T)]
  return(subset)
}

# Select data from a single year
year_subset <- function(year_id, table){
  sample_ids <- colnames(table)
  year_names <- substr(sample_ids, start = 9, stop = 10)
  subset <- table[, grep(year_id, year_names)]
  return(subset)
}

# Create a date vector from a vector of sample names
extract_date <- function(sample_ids){
  dates <- as.Date(substr(sample_ids, start = 4, stop = 10), format = "%d%B%y")
  return(dates)
}

# Remove replicate samples
remove_reps <- function(table){
  sample_ids <- colnames(table)
  rep_ids <- substr(sample_ids, start = 12, stop = 13)
  subset <- table[, which(rep_ids == "" | rep_ids == "R1")]
  return(subset)
}

# Select data from a taxonomic group
grab_group <- function(group, level, table, taxonomy){
  column <- match(level, colnames(taxonomy))
  search <- grep(group, taxonomy[, column])
  subset <- table[search, ]
  names <- c()
  for(i in 1:length(search)){
    names[i] <- paste(taxonomy[search[i], ], collapse = ";")
  }
  rownames(subset) <- make.unique(names)
  return(subset)
}

# Create a table with groups at a higher taxonomic level than OTUs
combine_otus <- function(level, table, taxonomy){
  key <- match(level, colnames(taxonomy))
  column <- c()
  for(i in 1:dim(taxonomy)[1]){
    column[i] <- paste(taxonomy[i, 1:key], collapse = ";")
  }
  unique_groups <- unique(column)
  new_table <- rep(NA, dim(table)[2])
  for(i in 1:length(unique_groups)){
    members <- which(column == unique_groups[i])
    if(length(members) > 1){
      member_abun <- colSums(table[members, ])
      new_table <- rbind(new_table, member_abun)
    }else{
      member_abun <- table[members,]
      new_table <- rbind(new_table, member_abun)
    }
  }
  new_table <- new_table[2:dim(new_table)[1], ]
  rownames(new_table) <- unique_groups
  colnames(new_table) <- colnames(table)
  new_table <- data.frame(new_table)
  return(new_table)
}

# Reduce full taxonomic name to only the last classified value in tables at higher taxonomic levels than OTUs
reduce_names <- function(table){
  long_names <- rownames(table)
  split_names <- strsplit(long_names, split = ";")
  short_names <- c()
  for(i in 1:length(split_names)){
    taxa <- split_names[[i]]
    unknown <- grep("unclassified", taxa)
    if(length(unknown) == 0){
      blanks <- grep("__$", taxa)
      short_names[i] <- taxa[length(taxa) - length(blanks)]
    }else if(length(taxa) == length(unknown)){
      short_names[i] <- "unclassified"
    }else{
      short_names[i] <- taxa[length(taxa) - length(unknown)]
    }
  }
  rownames(table) <- make.unique(short_names)
  return(table)
}


# Convert a vector of relative abundance data for a taxon into z-score normalized form.
zscore <- function(table){
  ztable <- matrix(NA, ncol=length(table[1,]), nrow=length(table[,1]))
  for (i in 1:length(table[,1])){
    otu <- as.numeric(table[i,])
    zotu <- (otu - mean(otu)) / sd(otu)
    ztable[i,] <- zotu
  }
  rownames(ztable) <- rownames(table)
  colnames(ztable) <- colnames(table)
  return(ztable)
}
