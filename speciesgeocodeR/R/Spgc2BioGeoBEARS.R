Spgc2BioGeoBEARS <- function(x, phyl = NULL, file = NULL, true.areas = T, 
                             true.species = T) {
  
  # reformat input
  dat <- x$spec_table
  rownames(dat) <- dat[, 1]
  dat <- data.frame(dat[, -1])
  
  # set occurrences to logical
  dat[dat >= 1] <- 1
  
  # check for redundant areas
  if (true.areas) {
    if (sum(colSums(dat) == 0) > 0) {
      warning(sprintf("%s areas with zeros species removed", 
                      sum(!colSums(dat) > 0)))
    }
    
    dat <- data.frame(dat[, colSums(dat) > 0])
  }
  
  # check for redundant species
  if (true.species) {
    if (sum(rowSums(dat) == 0) > 0) {
      warning(sprintf("%s species with zeros occurrences removed", 
                      sum(!rowSums(dat) > 0)))
    }
    
    dat <- data.frame(dat[rowSums(dat) > 0, ])
  }
  
  # treat phylogeny if provided
  if (!is.null(phyl)) {
    tr <- phyl
    data <- geiger::treedata(tr, dat, warnings = F)
    tr <- data$phy
    dat <- data$data
    ape::write.tree(tr, file = paste(file, ".tre", sep = ""))
  }
  
  # prepare output files
  length_dat <- length(dat[, 1])
  if (length(dat[1, ]) > 26) {
    stop(sprintf("%s areas; More than 26 exceed function limits", 
                 length(dat[1, ])))
  }
  letters_data <- paste("(", paste(letters[1:length(dat[1, ])], 
                                   collapse = " "), ")", sep = "")
  command_string <- paste(length_dat, length(dat[1, ]), letters_data, 
                          sep = "\t")
  
  # write output file
  if(!is.null(file)){
    write(command_string, file)
    write.table(dat, file = file, append = TRUE, quote = FALSE, sep = "", 
                col.names = FALSE, row.names = paste(rownames(dat), "\t"))
    
    names(dat) <- c(paste(letters[1:length(names(dat))], names(dat), 
                          sep = "_"))
  }
  
  return(list(BioGeoBEARS_command = command_string, BioGeoBEARS_matrix = dat))
} 