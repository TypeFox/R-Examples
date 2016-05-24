#' snp2gen function
#'
#' Convert genotypes in a snp matrix to genepop genotype files.
#'
#' Kevin Keenan, QUB, 2014
snp2gen <- function(infile = NULL, prefix_length = 2){
  #infile <- "snpmat.txt"
  #prefix_length = 6
  if(is.null(infile)){
    stop("Please provide an input file!")
  }
  fastScan <- function(fname) {
    s <- file.info(fname)$size
    buf <- readChar(fname, s, useBytes = TRUE)
    # replace Mac encoded line endings
    if(length(grep("\r", buf)) != 0L){
      buf <- gsub("\r", "\n", buf)
      buf <- gsub("\n\n", "\n", buf)
    }
    return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
  }
  # conver the vector to a character vector
  if(!is.null(dim(infile))){
    infile <- as.matrix(infile)
    dat <- sapply(1:nrow(infile), function(i){
      paste(infile[i,], collapse = "\t")
    })
  } else {
    dat <- fastScan(infile)
    # deal with empty last lines in files
    if(dat[length(dat)] == ""){
      dat <- dat[-length(dat)]
    }
  }
  # now manipulate the dat character vector
  inds <- strsplit(dat[1], split = "\\s+")[[1]][-1]
  dat <- dat[-1]
  locs <- gsub("\\s+.*$", "", dat)
  # strip locus names and add a leading tab
  dat <- gsub("^\\S+\\s+", "", dat)
  # add a trailing tab
  dat <- paste0(dat, "\t")
  # transpose genotypes
  dat <- apply(do.call(cbind, strsplit(dat, split = "\\s+")), 1,
               paste, collapse = "\t")
  # create and index factor
  splitNames <- lapply(inds, function(x){
    return(strsplit(x, split = "")[[1]])
  })
  pop_pre <- sapply(splitNames, function(x){
    return(paste(x[1:prefix_length], collapse = ""))
  })
  prefixes <- unique(pop_pre)
  pop_idx <- lapply(prefixes, function(x){
    which(x == pop_pre)
  })
  npops <- length(prefixes)
  pop_sizes <- sapply(pop_idx, length)
  # define all possible genotypes
  bases <- c("A", "C", "G", "T")
  genos <- c("01", "02", "03", "04")
  nuc_gts <- apply(expand.grid(bases, bases), 1, paste,
                   collapse = "")
  gp_gts <- apply(expand.grid(genos, genos), 1, paste,
                  collapse = "")
  # define a replacement function
  gt_replace <- function(bases, gp, dat){
    gsub(bases, gp, dat)
  }
  # convert bases to genepop format
  for(i in 1:length(nuc_gts)){
    dat <- gsub(nuc_gts[i], gp_gts[i], dat)
  }
  # replace missing data
  dat <- gsub("--", "0000", dat)
  # strip leading an trailing whitespace
  dat <- gsub("^\\s+|\\s+$", "", dat)
  # paste individual names to each string
  dat <- paste(paste0(inds, " ,"), dat, sep = "\t")
  # extract populations
  dat <- lapply(pop_idx, function(x){
    return(c("POP", dat[x]))
  })
  dat <- c("snp2gen-converted", paste(locs, collapse = ", "), unlist(dat))
  dat <- paste0(dat, collapse = "\n")
  writeLines(dat, "snp2gen_converted.gen")
}