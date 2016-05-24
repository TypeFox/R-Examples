.valid.data <- function(data, is.OTU=TRUE, 
                        export.sampleIDs=FALSE) {
  if ( class(data) != "list" ) {
      stop("please provide ecology data sets as list, see ?RAM.input.formatting for details")
  } 
  if ( length(data) != length(names(data)) ) {
        stop("each datasets in the list should have a given name, see ?RAM.input.formatting for details")
  }
  samples <- list()
  for ( i in 1:length(data) ) {
    label <- names(data)[i]
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    if ( is.OTU ) {
      valid.OTU(elem)
      samples[[label]] <- names(elem)[-ncol(elem)]
    } else {
      samples[[label]] <- rownames(elem)
      if ( "taxonomy" %in% colnames(elem) ) {
         stop("are you sure this is NOT an OTU table?")
      }
    }
  }
  sampleids <- samples[[1]]
  if ( length(data) != 1 ) {
    for ( i in 2:length(data) ) {
      if ( !identical(sampleids, samples[[i]]) ) {
          stop("samples in each data set don't match, check sample# and order in your data sets before performing any analysis")
      }
    }
  } 
  if (export.sampleIDs) { return(sampleids) }
}


valid.OTU <- function(otu1, otu2=NULL) {
  
  given.both <- !is.null(otu2)

  # tests for otu1
  # ensure object is a data.frame
  if (!is.data.frame(otu1)) {
    stop("the given object for otu1 is not a data frame.")
  }
  
  tax <- dim(otu1)[2]
  
  # ensure the last column contains taxonomic data
  if (names(otu1)[tax] != "taxonomy") {
    stop("no 'taxonomy' column detected in otu1. Please ensure the last column of the table is titled 'taxonomy' (capitalization matters) and try again.")
  }
  
  # check if all data is numeric other than tax columns
  if (!all(apply(otu1[ ,-tax, drop=FALSE], is.numeric, MARGIN=1))) {
    stop("OTU data other than taxonomic information is not numeric in otu1. Please fix the data and try again.")
  }
  
  # check that all counts are non-negative
  if (any(otu1[ , -tax] < 0)) {
    stop("negative counts detected in otu1; cannot continue. Please ensure all counts are positive and try again.")
  }
  
  if (any(colSums(otu1[ ,-tax, drop=FALSE]) == 0)) {
    warning("some samples in otu1 have zero counts. Some functions may behave bizarrely as a result.")
  }
  
  missing.tax <- which(otu1[ , tax, drop=FALSE] == "")
  
  if (length(missing.tax) != 0) {
    missing.tax.percent <- 100 * (length(missing.tax) / dim(otu1)[1]) 
    warning(paste("taxonomic data is missing for ", format(missing.tax.percent, digits=3),
                  "% of OTUs in otu1.", sep=""))
  }
  
  # tests for otu2 data
  if (given.both) {
    
    if (!identical(names(otu1), names(otu2))) {
      stop("the samples for otu1 and otu2 do not match. Please ensure you have matching otu1 and otu2 data.")
    }
    
    # ensure object is a data.frame
    if (!is.data.frame(otu2)) {
      stop("the given object for otu2 is not a data frame.")
    }
    
    tax <- dim(otu2)[2]
    
    # ensure the last column contains taxonomic data
    if (names(otu2)[tax] != "taxonomy") {
      stop("no 'taxonomy' column detected in otu2. Please ensure the last column of the table is titled 'taxonomy' (capitalization matters) and try again.")
    }
    
    # check if all data is numeric other than tax columns
    if (!all(apply(otu2[ ,-tax, drop=FALSE], is.numeric, MARGIN=1))) {
      stop("OTU data other than taxonomic information is not numeric in otu2. Please fix the data and try again.")
    }
    
    # check that all counts are non-negative
    if (any(otu2[ , -tax, drop=FALSE] < 0)) {
      stop("negative counts detected in otu2; cannot continue. Please ensure all counts are positive and try again.")
    }
    
    if (any(colSums(otu2[ ,-tax, drop=FALSE]) == 0)) {
      warning("some samples in otu2 have zero counts. Some functions may behave bizarrely as a result.")
    }
    
    missing.tax <- which(otu2[ , tax, drop=FALSE] == "")
    
    if (length(missing.tax) != 0) {
      missing.tax.percent <- 100 * (length(missing.tax) / dim(otu2)[1]) 
      warning(paste("taxonomic data is missing for ", format(missing.tax.percent, digits=3),
                    "% of OTUs in otu2.", sep=""))
    }
  }
}

.valid.factors <- function(meta, factors, min.factors=1, max.factors=Inf) {
  
  if (max.factors < min.factors || max.factors < 1 || min.factors < 0) {
    stop("something has gone wrong internally; please contact the package maintainer.")
  }

  if (!is.data.frame(meta)) {
    stop("'meta' must be a data frame.")
  }
  
  if (!is.character(factors)) {
    stop("'factors' must be a character vector; see ?RAM.factors for help.")
  }
  
  factors.len <- length(factors)
  
  if (factors.len < min.factors) {
    stop(paste("only", factors.len, "factor(s) were given where", min.factors, "were needed."))
  }
  
  if (length(names(factors)) != factors.len) {
    warning("not all factors are named; converting the names of ALL factors to column names; see ?RAM.factors for help.")
    names(factors) <- factors
  }
  
  if (factors.len > max.factors) {
    warning(paste(factors.len, "factors were given, this function supports up to",
                  max.factors, "factor(s); ignoring the others."))
    factors <- factors[1:max.factors]
    names(factors) <- names(factors)[1:max.factors]
  }
  
  if (!all(factors %in% names(meta))) {
    stop("not all factors were found in 'meta', please check your spelling and try again.")
  }
  
  output <- meta[ , factors, drop=FALSE]
  names(output) <- names(factors)
  
  output
}

.valid.plot.settings <- function(file=NULL, ext=NULL) {
  # ignore NULL values (the default when unused)
  if (is.null(file) && is.null(ext)) {
    return()
  }
  
  # if only one of {ext, file} is specified, raise an error
  if (xor(is.null(ext), is.null(file))) {
    stop("please specify either both 'file' and 'ext', or neither.")
  }
  
  # validate the extension
  .valid.ext(ext)
  
  # validate file path; prepend the device with period (for valid file name)
  file <- .ensure.filepath(file, ext)
  
  file
}

.valid.ext <- function(ext) {
  valid.exts <- c("pdf", "png", "tiff", "svg", "bmp", "jpeg")
  
  if (!ext %in% valid.exts) {
    stop("invalid extension given. See ?RAM.plotting for more details.")
  }
}

.valid.rank <- function(rank) {
  if (!is.character(rank) || !identical(length(rank), 1L)) {
    stop("rank must be a character vector of length 1; see ?RAM.rank.formatting for help.")
  }
  
  tax.classes <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  tax.classes.short <- c("k", "p", "c", "o", "f", "g", "s")
  tax.classes.pattern <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  tax.classes.all <- c(tax.classes, tax.classes.short, tax.classes.pattern)
  
  if (!(rank %in% tax.classes.all)) {
    stop("invalid format for rank. See ?RAM.rank.formatting for help.")
  }
  
  invisible()
}

.valid.meta <- function(otu1, otu2=NULL, meta) {
  valid.OTU(otu1, otu2)
  
  # break when meta == NULL, since meta is an optional argument in some 
  # functions
  if (is.null(meta)) { return(invisible(TRUE)) }
  
  if (!is.data.frame(meta)) {
    stop("'meta' must be a data frame.")
  }
  
  otu1.samples <- names(otu1)[-dim(otu1)[2]]
  meta.samples <- rownames(meta)
  
  if (!identical(otu1.samples, meta.samples)) {
    in.otu1 <- setdiff(otu1.samples, meta.samples)
    in.meta <- setdiff(meta.samples, otu1.samples)
    
    error.msg <- paste("Sample names for otu1 and meta differ. They may be out of order.\nIn otu1, not in meta:\n", 
                       paste(in.otu1, collapse=" "), "\nIn meta, not in otu1:\n", 
                       paste(in.meta, collapse=" "))
    
    stop(error.msg)
  }
  
  if (!is.null(otu2)) {
    otu2.samples <- names(otu2)[-dim(otu2)[2]]
    
    if (!identical(otu2.samples, meta.samples)) {
      in.otu2 <- setdiff(otu2.samples, meta.samples)
      in.meta <- setdiff(meta.samples, otu2.samples)
      
      
      error.msg <- paste("Sample names for otu2 and meta differ. They may be out of order.\nIn otu2, not in meta:\n", 
                         paste(in.otu2, collapse=" "), "\nIn meta, not in otu2:\n", 
                         paste(in.meta, collapse=" "))
      
      stop(error.msg)
    }
  }
  
  invisible(TRUE)
}


.valid.labels <- function(num.otus, labels) {
  
  if (!is.numeric(num.otus) || length(num.otus) != 1L) {
    stop("internal error: 'num.otus' must be a numeric vector of length one. If you can see this message please contact the package maintainer.")
  }
  
  if (!is.character(labels)) {
    stop("'labels' must be a character vector.")
  }
  
  num.labels <- length(labels)
  
  if (num.labels < num.otus) {
    stop(paste("only", num.labels, "label was given for", num.otus, "OTU tables."))
  }
  
  if (num.labels != num.otus) {
    warning(paste("there were", num.otus, "OTU table(s) given, but", 
                  num.labels, "label(s) given. Other labels will be ignored."))
  }
}


.valid.factor <- function(meta, meta.factor){
  if ( !any(meta.factor %in% names(meta)) ) {
     stop("none of the variables are in the metadata")
  } else {
    vec <- vector()
    for ( i in meta.factor) {
      if ( !i %in% names(meta) ) {
        vec <- c(vec, i)
      }
    }
  }
  vec <- unique(vec)
  if ( length(vec) !=0  ) {
    warning(paste("the following variables are not in the metadata: ", paste(vec, collapse=", "), sep=""))
    factors <- meta.factor[-which(meta.factor %in% vec)]
  } else {
    factors <- meta.factor
  }
  return(factors)
}     
