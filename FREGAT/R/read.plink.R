# A function from package 'snpStats' v.1.20.0 (c) 2015 David Clayton

read.plink <- function(bed, bim, fam, na.strings=c("0", "-9"), sep=".",
                       select.subjects=NULL, select.snps=NULL) {
  lb <- nchar(bed);
  ext <- substr(bed, lb-3, lb);
  if (ext==".bed") {
    stub <- substr(bed, 1, lb-4)
  } else {
    stub <- bed
    bed <- paste(bed, ".bed", sep="")
  }
  if (missing(bim))
    bim <- paste(stub, ".bim", sep="")
  if (missing(fam))
    fam <- paste(stub, ".fam", sep="")
  df.bim <- read.table(bim, comment.char="", as.is=TRUE, na.strings=na.strings)
  snps <- as.character(df.bim[,2])
  if (!is.null(select.snps)) {
    if (is.character(select.snps)) {
      select.snps <- match(select.snps, snps)
      if (any(is.na(select.snps)))
        stop("unrecognised snp selected")
    }
    else if (is.numeric(select.snps)) {
      if (min(select.snps)<1 || max(select.snps)>nrow(df.bim))
        stop("selected snp subscript out of range")
      select.snps <- as.integer(select.snps)
    }
    else
      stop("unrecognised snp selection")
    if (any(duplicated(select.snps)))
      stop("snp selection contains duplicates")
    select.snps <- sort(select.snps)
    df.bim <- df.bim[select.snps,]
    snps <- snps[select.snps]
  }     
  if (ncol(df.bim)==6) {
    names(df.bim) <- c("chromosome", "snp.name", "cM", "position",
                       "allele.1", "allele.2")
  }
  else {
    warning("non-standard .bim file")
  }
  rownames(df.bim) <- snps
  if (any(duplicated(snps)))
    stop("duplicated snp name(s)")
  df.fam <- read.table(fam, comment.char="", as.is=TRUE, na.strings=na.strings)
  if (!is.null(select.subjects)) {
    if (is.numeric(select.subjects)) {
      if (min(select.snps)<1 || max(select.snps)>nrow(df.fam))
        stop("selected subject subscript out of range")
      select.subjects <- as.integer(select.subjects)
    }
    else
      stop("unrecognised subject selection")
    if (any(duplicated(select.subjects)))
      stop("subject selection contains duplicates")
    select.subjects <- sort(select.subjects)
    df.fam <- df.bim[select.subjects,]
  }     
  ped <- as.character(df.fam[,1])
  mem <- as.character(df.fam[,2])
  if (any(duplicated(ped))) {
    if (any(duplicated(mem))) {
      id <- paste(ped, mem, sep=sep)
      if (any(duplicated(id)))
        stop("couldn't create unique subject identifiers")
    }
    else
      id <- mem
  }
  else
    id <- ped
  names(df.fam) <-  c("pedigree", "member", "father", "mother", "sex",
                      "affected")
  rownames(df.fam) <- id
  gt <- .Call("readbed", bed, id, snps, select.subjects, select.snps,
              PACKAGE="FREGAT")
  list(genotypes = SnpMatrix2numeric(gt), fam = df.fam, map = df.bim)
}

# same code, modified

read.plink.names <- function(bed, na.strings=c("0", "-9"), sep=".", ...) {
# dots to prevent 'unused argument' error

  lb <- nchar(bed);
  ext <- substr(bed, lb-3, lb);
  if (ext==".bed") {
    stub <- substr(bed, 1, lb-4)
  } else {
    stub <- bed
    bed <- paste(bed, ".bed", sep="")
  }
  bim <- paste(stub, ".bim", sep="")
  fam <- paste(stub, ".fam", sep="")

  if (!file.exists(bim)) stop("Cannot find 'genodata' file")

  df.bim <- read.table(bim, comment.char="", as.is=TRUE, na.strings=na.strings)
  snps <- as.character(df.bim[,2])

  df.fam <- read.table(fam, comment.char="", as.is=TRUE, na.strings=na.strings)
  ped <- as.character(df.fam[,1])
  mem <- as.character(df.fam[,2])
  if (any(duplicated(ped))) {
    if (any(duplicated(mem))) {
      id <- paste(ped, mem, sep=sep)
      if (any(duplicated(id)))
        stop("couldn't create unique subject identifiers")
    }
    else
      id <- mem
  }
  else
    id <- ped
  list(bed = bed, snvnames = snps, idnames = id)
}


read.plink.region <- function(bed, snps, id, select.snps) {
  if (is.logical(select.snps)) {
    select.snps <- which(select.snps)
  }
  if (min(select.snps)<1 || max(select.snps)>length(snps))
    stop("selected snp subscript out of range")
  select.snps <- as.integer(select.snps)
  select.snps <- sort(select.snps)
  snps <- snps[select.snps]
  .Call("readbed", bed, id, snps, NULL, select.snps,
              PACKAGE="FREGAT")
}
