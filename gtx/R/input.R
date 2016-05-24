read.snpdata.mach <- function(fileroot, tol.af = 0.01, phenotypes = NULL, isuffix = ".mlinfo", dsuffix = ".mldose") {
  if (substr(isuffix, nchar(isuffix) - 2, nchar(isuffix)) == ".gz") {
    mlinfo <- read.table(gzfile(paste(fileroot, isuffix, sep = "")), header = TRUE, colClasses = "character")
  } else {
    mlinfo <- read.table(paste(fileroot, isuffix, sep = ""), header = TRUE, colClasses = "character")
  }
  ## read info as characters to prevent small .mlinfo files with Al1 or Al2 all "T" alleles being coerced to logical
  names(mlinfo) <- sub("^AL", "Al", names(mlinfo)) # some old files have AL1,AL2
  stopifnot(all(c("SNP", "Al1", "Al2", "Freq1") %in% names(mlinfo)))
  for (colname in intersect(c("Freq1", "MAF", "Quality", "Rsq"), names(mlinfo))) mlinfo[[colname]] <- as.double(mlinfo[[colname]]) # if present, coerce these columns to double
  ## we could use something like try(as.numeric(... to make this more general
  if (substr(dsuffix, nchar(dsuffix) - 2, nchar(dsuffix)) == ".gz") {
    mldose <- read.table(gzfile(paste(fileroot, dsuffix, sep = "")), header = FALSE, col.names = c("MACHID", "DATATYPE", paste(mlinfo$SNP, mlinfo$Al1, sep = "_")), as.is = TRUE)
  } else {
    mldose <- read.table(paste(fileroot, dsuffix, sep = ""), header = FALSE, col.names = c("MACHID", "DATATYPE", paste(mlinfo$SNP, mlinfo$Al1, sep = "_")), as.is = TRUE)
  }
  if (nrow(mlinfo) > 0 && any(abs(apply(mldose[ , -2:-1, drop = FALSE], 2, mean)/2 - mlinfo$Freq1) > tol.af)) stop("mlinfo/mldose allele frequency mismatch") # check mean dose corresponds to Freq1 from mlinfo file; we allow empty mlinfo/mldose files
  ## we should check for MACH-style names with "->" separator
  mldose$FID <- sapply(mldose$MACHID, function(ii) unlist(strsplit(ii, "->"))[1])
  mldose$IID <- sapply(mldose$MACHID, function(ii) unlist(strsplit(ii, "->"))[2]) # assumes two part names
  snpinfo <- subset(mlinfo, select = c("SNP", "Al1", "Al2", "Freq1"))
  names(snpinfo) <- c("snp", "coded.allele", "noncoded.allele", "coded.freq")
  ## note Al1 is the coded.allele in the sense that mldose contains does of Al1
  ## however mach2qtl inverts effect directions
  snpinfo <- cbind(snpinfo, subset(mlinfo, select = setdiff(names(mlinfo), c("SNP", "Al1", "Al2", "Freq1")))) # keep any other columns, user may have to coerce additional mini-mac columns back to double
  if (!is.null(phenotypes)) {
    stopifnot(is.data.frame(phenotypes))
    stopifnot("MACHID" %in% names(phenotypes))
    mldose <- cbind(mldose, subset(phenotypes, select = setdiff(names(phenotypes), "MACHID"))[match(mldose$MACHID, phenotypes$MACHID), ])
  }
  snpdata <- list(snpinfo = snpinfo, data = mldose)
  class(snpdata) <- "snpdata"
  return(snpdata)
}

read.snpdata.minimac <- function(fileroot, tol.af = 0.01, phenotypes = NULL, isuffix = ".info.gz", dsuffix = ".dose.gz") {
  return(read.snpdata.mach(fileroot, tol.af, phenotypes, isuffix, dsuffix))
}

read.snpdata.plink <- function(fileroot, tol.af = 0.01, phenotypes = NULL) {
    frq <- read.table(paste(fileroot, "frq", sep = "."), header = TRUE, colClasses = "character") # read as characters to prevent small .frq files with A1 or A2 all "T" alleles being coerced to logical
  stopifnot(all(c("SNP", "A1", "A2", "MAF") %in% names(frq)))
  for (colname in intersect(c("MAF", "NCHROBS"), names(frq))) frq[[colname]] <- as.double(frq[[colname]]) # if present, coerce these columns to double
  raw <- read.table(paste(fileroot, "raw", sep = "."), header = TRUE, as.is = TRUE, check.names = FALSE) # do not check.names since metabochip names like "chr1:11717263" will be broken
  stopifnot(all(paste(frq$SNP, frq$A1, sep = "_") == names(raw)[-6:-1])) # assume 6 columns before genotypes
  if (nrow(frq) > 0 && any(abs(frq$MAF - apply(raw[ , -6:-1], 2, mean, na.rm = TRUE)/2) > tol.af)) stop("frq/raw allele frequency mismatch") # check mean dose corresponds to MAF from frq file
  if (ncol(raw) > 6) for (colid in 7:ncol(raw)) raw[is.na(raw[ , colid]), colid] <- mean(raw[ , colid], na.rm = TRUE) # replace NA genotypes with column mean
  snpinfo <- subset(frq, select = c("SNP", "A1", "A2", "MAF"))
  names(snpinfo) <- c("snp", "coded.allele", "noncoded.allele", "coded.freq")                    
  snpinfo <- cbind(snpinfo, subset(frq, select = setdiff(names(frq), c("SNP", "A1", "A2", "MAF")))) # keep any other columns
  if (!is.null(phenotypes)) {
    stopifnot(is.data.frame(phenotypes))
    stopifnot(all(c("FID", "IID") %in% names(phenotypes)))
    raw <- cbind(raw, subset(phenotypes, select = setdiff(names(phenotypes), c("FID", "IID")))[match(paste(raw$FID, raw$IID, sep = "->"), paste(phenotypes$FID, phenotypes$IID, sep = "->")), ])
  }
  snpdata <- list(snpinfo = snpinfo, data = raw)
  class(snpdata) <- "snpdata"
  return(snpdata)
}

read.snpdata.impute <- function(samplefile, genofile, phenotypes = NULL) {
  ## hack to work around 0 0 0 P ... line
  samples <- read.table(samplefile, skip = 2, header = FALSE, as.is = TRUE, 
                        col.names = names(read.table(samplefile, header = TRUE, nrows = 0)))
  nsamples <- nrow(samples)
  ## hope we have right number of colClasses
  gts <- read.table(genofile, header = FALSE,
                    colClasses = c(rep("character", 5), rep("numeric", nrow(samples)*3)))
  stopifnot(ncol(gts) == 5 + nsamples*3)
  nsnps <- nrow(gts)
  snpinfo <- gts[ , 1:5]
  names(snpinfo) <- c("chr", "snp", "pos", "noncoded.allele", "coded.allele")
  i0 <- seq(from = 1, by = 3, length.out = nsamples)
  i1 <- seq(from = 2, by = 3, length.out = nsamples)
  i2 <- seq(from = 3, by = 3, length.out = nsamples)
  for (snpidx in 1:nsnps) {
    pvec <- t(gts[snpidx, 6:ncol(gts)]) # this traverses memory with large stride so only do once
    p0 <- pvec[i0, ]
    p1 <- pvec[i1, ]
    p2 <- pvec[i2, ]
    psum <- p0+p1+p2
    dose <- ifelse(psum > 0, (p1+2*p2)/psum, NA)
    mdose <- mean(dose, na.rm = TRUE)
    snpinfo$coded.freq[snpidx] <- mdose/2
    dose[is.na(dose)] <- mdose
    samples[[paste(snpinfo$snp[snpidx], snpinfo$coded.allele[snpidx], sep = "_")]] <- dose
  }
  # check on dosages on [0,2]?
  if (!is.null(phenotypes)) {
    stopifnot(is.data.frame(phenotypes))
    stopifnot(all(c("FID", "IID") %in% names(phenotypes)))
    samples <- cbind(samples, subset(phenotypes, select = setdiff(names(phenotypes), c("FID", "IID")))[match(paste(samples[ , 1], samples[ , 2], sep = "->"), paste(phenotypes$FID, phenotypes$IID, sep = "->")), ])
  }
  snpdata <- list(snpinfo = snpinfo, data = samples)
  class(snpdata) <- "snpdata"
  return(snpdata)
}

allelesAB <- function(A1, A2, sep = "/") return(ifelse(A1 < A2, paste(A1, A2, sep = sep), paste(A2, A1, sep = sep)))

align.snpdata.coding <- function(params, snpdata, ploidy = 2, missing.snp = "fail") {
  stopifnot(is.data.frame(params))
  stopifnot(all(c("snp", "coded.allele", "noncoded.allele") %in% names(params)))
  ## snpdata is list with snpdata$snpinfo and snpdata$data returned by read.snpdata.plink or read.snpdata.mldose
  stopifnot(all(c("snp", "coded.allele", "noncoded.allele", "coded.freq") %in% names(snpdata$snpinfo)))
  if (missing.snp == "okay") {
    snpinfo <- data.frame(snp = setdiff(params$snp, snpdata$snpinfo$snp))
    if (nrow(snpinfo) > 0) {
      snpinfo$coded.allele <- params$coded.allele[match(snpinfo$snp, params$snp)]
      snpinfo$noncoded.allele <- params$noncoded.allele[match(snpinfo$snp, params$snp)]
      snpinfo$coded.freq <- 0
      for (extracol in setdiff(names(snpdata$snpinfo), names(snpinfo))) snpinfo[[extracol]] <- NA
      snpdata$snpinfo <- rbind(snpdata$snpinfo, snpinfo)
      for (codeB in paste(snpinfo$snp, snpinfo$coded.allele, sep = "_")) snpdata$data[[codeB]] <- 0
    }
  }
  stopifnot(all(params$snp %in% snpdata$snpinfo$snp))
  ### non-robust workaround to deal with PLINK calling coded allele "0"
  for (idx in which(snpdata$snpinfo$snp %in% params$snp & snpdata$snpinfo$coded.allele == "0")) {
    jdx <- match(snpdata$snpinfo$snp[idx], params$snp)
    if (snpdata$snpinfo$noncoded.allele[idx] == params$coded.allele[jdx]) {
      names(snpdata$data) <- sub(paste("^", snpdata$snpinfo$snp[idx], "_", snpdata$snpinfo$coded.allele[idx], "$", sep = ""), paste(params$snp[jdx], params$noncoded.allele[jdx], sep = "_"), names(snpdata$data))
      snpdata$snpinfo$coded.allele[idx] <- params$noncoded.allele[jdx]
    } else if (snpdata$snpinfo$noncoded.allele[idx] == params$noncoded.allele[jdx]) {
      names(snpdata$data) <- sub(paste("^", snpdata$snpinfo$snp[idx], "_", snpdata$snpinfo$coded.allele[idx], "$", sep = ""), paste(params$snp[jdx], params$coded.allele[jdx], sep = "_"), names(snpdata$data))
      snpdata$snpinfo$coded.allele[idx] <- params$coded.allele[jdx]
    } else {
      stop(paste("cannot match alleles for", snpdata$snpinfo$snp[idx]))
    }
  }
  stopifnot(all(allelesAB(params$coded.allele, params$noncoded.allele) == allelesAB(snpdata$snpinfo$coded.allele, snpdata$snpinfo$noncoded.allele)[match(params$snp, snpdata$snpinfo$snp)])) # check alleles match even if inverted 
  params$codeB <- paste(params$snp, params$coded.allele, sep = "_")
  params$codeA <- paste(params$snp, params$noncoded.allele, sep = "_")
  ## the next section recodes columsn of snpdata$data to coded allele of params if necessary
  uset <- unique(subset(params, select = c("codeB", "codeA")))
  for (idx in which(!uset$codeB %in% names(snpdata$data))) {
    if (uset$codeA[idx] %in% names(snpdata$data)) {
      snpdata$data[[uset$codeB[idx]]] <- ploidy - snpdata$data[[uset$codeA[idx]]]
    } else {
      stop("cannot find dose for ", uset$codeB[idx], "/", uset$codeA[idx])
    }
  }
  stopifnot(all(params$codeB %in% names(snpdata$data)))
  rm(uset)
  for (colid in which(names(snpdata$data) %in% params$codeB)) {
    if (any(is.na(snpdata$data[ , colid]))) stop("missing dosages for", names(snpdata$data)[colid])
    if (any(snpdata$data[ , colid] < 0 | snpdata$data[ , colid] > ploidy)) stop("dosages outside [0,ploidy] for", names(snpdata$data)[colid])
  }
  params$data.coded.freq <- apply(subset(snpdata$data, select = params$codeB), 2, mean)/ploidy
  ## do some kind of check...
  return(list(params = params, snpdata = snpdata))
}
