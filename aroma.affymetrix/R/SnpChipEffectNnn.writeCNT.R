setMethodS3("writeCNT", "SnpChipEffectFile", function(this, reference, filename=NULL, path=".", tags=NULL, fields=c("Log2Ratio", "FreqB"), chromosomes=NULL, digits=3, ..., overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'fields':
  fields <- match.arg(fields, several.ok=TRUE);

  # Arguments 'tags':
  if (is.null(tags)) {
    if (!all(fields == c("Log2Ratio", "FreqB")))
      tags <- c(tags, fields);
  }

  # Arguments 'filename':
  if (is.null(filename)) {
    fullname <- getFullName(this);
    fullname <- gsub(",chipEffects", "", fullname, fixed=TRUE);
    fullname <- paste(c(fullname, tags), collapse=",");
    filename <- sprintf("%s.cnt", fullname);
  }

  # Arguments 'filename' & 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=!overwrite);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Writing CNT file");
  verbose && cat(verbose, "Pathname: ", pathname);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  cdf <- getCdf(this);
  gi <- getGenomeInformation(cdf);
  allChromosomes <- getChromosomes(gi);
  if (is.null(chromosomes)) {
    chromosomes <- allChromosomes;
  } else {
    chromosomes <- intersect(chromosomes, allChromosomes);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Write CNT file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  con <- file(pathname, open="w");
  on.exit(close(con));

  verbose && enter(verbose, "Extract CNT data");
  data <- extractCNT(this, reference=reference, fields=fields, chromosomes=chromosomes, ..., verbose=less(verbose,5));
  if (nrow(data) > 0) {
    # Round signals
    cc <- which(sapply(data, FUN=is.double));
    data[,cc] <- round(data[,cc], digits=digits);
  }
  verbose && str(verbose, data);
  verbose && exit(verbose);

  # Write CNT header
  verbose && enter(verbose, "Writing CNT header");
  cat(file=con, "[Header]\n");
  chipType <- getChipType(cdf);
  chipType <- gsub(",monocell", "", chipType, fixed=TRUE);
  cat(file=con, sprintf("ChipType1=%s\n", chipType));
  cat(file=con, "[ColumnName]\n");
  cat(file=con, paste(colnames(data), collapse="\t"), "\n", sep="");
  verbose && exit(verbose);

  # Write CNT data
  verbose && enter(verbose, "Writing CNT data");
  cat(file=con, "[Data]\n");
  write.table(file=con, data, sep="\t", row.names=FALSE, col.names=FALSE, 
                                                              quote=FALSE);
  # Not needed anymore
  data <- NULL;
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(pathname);
}, protected=TRUE)




setMethodS3("writeCNT", "SnpChipEffectSet", function(this, reference, filename=NULL, path=".", tags=NULL, fields=c("Log2Ratio", "FreqB"), chromosomes=NULL, digits=3, ..., overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'fields':
  fields <- match.arg(fields, several.ok=TRUE);

  # Arguments 'tags':
  if (is.null(tags)) {
    if (!all(fields == c("Log2Ratio", "FreqB")))
      tags <- c(tags, fields);
  }

  # Arguments 'filename':
  if (is.null(filename)) {
    fullname <- getFullName(this);
    fullname <- paste(c(fullname, tags), collapse=",");
    filename <- sprintf("%s.cnt", fullname);
  }

  # Arguments 'filename' & 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=!overwrite);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Writing CNT file");
  verbose && cat(verbose, "Pathname: ", pathname);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  cdf <- getCdf(this);
  gi <- getGenomeInformation(cdf);
  allChromosomes <- getChromosomes(gi);
  if (is.null(chromosomes)) {
    chromosomes <- allChromosomes;
  } else {
    chromosomes <- intersect(chromosomes, allChromosomes);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Write CNT file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  con <- file(pathname, open="w");
  on.exit(close(con));

  # Write CNT header
  verbose && enter(verbose, "Writing CNT header");
  cat(file=con, "[Header]\n");
  chipType <- getChipType(cdf);
  chipType <- gsub(",monocell", "", chipType, fixed=TRUE);
  cat(file=con, sprintf("ChipType1=%s\n", chipType));
  cat(file=con, "[ColumnName]\n");
  data <- extractCNT(this, reference=reference, fields=fields, chromosomes=99, ...);
  cat(file=con, paste(colnames(data), collapse="\t"), "\n", sep="");
  # Not needed anymore
  data <- NULL;
  verbose && exit(verbose);

  # Write CNT data
  verbose && enter(verbose, "Writing CNT data");
  cat(file=con, "[Data]\n");
  for (chr in chromosomes) { 
    verbose && enter(verbose, sprintf("Chromosome %d", chr));
    data <- extractCNT(this, reference=reference, fields=fields, chromosomes=chr, ...);

    if (nrow(data) > 0) {
      # Round signals
      cc <- which(sapply(data, FUN=is.double));
      data[,cc] <- round(data[,cc], digits=digits);

      # Write to file
      write.table(file=con, data, sep="\t", row.names=FALSE, col.names=FALSE, 
                                                                quote=FALSE);
    }
    # Not needed anymore
    data <- NULL;
    verbose && exit(verbose);
  } # for (chr ...)
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(pathname);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2008-12-29
# o Added extractCNT() and writeCNT().
# o Created.
############################################################################
