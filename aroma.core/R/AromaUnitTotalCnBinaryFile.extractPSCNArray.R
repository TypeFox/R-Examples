setMethodS3("extractPSCNMatrix", "AromaUnitTotalCnBinaryFile", function(dfTCN, dfBAF="*", units=NULL, ..., verbose=FALSE) {
  # Argument 'dfTCN':
  dfTCN <- Arguments$getInstanceOf(dfTCN, "AromaUnitTotalCnBinaryFile");
  chipType <- getChipType(dfTCN);
  nbrOfUnits <- nbrOfUnits(dfTCN);

  # Argument 'dfBAF':
  if (identical(dfBAF, "*")) {
    # Automagically locate the BAF file
    path <- getPath(dfTCN);
    filename <- getFilename(dfTCN);
    filename <- gsub(",total", ",fracB", filename, fixed=TRUE);
    dfBAF <- AromaUnitFracBCnBinaryFile(filename, path=path);
  }

  dfBAF <- Arguments$getInstanceOf(dfBAF, "AromaUnitFracBCnBinaryFile");
  stopifnot(nbrOfUnits(dfBAF) == nbrOfUnits);
  stopifnot(getChipType(dfBAF) == chipType);

  # Argument 'units':
  if (is.null(units)) {
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits);
    nbrOfUnits <- length(units);
  }

  # Extract signals
  theta <- dfTCN[units,1,drop=TRUE];
  beta <- dfBAF[units,1,drop=TRUE];

  dimnames <- list(NULL, c("total", "fracB"));
  data <- matrix(c(theta, beta), nrow=nbrOfUnits, ncol=2L,
                             byrow=FALSE, dimnames=dimnames);

  data;
}) # extractPSCNMatrix()

setMethodS3("extractPSCNArray", "AromaUnitTotalCnBinaryFile", function(dfTCN, ..., drop=FALSE) {
  data <- extractPSCNMatrix(dfTCN=dfTCN, ...);
  if (!drop) {
    dimnames <- dimnames(data);
    dimnames[[3]] <- getName(dfTCN);
    dim(data) <- c(dim(data), 1L);
    dimnames(data) <- dimnames;
  }
  data;
})


setMethodS3("extractPSCNArray", "AromaUnitTotalCnBinarySet", function(dsTCN, dsBAF="*", ..., verbose=FALSE) {
  # Argument 'dsTCN':
  dsTCN <- Arguments$getInstanceOf(dsTCN, "AromaUnitTotalCnBinarySet");

  # Argument 'dsBAF':
  if (identical(dsBAF, "*")) {
  } else {
    dsBAF <- Arguments$getInstanceOf(dsBAF, "AromaUnitFracBCnBinaryFile");
    stopifnot(length(dsBAF) == length(dsTCN));
    stopifnot(all(getNames(dsBAF) == getNames(dsTCN)));
  }

  dfBAF <- "*";
  data <- NULL;
  nbrOfArrays <- length(dsTCN);
  for (ii in seq_len(nbrOfArrays)) {
    dfTCN <- dsTCN[[ii]];
    if (!identical(dsBAF, "*")) {
      dfBAF <- dsBAF[[ii]];
    }
    dataII <- extractPSCNArray(dfTCN, dfBAF, ..., drop=FALSE, verbose=verbose);
    if (is.null(data)) {
      naValue <- as.double(NA);
      dim <- dim(dataII);
      dimnames <- dimnames(dataII);
      dim[3] <- nbrOfArrays;
      dimnames[[3]] <- getNames(dsTCN);
      data <- array(naValue, dim=dim, dimnames=dimnames);
    }
    data[,,ii] <- dataII;
  } # for (ii ...)

  data;
})

###########################################################################
# HISTORY:
# 2011-11-11
# o Added extractPSCNArray() for AromaUnitTotalCnBinary{File|Set}.
# o Added extractPSCNMatrix().
###########################################################################
