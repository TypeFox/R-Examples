textRegion <- function(region=NULL, ...) {
  # Arguments 'region':
  if (is.null(region)) {
    query <- "Set region in Mb: ";
  } else {
    region <- Arguments$getDoubles(region, range=c(0,999), length=c(2,2));
    currRegion <- sprintf("%.2f-%.2f Mb", region[1], region[2]);
    query <- sprintf("Set region (now %s): ", currRegion);
  }

  newRegion <- NULL;
  while (length(newRegion) != 2) {
    newRegion <- textNumericVector(query=query, default=region);
  }

  newRegion;
} # textRegion()


textNumericVector <- function(query="Enter numeric vector: ", default=NULL, ...) {
  # Arguments 'default':
  if (!is.null(default)) {
    default <- Arguments$getNumerics(default, ...);
  }

  ans <- NULL;
  while (length(ans) == 0) {
    cat(query);
    ans <- textCharacterVector(query=query, default=default, ...);
    ans <- as.double(ans);
    ans <- ans[is.finite(ans)];
    ans0 <- ans;
    ans <- NULL;
    tryCatch({
      ans <- Arguments$getNumerics(ans0, ...);
    }, error = function(ex) {
      print(ex$message);
    });
  }

  ans;
} # textNumericVector()


textCharacterVector <- function(query="Enter character vector: ", default=NULL, ...) {
  # Arguments 'default':
  if (!is.null(default)) {
    default <- Arguments$getCharacters(default, ...);
  }

  ans <- NULL;
  while (length(ans) == 0) {
    cat(query);
    ans <- readline();
    if (identical(ans, "@@@")) {
      return(default);
    }
    ans <- strsplit(ans, split="[,]")[[1]];
    # Return default value
    ans <- sapply(ans, FUN=trim);
    ans0 <- ans;
    ans <- NULL;
    tryCatch({
      ans <- Arguments$getCharacters(ans0, ...);
    }, error = function(ex) {
      print(ex$message);
    });
  }

  ans;
} # textCharacterVector()


textPlotAnnotations <- function(pa=NULL, ...) {
  paDefault <- list(
    xScale = 1e-6,
    Clim = c(0,6),
    width = 640,
    height = 240,
    rawCol = "#999999",
    smoothCol = "NA",
    smoothCex = 1.5,
    showSampleName = FALSE,
    showChr = FALSE,
    xlab = "",
    ylab = ""
  );

  if (identical(pa, "default")) {
    return(paDefault);
  }

  if (is.null(pa)) {
    pa <- paDefault;
  }

  pa0 <- pa;

  while (TRUE) {
    str(pa);
    names <- names(pa);
    names(names) <- seq(along=names);
    choices <- c(names, r="RESET", s="SAVE", c="CANCEL");
    choice <- textMenu(choices, title="Select attribute:", value=FALSE);
    choiceStr <- names(choices)[choice];
    if (choiceStr == "c") {
      return(pa0);
    } else if (choiceStr == "s") {
      return(pa);
    } else if (choiceStr == "r") {
      pa <- paDefault;
    } else {
      name <- choices[choice];
      default <- pa[[name]];
      cat(sprintf("Current value (length=%d):\n", length(default)));
      print(default);
      length <- rep(length(default),2);
      if (is.character(default)) {
        query <- "Enter value: ";
        ans <- textCharacterVector(query=query, default=default, length=length);
      } else if (is.numeric(default)) {
        query <- "Enter a numeric vector: ";
        ans <- textNumericVector(query=query, default=default, length=length);
      } else if (is.logical(default)) {
        query <- "Enter logical value: ";
        ans <- textCharacterVector(query=query, default=default, length=length);
        ans <- as.logical(ans);
      }
      cat("New value:\n");
      print(ans);
      pa[[name]] <- ans;
    }
  } # while ()

  pa;
} # textPlotAnnotations()


extractTcgaSignals <- function(dsList, chromosome, region=NULL, kappa=1, ...) {
  ugp <- getAromaUgpFile(dsList$total);
  units0 <- getUnitsOnChromosome(ugp, chromosome=chromosome, region=region);

  # Extract allele B fractions
  ds <- dsList$fracB;
  dfT <- getFile(ds, 1);
  dfN <- getFile(ds, 2);

  # Identify SNPs only
  units <- units0;
str(units);
  betaN <- extractRawAlleleBFractions(dfN, chromosome=chromosome, units=units);
  units <- units[!is.na(getSignals(betaN))];
str(units);

  betaN <- extractRawAlleleBFractions(dfN, chromosome=chromosome, units=units);
  betaT <- extractRawAlleleBFractions(dfT, chromosome=chromosome, units=units);

  sampleName <- getName(dfN);
  betaN$name <- betaT$name <- sampleName;
  betaN$ylim <- betaT$ylim <- c(-0.1,1.1);

  # Extract total CNs
  ds <- dsList$total;
  dfT <- getFile(ds, 1);
  dfN <- getFile(ds, 2);
  thetaN <- extractRawCopyNumbers(dfN, chromosome=chromosome, units=units, logBase=NULL);
#print(thetaN$.yLogBase);
  thetaT <- extractRawCopyNumbers(dfT, chromosome=chromosome, units=units, logBase=NULL);
#print(thetaT$.yLogBase);
  thetaT$y <- 2 * thetaT$y/thetaN$y;
  C <- thetaT;

  C$name <- sampleName;
  C$col <- "black";
  C$ylim <- c(0,6);

  #Purify
  c <- C$y;
  bT <- betaT$y;
  cTA <- (1-bT)*c;
  cTB <- bT*c;
  bN <- betaN$y;
  cNA <- (1-bN)*2;
  cNB <- bN*2;
  cTA <- (cTA - (1-kappa)*cNA)/kappa;
  cTB <- (cTB - (1-kappa)*cNB)/kappa;
  cT <- cTA+cTB;
  bT <- cTB/cT;

  Cp <- clone(C);
  Cp$y <- cT;
  betaTp <- clone(betaT);
  betaTp$y <- bT;

  # Call genotypes
  muN <- callGenotypes(betaN);

  # TumorBoost
  betaTN <- normalizeTumorBoost(betaT, betaN=betaN, muN=muN);
  betaTpN <- normalizeTumorBoost(betaTp, betaN=betaN, muN=muN);

  # All beta tracks
  betaList <- list(betaT=betaT, betaN=betaN, betaTN=betaTN, 
                   betaTp=betaTp, betaTpN=betaTpN);

  # Add colors
  colorMap <- c(AA="red", AB="black", BB="red", "NA"="#999999");
  colorMap <- c(AA="#999999", AB="#000000", BB="#999999", "NA"="#cccccc");
  col <- getColors(muN, colorMap=colorMap);
  betaList <- lapply(betaList, FUN=function(beta) {
    beta$col <- col;
    addLocusFields(beta, "col");
    beta$yAt <- c(0,1/2,1);
    beta;
  })

  # Add heterozygous tracks
  hets <- which(isHeterozygous(muN));
  betaHetList <- lapply(betaList, FUN=extractSubset, hets);
  names(betaHetList) <- paste(names(betaHetList), "Het", sep="");
  betaList <- c(betaList, betaHetList);


  # Add rho tracks to each beta track
  rhoList <- lapply(betaList, FUN=function(beta) {
    rho <- extractRawMirroredAlleleBFractions(beta);
    rho$y <- 2*rho$y;
    rho$ylim <- c(0,1.1);
    rho;
  });
  names(rhoList) <- gsub("beta", "rho", names(rhoList));
  
  res <- list(C=C, Cp=Cp, muN=muN);
  res <- c(res, betaList, rhoList);

  res;
} # extractTcgaSignals()


extractTcgaPairs <- function(ds, ...) {
  library("aroma.tcga");
  library("gsubfn");
  ps <- BiospecimenCoreResource$getBarcodePatterns();
  pattern <- ps$aliqoutBarcode;
  setFullNamesTranslator(ds, function(names, ...) {
    strapply(names, pattern=pattern, FUN=function(...) {
      x <- unlist(list(...)); n <- length(x);
      sprintf("%s,%s", x[2], x[6]);
    });
  });

  # tumor-normal pairs
  names <- getNames(ds);
  t <- table(names);
  keep <- names[t == 2];
  dsT <- extract(ds, unique(keep));
  names <- getFullNames(dsT);
  dsT <- extract(dsT, order(names));

  dsT;
} # extractTcgaPairs()
