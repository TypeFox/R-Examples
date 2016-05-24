library("aroma.affymetrix");
log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSetName <- "Affymetrix_2006-TumorNormal";
tags <- "ACC,-XY,BPN,-XY,RMA,FLN,-XY";
chipType <- "Mapping250K_Nsp";

pairs <- matrix(c(
  "CRL-2325D", "CRL-2324D",
  "CRL-5957D", "CRL-5868D",
  "CCL-256.1D", "CCL-256D",
  "CRL-2319D", "CRL-2320D",
  "CRL-2362D", "CRL-2321D",
  "CRL-2337D", "CRL-2336D",
  "CRL-2339D", "CRL-2338D",
  "CRL-2341D", "CRL-2340D",
  "CRL-2346D", "CRL-2314D"
), ncol=2, byrow=TRUE);
colnames(pairs) <- c("normal", "tumor");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("cesN")) {
  cdf <- AffymetrixCdfFile$byName(chipType);
  gi <- getGenomeInformation(cdf);

  cdfM <- getMonocellCdf(cdf);
  cesN <- SnpChipEffectSet$byName(dataSetName, tags=tags, cdf=cdfM);
  print(cesN);

  # Reorder arrays according to 'pairs' matrix
  cesN <- cesN[indexOf(cesN, pairs)];
  names <- getNames(cesN);
  dim(names) <- dim(pairs);
  print(names);


  Blim <- c(0,1);
  Clim <- c(0,4);
  Clab <- expression(C[T]/C[N]);
  BTlab <- expression(beta[T]);
  BNlab <- expression(beta[N]);
  Flab <- expression(beta[T]-beta[N]);
  Flim <- c(-1/2,1/2);
}

# Identify SNPs
chromosome <- 9
units <- getUnitsOnChromosome(gi, chromosome);
pos <- getPositions(gi, units=units)/1e6;
o <- order(pos);
units <- units[o];
pos <- pos[o];
str(units);

isSNP <- (getUnitTypes(cdf, units=units) == 2);
units <- units[isSNP];
pos <- pos[isSNP];

data <- extractTotalAndFreqB(cesN, units=units);
dim(data) <- c(dim(data)[1:2], dim(names));
dimnames(data) <- list(NULL, c("total", "freqB"), NULL, colnames(pairs));
C <- data[,"total",,"tumor"] / data[,"total",,"normal"];
B <- data[,"freqB",,];
F <- B[,,"tumor"] - B[,,"normal"];

F2 <- 2*(1/2-abs(B[,,"normal"]-1/2));
dAA <- sqrt((B[,,"normal"]-0)^2);
dBB <- sqrt((B[,,"normal"]-1)^2);
#F[dAA < 0.2 | dBB < 0.2] <- NA;


# Plot (x,C) and (x,F) along genome
devSet(2);
subplots(3, ncol=1);
par(mar=c(2,4,1,1)+0.1);
par(ask=TRUE);
for (kk in seq_len(ncol(C))) {
  plot(pos, C[,kk], ylim=Clim, ylab=Clab);
  abline(h=1, lty=3, col="red");
  stext(side=3, pos=0, paste(pairs[kk,], collapse="/"));
  stext(side=3, pos=1, sprintf("Chr %d", chromosome));
  plot(pos, B[,kk,"tumor"], ylim=Blim, ylab=BTlab);
#  smoothScatter(pos, B[,kk,"tumor"], bandwidth=c(1,0.0001), ylim=Blim, ylab=BTlab);
  plot(pos, F[,kk], cex=F2, ylim=Flim, ylab=Flab);
#  smoothScatter(pos, F[,kk], bandwidth=c(1,0.0001), ylim=Flim, ylab=Flab);
  abline(h=0, lty=3, col="red");
}


devSet(3);
BNlab <- expression(beta[N]);
BTlab <- expression(beta[T]);
unitNames <- getUnitNames(cdf, units=units);
col <- seq_len(ncol(C));
par(ask=TRUE);
for (kk in seq_len(nrow(C))) {
  plot(NA, xlim=Blim, ylim=Blim, xlab=BNlab, ylab=BTlab);
  stext(side=3, pos=0, unitNames[kk]);
  points(B[kk,,], col=col, pch=19);
}
