library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

csR <- AffymetrixCelSet$byName("GSE8605", chipType="Mapping10K_Xba142");
print(csR);

cesN <- doASCRMAv2(csR, arrays=1:6, drop=FALSE, verbose=verbose)$cesN;
print(cesN);

ceR <- getAverage(cesN, verbose=verbose);
ce <- cesN[[1]];

ugp <- getAromaUgpFile(csR);
print(ugp);

for (chr in getChromosomes(ugp)) {
  units <- getUnitsOnChromosome(ugp, chr);
  pos <- getPositions(ugp, units=units) / 1e6;
  thetaR <- extractTotalAndFreqB(ceR, units=units)[,"total"];
  data <- extractTotalAndFreqB(ce, units=units);
  data[,"total"] <- 2*data[,"total"] / thetaR;

  chrTag <- sprintf("Chr%02d", chr);
  toPNG(getFullName(cesN), tags=c(getFullName(ce), chrTag, "TCN"), aspectRatio=1, {
    layout(matrix(1:2, ncol=1));
    par(mar=c(3,4,2,1)+0.1, pch=".");

    # TCN tracks
    cn <- RawCopyNumbers(data[,"total"], pos, chromosome=chr);
    plot(cn, col="gray", cex=0.8, ylim=c(0,4));
    cnS <- gaussianSmoothing(cn, xOut=seq(xMin(cn), xMax(cn), by=1/2), sd=1);
    points(cnS, col="black");
    stext(side=3, pos=0, getName(ce));
    stext(side=3, pos=1, chrTag);

    # BAF tracks
    plot(pos, data[,"freqB"], cex=3, ylim=c(0,1));
    box(col="blue");
    stext(side=3, pos=0, getTags(cesN, collapse=","));
  });
} # for (chr ...)
