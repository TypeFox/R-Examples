##########################################################################
# Data set:
# GSE8605/
#   Mapping10K_Xba142/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
library("aroma.affymetrix");
library("ACNE");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

cdf <- AffymetrixCdfFile$byChipType("Mapping10K_Xba142");
print(cdf);
gi <- getGenomeInformation(cdf);
print(gi);
si <- getSnpInformation(cdf);
print(si);

csR <- AffymetrixCelSet$byName("GSE8605", cdf=cdf);
print(csR);

acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2");
print(acc);
csC <- process(acc, verbose=log);
print(csC);

bpn <- BasePositionNormalization(csC, target="zero");
print(bpn);
csN <- process(bpn, verbose=log);
print(csN);

plm <- NmfSnpPlm(csN, mergeStrands=TRUE);
print(plm);
fit(plm, ram=0.05, verbose=log);

ces <- getChipEffectSet(plm);
fln <- FragmentLengthNormalization(ces, target="zero");
print(fln);
cesN <- process(fln, verbose=log);
print(cesN);

ceR <- getAverage(cesN, verbose=log);

ce <- getFile(cesN, 1);

fig <- sprintf("%s", getFullName(ce));
if (!devIsOpen(fig)) {
  devSet(fig);
  subplots(3*2, ncol=2);
  par(mar=c(3,4,2,1)+0.1, pch=".");

  for (chr in getChromosomes(gi)[1:3]) {
    units <- getUnitsOnChromosome(gi, chr);
    pos <- getPositions(gi, units=units) / 1e6; 
    thetaR <- extractTotalAndFreqB(ceR, units=units)[,"total"];
    data <- extractTotalAndFreqB(ce, units=units);
    data[,"total"] <- 2*data[,"total"] / thetaR;  
  
  
    cn <- RawCopyNumbers(data[,"total"], pos, chromosome=chr);
    plot(cn, col="gray", cex=0.8, ylim=c(0,4));
    cnS <- gaussianSmoothing(cn, xOut=seq(xMin(cn), xMax(cn), by=1/2), sd=1);
    points(cnS, col="black");
    stext(side=3, pos=0, getName(ce));
    stext(side=3, pos=1, sprintf("Chr%d", chr));

    plot(pos, data[,"freqB"], cex=3, ylim=c(0,1));
    box(col="blue");
    stext(side=3, pos=0, getTags(cesN, collapse=","));
  } # for (chr ...)

  devDone();
}
