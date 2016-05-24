library("ACNE");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "HapMap270,100K,CEU,5trios";
chipType <- "Mapping50K_Hind240";
res <- doACNE(dataSet, chipType=chipType, verbose=verbose);
print(res);

ds <- res$total;
dfR <- getAverageFile(ds, verbose=verbose);
df <- getFile(ds, 1);
baf <- getFile(res$fracB, 1);
ugp <- getAromaUgpFile(ds);

fig <- sprintf("%s", getFullName(df));
if (!devIsOpen(fig)) {
  devSet(fig, width=10, height=5);

  subplots(2*3, nrow=2, byrow=FALSE);
  par(mar=c(3,4,2,1)+0.1, pch=".");

  for (chr in 1:3) {
    units <- getUnitsOnChromosome(ugp, chr);
    pos <- getPositions(ugp, units=units);

    beta <- extractMatrix(baf, units=units, drop=TRUE);
    fracB <- RawAlleleBFractions(beta, pos, chromosome=chr);

    theta <- extractMatrix(df, units=units, drop=TRUE);
    thetaR <- extractMatrix(dfR, units=units, drop=TRUE);
    C <- 2 * theta/thetaR;
    cn <- RawCopyNumbers(C, pos, chromosome=chr);

    plot(cn, col="gray", cex=0.8, ylim=c(0,4));
    xOut <- seq(xMin(cn), xMax(cn), by=0.5e6);
    cnS <- gaussianSmoothing(cn, xOut=xOut, sd=1e6);
    points(cnS, col="black");
    stext(side=3, pos=0, getName(df));
    stext(side=3, pos=1, sprintf("Chr%d", chr));

    plot(fracB, ylim=c(0,1));
    box(col="blue");
    stext(side=3, pos=0, getTags(ds, collapse=","));
    stext(side=3, pos=1, sprintf("Chr%d", chr));
  } # for (chr ...)

  devDone();
}
