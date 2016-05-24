library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-50, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE13372,testset";
chipType <- "GenomeWideSNP_6,Full";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GC-content normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesN <- res$cesN;
gcn <- GcContentNormalization2(cesN, target="zero");
print(gcn);
cesN2 <- process(gcn, verbose=verbose);
print(cesN2);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GC-content effects before and after normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- getCdf(cesN);
ugp <- getAromaUgpFile(cdf);
units <- getUnitsOnChromosomes(ugp, 1:22);

# Down sample
units <- units[seq(from=1, to=length(units), by=4)];

# Before
toPNG(getFullName(gcn), tags=c("BeforeGCN"), {
  par(mar=c(4,3,1,1)+0.1, mgp=c(2.0,0.7,0));
  plotCovariateEffects(gcn, units=units, ylim=c(0,16), verbose=verbose);
});

# After
toPNG(getFullName(gcn), tags=c("AfterGCN"), {
  par(mar=c(4,3,1,1)+0.1, mgp=c(2.0,0.7,0));
  gcn2 <- GcContentNormalization2(cesN2, target="zero");
  plotCovariateEffects(gcn2, units=units, ylim=c(0,16), verbose=verbose);
});
