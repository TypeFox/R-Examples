library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, verbose=verbose);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extraction test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ces <- res$cesN;

# (a) Single array
ce <- ces[[1]];
Y <- extractMatrix(ce, verbose=verbose);
print(head(Y));

theta <- extractTheta(ce, verbose=verbose);
print(head(theta));

data <- extractDataFrame(ce, addNames=TRUE, verbose=verbose);
print(head(data));


# (b) Multiple arrays
Y <- extractMatrix(ces, verbose=verbose);
print(head(Y));

theta <- extractTheta(ces, verbose=verbose);
print(theta[1:6,,]);

data <- extractDataFrame(ces, addNames=TRUE, verbose=verbose);
print(head(data));

# Extract unit and group name and indices together with chip effects
dataT <- extractDataFrame(ces, units=1:50, addNames=TRUE);
print(head(dataT));
# Sanity check
stopifnot(all.equal(dataT[1:50,], data[1:50,]));

# SNPs only
unf <- getUnitNamesFile(ces);
print(unf);

units <- indexOf(unf, pattern="^SNP");
str(units);

theta <- extractTheta(ces, units=units[1:10]);
print(theta);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot theta:s
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
theta <- extractTheta(ces, units=units);
dimnames(theta)[[2]] <- c("A", "B");
ltheta <- log2(theta);
a <- (ltheta[,"A",]+ltheta[,"B",])/2;
m <- ltheta[,"A",]-ltheta[,"B",];

toPNG(getFullName(ces), tags=c("LogDiffvsLogSum"), aspectRatio=0.7, {
  Alim <- c(5,14);
  Mlim <- diff(Alim)*c(-1,1);
  Alab <- expression(1/2%*%log[2](theta[A]*theta[B]));
  Mlab <- expression(log[2](theta[A]/theta[B]));

  I <- dim(theta)[3];
  layout(matrix(1:I, ncol=3, byrow=TRUE));
  par(mar=c(3.8,4,3,1)+0.1);

  for (ii in 1:I) {
    name <- dimnames(theta)[[3]][ii];
    plot(a[,ii], m[,ii], pch=".", xlim=Alim, ylim=Mlim, xlab=Alab, ylab=Mlab, main=name);
  }
});
