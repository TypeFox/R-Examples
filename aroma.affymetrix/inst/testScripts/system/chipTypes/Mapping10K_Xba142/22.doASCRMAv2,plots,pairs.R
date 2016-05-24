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
# (a) AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);

cesN <- res$cesN;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extraction test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract unit and group name and indices together with chip effects
data <- extractDataFrame(cesN, units=1:50, addNames=TRUE, verbose=verbose);
print(data[1:50,1:8]);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract (thetaA, thetaB)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
theta <- extractTheta(cesN, groups=1:2, verbose=verbose);
dimnames(theta)[[2]] <- c("A", "B");
thetaA <- theta[,"A",];
thetaB <- theta[,"B",];
theta <- thetaA + thetaB;
freqB <- thetaB / theta;


tlim <- c(8,15);
Blim <- c(0,1);
tlab <- expression(log[2](theta));
Blab <- expression(theta[B]/theta);

nbrOfArrays <- ncol(theta);

toPNG(getFullName(cesN), tags="BAFvsTheta", width=1024, {
  layout(matrix(1:nbrOfArrays, ncol=2, byrow=TRUE));
  par(mar=c(3.8,4,3,1)+0.1);
  for (ii in seq_len(nbrOfArrays)) {
    name <- colnames(theta)[ii];
    smoothScatter(log2(theta[,ii]), freqB[,ii], 
                  xlim=tlim, ylim=Blim, xlab=tlab, ylab=Blab, main=name);
  }
});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract (log2(total), freqB)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data <- extractTotalAndFreqB(cesN, verbose=verbose);
data[,1,] <- log2(data[,1,]);

tlim <- c(9,15);

upperPanel <- function(x,y, pch=".", ...) {
  xy <- data[,1,c(x,y)];
  xy <- (xy - tlim[1])/diff(tlim);
  abline(a=0,b=1, col="#999999");
  points(xy, pch=pch);
}

lowerPanel <- function(x,y, pch=".", ...) {
  xy <- data[,2,c(x,y)];
  abline(a=0,b=1, col="#999999");
  abline(a=1,b=-1, col="#999999");
  points(xy, pch=pch);
}

diagPanel <- function(x, pch=".", ...) {
  xy <- data[,,x];
  xy[,1] <- (xy[,1] - tlim[1])/diff(tlim);
  smoothScatter(xy, add=TRUE);
  abline(h=1/2, col="#999999");
}

pairs <- matrix(1:nbrOfArrays, nrow=1, ncol=nbrOfArrays);
colnames(pairs) <- dimnames(data)[[length(dim)]];

toPNG(getFullName(cesN), tags="pairs", width=1024, {
  pairs(pairs, upper.panel=upperPanel, lower.panel=lowerPanel, 
               diag.panel=diagPanel, xlim=c(0,1), ylim=c(0,1));
});
