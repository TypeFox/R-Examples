#CalMaTe's Test Script
# - - - - - - - - - - - - - - - - - - - - - - -
# load the test data
#  \item{data}{An Jx2xI @numeric array, where J is the number of SNPs,
#          2 is the number of alleles, and I is the number of samples.}
#
# In this case, the variable data has the J=100 SNPs in the I=40 samples
# used in the manuscript
# - - - - - - - - - - - - - - - - - - - - - - -
library("calmate");

# Load example (thetaA,thetaB) signals
path <- system.file("exData", package="calmate");
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path);

# Scaling the initial data to CN scale
thetaR <- matrixStats::rowMedians(theta[,"A",] + theta[,"B",], na.rm=TRUE);
C <- 2*theta/thetaR;

# Transform to (total,fracB) signals
data <- thetaAB2TotalAndFracB(theta);

# It returns a list where the elements of the list are the SNPs and
# their total copy number and fracB
dataC <- calmateByTotalAndFracB(data);

# Allele-specific copy numbers after CalMaTe
CC <- C;
CC[,"A",] <- dataC[,"total",]*(1-dataC[,"fracB",]);
CC[,"B",] <- dataC[,"total",] - CC[,"A",];


if (interactive()) {
  devNew(type="x11", aspectRatio=1.9);
} else {
  devNew(type="png", "test-testScriptCalMaTe.png", aspectRatio=1.9);
}

subplots(2, ncol=1, byrow=FALSE);
par(mar=c(3,3,1,1)+0.1, mgp=c(1.8,0.7,0));

# Comparing allele-specific copy numbers before and after CalMaTe
# calibration for Sample #3
ii <- 3;
Clim <- c(-0.2,3);
plot(C[,,ii], xlim=Clim, ylim=Clim);
points(CC[,,ii], col="blue");

# Comparing BAFs before and after CalMaTe calibration for Sample #3.
ii <- 3;
plot(data[,"fracB",ii], ylab="BAF", ylim=c(0,1));
points(dataC[,"fracB",ii], col="blue");

devDone();
