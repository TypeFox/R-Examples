### Generate an object and store the parameters
reliabilityExample <- list();
reliabilityExample$nrOfItems <- 10;
reliabilityExample$trueScoreMean <- 20;
reliabilityExample$trueScoreVariance <- 4;
reliabilityExample$nrOfParticipants <- 250;
reliabilityExample$nrOfMoments <- 2;
reliabilityExample$itemErrorVariance <- 9;
reliabilityExample$transientErrorVariance <- 3;
reliabilityExample$itemVariancesMultiplierBasis <- .5;
reliabilityExample$itemMeanAdditionBasis <- 10;
reliabilityExample$seed <- 19811026;

### Set random seed
set.seed(reliabilityExample$seed);

### Determine itemvariance multiplyer for each item; get a random number
### from the normal distribution, then use exp to convert it to a sensible
### multiplier.
reliabilityExample$itemVarianceMultipliers <- 
  exp(rnorm(reliabilityExample$nrOfItems, sd=reliabilityExample$itemVariancesMultiplierBasis));

### Determine mean adjustment per item
reliabilityExample$itemMeanDeviations <- 
  rnorm(reliabilityExample$nrOfItems, sd=reliabilityExample$itemMeanAdditionBasis);

### Generate the dataframe and each participants' true score
reliabilityExample$dat <- data.frame(trueScore = 
                                       rnorm(reliabilityExample$nrOfParticipants,
                                             mean=reliabilityExample$trueScoreMean,
                                             sd=sqrt(reliabilityExample$trueScoreVariance)));

### Generate items nested within test administrations
for (currentMoment in 0:(reliabilityExample$nrOfMoments-1)) {
  ### Determine transient error for this moment for each participant
  currentTransientErrors <- rnorm(reliabilityExample$nrOfParticipants,
                                  sd=sqrt(reliabilityExample$transientErrorVariance));
  for (currentItem in 1:reliabilityExample$nrOfItems) {
    ### Determine item scores for this item for this moment for each participant;
    ### conform page 93 of Green (2003, "Test-retest alpha"), this is the participants'
    ###   true score plus
    ###   a random measurement error with a given variance plus
    ###   a transient error for this moment for this participant
    reliabilityExample$dat[[paste0("t", currentMoment, "_item", currentItem)]] <-
      reliabilityExample$itemMeanDeviations[currentItem] +
      reliabilityExample$itemVarianceMultipliers[currentItem] *
      (reliabilityExample$dat$trueScore +
         rnorm(reliabilityExample$nrOfParticipants, sd = sqrt(reliabilityExample$itemErrorVariance)) +
         currentTransientErrors);
  }
}

### Remove first column (with true scores) and store
### dataframe in the somewhat easier-to-use variable 'dat'
dat <- reliabilityExample$dat[2:(1+reliabilityExample$nrOfItems)];

### Show summaries (make sure to use the 'describe' function from
### the 'psych' package; the 'Hmisc' package also has one)
round(psych::describe(dat), 2);

### Order correlation matrix
round(cor(dat), 2);

### Extract administration at time 1
dat.time1 <- dat[1:reliabilityExample$nrOfItems];

### Check variables again to verify extraction
names(dat.time1);

### Order reliabilities for the administration at time1
### (remove 'ci=FALSE' or set it to true to compute confidence
###  intervals, as well)
scaleReliability(dat.time1, ci=FALSE);

### Order more comprehensive scale diagnostics
scaleDiags <- scaleDiagnosis(dat.time1);

### Show correlation matrices within each measurement moment
for (currentMoment in 0:(reliabilityExample$nrOfMoments-1)) {
  print(round(cor(reliabilityExample$dat[, paste0("t", currentMoment, "_item", 1:reliabilityExample$nrOfItems)]), digits=2));
}

### Show complete correlation matrix
print(cor(reliabilityExample$dat[, 2:ncol(reliabilityExample$dat)]), digits=2);

### And covariance matrix
cov(reliabilityExample$dat[, 2:ncol(reliabilityExample$dat)]);

### Compute reliability estimates for each measurement moment
for (currentMoment in 0:(reliabilityExample$nrOfMoments-1)) {
  print(scaleReliability(reliabilityExample$dat,
                         paste0("t", currentMoment, "_item", 1:reliabilityExample$nrOfItems),
                         ci=FALSE));
}

### Split dataframe into separate dataframes for each measurement moment
reliabilityExample$dat.split <- list();
for (currentMoment in 0:(reliabilityExample$nrOfMoments-1)) {
  reliabilityExample$dat.split[[paste0("t", currentMoment)]] <-
    reliabilityExample$dat[, paste0("t", currentMoment, "_item", 1:reliabilityExample$nrOfItems)];
}

### Show item-time covariances for the first two measurements
print(cov(x=reliabilityExample$dat.split$t0,
          y=reliabilityExample$dat.split$t1), digits=2);

### Show test-retest alpha
print(testRetestAlpha(reliabilityExample$dat[, 2:ncol(reliabilityExample$dat)]));

### Show test-retest CES
print(testRetestCES(reliabilityExample$dat[, 2:ncol(reliabilityExample$dat)]));

### Show test-retest CES when subscales have uneven numbers of items
start_time1 <- 2;
end_time1 <- start_time1 + reliabilityExample$nrOfItems - 1;
start_time2 <- end_time1 + 1;
end_time2 <- start_time2 + reliabilityExample$nrOfItems - 1;
itemSelection <- c(start_time1:(end_time1-1), start_time2:(end_time2-1));

print(testRetestCES(reliabilityExample$dat[, itemSelection]));

### Compute both at the same time
print(testRetestReliability(reliabilityExample$dat[, itemSelection]));

### Or for the complete scales
print(testRetestReliability(reliabilityExample$dat[, 2:ncol(reliabilityExample$dat)]));

### Once more with less items
itemSelection <- c(start_time1:(end_time1-5), start_time2:(end_time2-5));
print(testRetestReliability(reliabilityExample$dat[, itemSelection]));

### For comparison, the single administration reliability measures:
itemSelection <- c(start_time1:(end_time1-5));
print(scaleReliability(reliabilityExample$dat[, itemSelection], ci=FALSE));

### And for time 2:
itemSelection <- c(start_time2:(end_time2-5));
print(scaleReliability(reliabilityExample$dat[, itemSelection], ci=FALSE));

### Uncomment these line to store the data, respectively as .csv or as an R dataframe object.
### Don't forget to change the paths!

# write.csv(reliabilityExample$dat, row.names = FALSE,
#           file="B:/Data/statistics/R/library/userfriendlyscience/data/testRetestSimData.csv");
# testRetestSimData <- reliabilityExample$dat;
# save(testRetestSimData,
#     file="B:/Data/statistics/R/library/userfriendlyscience/data/testRetestSimData.rda")

### To save a version without the first column (i.e. without the true score), so that it
### can be loaded directly into the functions when providing no arguments, use:

# write.csv(reliabilityExample$dat[, 2:ncol(reliabilityExample$dat)], row.names = FALSE,
#           file="B:/Data/research/cronbach's alpha - reliability and validity/osf/exampleData.csv");

### Or, to save less items, use itemSelection:

# write.csv(reliabilityExample$dat[, itemSelection], row.names = FALSE,
#           file="B:/Data/research/cronbach's alpha - reliability and validity/osf/selectedExampleData.csv");