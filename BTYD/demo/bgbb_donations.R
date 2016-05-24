## Authors: Lukasz Dziurzynski, Edward Wadsworth

data(donationsSummary)

## Get the calibration period recency-frequency matrix from the donation data:
rf.matrix <- donationsSummary$rf.matrix
rf.matrix

## Estimate parameters for the BG/BB model from the recency-frequency matrix:
par.start <- c(0.5, 1, 0.5, 1)
params <- bgbb.EstimateParameters(rf.matrix, par.start)
params

## Check log-likelihood of the params:
bgbb.rf.matrix.LL(params, rf.matrix)

## Plot the comparison of actual and expected calibration period frequencies:
bgbb.PlotFrequencyInCalibration(params, rf.matrix, censor=7, plotZero=TRUE)

n.star <- 5                        # Number of transaction opportunities in the holdout period
x.star <- donationsSummary$x.star  # Transactions made by each calibration period bin in the holdout period

## Plot the comparison of actual and conditional expected holdout period frequencies,
## binned according to calibration period frequencies:
bgbb.PlotFreqVsConditionalExpectedFrequency(params, n.star, rf.matrix, x.star)

## Plot the comparison of actual and conditional expected holdout period frequencies,
## binned according to calibration period recencies:
bgbb.PlotRecVsConditionalExpectedFrequency(params, n.star, rf.matrix, x.star)

inc.annual.trans <- donationsSummary$annual.trans           # incremental annual transactions
cum.annual.trans <- cumsum(donationsSummary$annual.trans)   # cumulative annual transactions
## set appropriate x-axis tickmarks:
x.tickmarks.yrs.all <- c( "'96","'97","'98","'99","'00","'01","'02","'03","'04","'05","'06" )

## Plot the comparison of actual and expected total cumulative transactions across
## both the calibration and holdout periods:
bgbb.PlotTrackingCum(params, rf.matrix, cum.annual.trans, xticklab=x.tickmarks.yrs.all)

## Plot the comparison of actual and expected total incremental transactions across
## both the calibration and holdout periods:
bgbb.PlotTrackingInc(params, rf.matrix, inc.annual.trans, xticklab=x.tickmarks.yrs.all)
