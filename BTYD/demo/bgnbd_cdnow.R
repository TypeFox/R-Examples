## Authors: Daniel McCarthy, Lukasz Dziurzynski, Edward Wadsworth

data(cdnowSummary)

## Get the calibration period customer-by-sufficient-statistic matrix from the cdnow data:
cbs <- cdnowSummary$cbs

## Estimate parameters for the BG/NBD model from the CBS:
par.start <- c(1,3,1,3)
params <- bgnbd.EstimateParameters(cbs, par.start)
params


## Check log-likelihood of the params:
bgnbd.cbs.LL(params, cbs)

## Plot the comparison of actual and expected calibration period frequencies:
bgnbd.PlotFrequencyInCalibration(params, cbs, censor=7, plotZero=TRUE)

T.star <- 39                # Length of holdout period
x.star <- cbs[,"x.star"]    # Transactions made by each customer in the holdout period

## Plot the comparison of actual and conditional expected holdout period frequencies,
## binned according to calibration period frequencies:
bgnbd.PlotFreqVsConditionalExpectedFrequency(params, T.star, cbs, x.star, censor=7)

## Plot the comparison of actual and conditional expected holdout period frequencies,
## binned according to calibration period recencies:
bgnbd.PlotRecVsConditionalExpectedFrequency(params, cbs, T.star, x.star)

