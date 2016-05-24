#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
##        
## Time-stamp: <[plotCurveEstimate.R] by DSB Mit 26/01/2011 14:19 (CET)>
##
## Description:
## Plot predictor curve estimates. 
##
## History:
## 26/01/2011   - add this header
##              - add options "partialResids" and "hpd"
#####################################################################################

`plotCurveEstimate` <-
function (model, termName, plevel = 0.95, slevel = plevel, plot = TRUE, legendPos = "topleft",
          rug=FALSE, partialResids=TRUE, hpd=TRUE, ..., main = NULL)
    UseMethod ("plotCurveEstimate")     # define generic function: BayesMfp and BmaSamples methods following

