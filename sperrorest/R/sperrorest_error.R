
#' Default error function
#'
#' Calculate a variety of accuracy measures from observations and predictions of numerical and categorical response variables.
#'
#' @name err.default
#' @param obs factor, logical, or numeric vector with observations
#' @param pred factor, logical, or numeric vector with predictions. Must be of same type as \code{obs} with the exception that \code{pred} may be numeric if \code{obs} is \code{factor} or \code{logical} ('soft' classification).
#' @return A list with (currently) the following components, depending on the type of prediction problem:
#' \item{'hard' classification}{misclassification error, overall accuracy; if two classes, sensitivity, specificity, positive predictive value (PPV), negative predictive value (NPV), kappa}
#' \item{'soft' classification}{area under the ROC curve, error and accuracy at a obs>0.5 dichotomization, false-positive rate (FPR; 1-specificity) at 70, 80 and 90 percent sensitivity, true-positive rate (sensitivity) at 80, 90 and 95 percent specificity}
#' \item{regression}{bias, standard deviation, mean squared error, MAD (\code{\link{mad}}), median, interquartile range (\code{\link{IQR}}) of residuals}
#' @note \code{NA} values are currently not handled by this function, i.e. they will result in an error.
#' @seealso \pkg{ROCR}
#' @examples
#' obs = rnorm(1000)
#' # Two mock (soft) classification examples:
#' err.default( obs>0, rnorm(1000) ) # just noise
#' err.default( obs>0, obs + rnorm(1000) ) # some discrimination
#' # Three mock regression examples:
#' err.default( obs, rnorm(1000) ) # just noise, but no bias
#' err.default( obs, obs + rnorm(1000) ) # some association, no bias
#' err.default( obs, obs + 1 ) # perfect correlation, but with bias
#' @export
err.default = function(obs, pred) {
    # The following wrapper functions are used in order to avoid warning messages:
    mmin = function(x,...) {
        if (length(x) == 0) x = Inf
        min(x,...)
    }
    mmax = function(x,...) {
        if (length(x) == 0) x = -Inf
        max(x,...)
    }

    # Convert logical to factor if necessary:
    if (is.logical(obs)) obs = factor(obs, levels = c("FALSE","TRUE"))
    if (is.logical(pred)) pred = factor(pred, levels = c("FALSE","TRUE"))

    # Classification problem:
    if (is.factor(obs)) {
        if (is.factor(pred)) { # "hard" classification:
            pred = as.character(pred)
            pred = factor(pred, levels = levels(obs))
            err = list(
                    error = mean(obs != pred),
                    accuracy = mean(obs == pred) )
            if (nlevels(obs) == 2) {
                npos = sum( obs == levels(obs)[2] )
                nneg = sum( obs == levels(obs)[1] )
                ntruepos = sum( (obs == levels(obs)[2]) & (pred == levels(obs)[2]) )
                ntrueneg = sum( (obs == levels(obs)[1]) & (pred == levels(obs)[1]) )
                err$sensitivity = ntruepos / npos
                err$specificity = ntrueneg / nneg
                npospred = sum( pred == levels(obs)[2] )
                nnegpred = sum( pred == levels(obs)[1] )
                err$ppv = ntruepos / npospred
                err$npv = ntrueneg / nnegpred
                n = length(obs)
                pexp = (npos/n)*(npospred/n) + (nneg/n)*(nnegpred/n)
                if (pexp == 1) { err$kappa = NA
                } else err$kappa = (err$accuracy - pexp) / (1 - pexp)
            }
        } else { # "soft" classification:
            # Calculate area under the ROC curve:
            require(ROCR)
            predobj = prediction( pred, obs )
            auroc = performance( predobj, measure="auc" )@y.values[[1]]
            err = list( auroc = auroc )
            
            pos = levels(obs)[2]
            neg = levels(obs)[1]

            err$error = mean( (obs==pos) != (pred>=0.5) )
            err$accuracy = 1 - err$error

            # internal functions for calculating false positive rate (1-specificity)
            # and true positive rate (sensitivity) for a given decision threshold:
            tpr = function(o,p,np,t)  sum(o & (p >= t)) / np
            fpr = function(o,p,nn,t)  sum(!o & (p >= t)) / nn

            # Number of positive and negative predictions:
            npos = sum(obs == pos)
            nneg = sum(obs != pos)
            
            err$sensitivity = tpr(obs==pos, pred, npos, 0.5)
            err$specificity = 1 - fpr(obs==pos, pred, nneg, 0.5)

            thrs = unique(pred)
            if (length(thrs) > 500) thrs = seq(min(pred)+0.0001,max(pred)+0.0001,length=500)

            thr = mmax( thrs[ sapply( thrs, function(x) tpr(obs==pos,pred,npos,x) >= 0.70 ) ] )
            err$fpr70 = fpr(obs==pos, pred, nneg, thr)
            thr = mmax( thrs[ sapply( thrs, function(x) tpr(obs==pos,pred,npos,x) >= 0.80 ) ] )
            err$fpr80 = fpr(obs==pos, pred, nneg, thr)
            thr = mmax( thrs[ sapply( thrs, function(x) tpr(obs==pos,pred,npos,x) >= 0.90 ) ] )
            err$fpr90 = fpr(obs==pos, pred, nneg, thr)

            thr = mmin( thrs[ sapply( thrs, function(x) fpr(obs==pos,pred,nneg,x) <= 0.20 ) ] )
            err$tpr80 = tpr(obs==pos, pred, npos, thr)
            thr = mmin( thrs[ sapply( thrs, function(x) fpr(obs==pos,pred,nneg,x) <= 0.10 ) ] )
            err$tpr90 = tpr(obs==pos, pred, npos, thr)
            thr = mmin( thrs[ sapply( thrs, function(x) fpr(obs==pos,pred,nneg,x) <= 0.05 ) ] )
            err$tpr95 = tpr(obs==pos, pred, npos, thr)
        }
        err$events = sum( obs == levels(obs)[2] )
    } else {
    # Regression problem:
        err = list(
            bias = mean( obs - pred ),
            stddev = sd( obs - pred ),
            mse = mean( (obs - pred)^2 ),
            mad  = mad( obs - pred ),
            median = median( obs - pred ),
            iqr = IQR( obs - pred ) )
    }
    # Number of observations available:
    err$count = length(obs)

    return(err)
}
