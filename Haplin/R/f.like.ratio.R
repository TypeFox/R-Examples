f.like.ratio <- function(res.0, res, data, info){
##
## PERFORMS A LIKELIHOOD RATIO TEST BETWEEN THE TWO RESULTS IN res.0 AND res.
## res.0 AND res ARE OF TYPE tri.glm, pred.0 and pred ARE res0$pred and res$pred.
## IMPORTANT: THE LIKELIHOODS ARE COMPUTED ONLY UP TO A CONSTANT. SO COMPARISONS
## SHOULD BE MADE ONLY ON NESTED MODELS COMPUTED FROM THE SAME (SINGLE) DATA SET!
##
## THE TEST CORRECTS FOR EM UNCERTAINTY. SEE COMMENTS IN f.final.loglike
##
#
## COMPUTE df's (THIS SHOLD BE DONE IN A MORE DIRECT FASHION, BUT OK....)
.anova <- anova(res.0$result, res$result, test = "Chisq")
.df <- .anova$Df[2]
## COMPUTE FINAL LIKELIHOOD FOR TWO RESULTS, TAKING EM INTO ACCOUNT
.loglike.0 <- f.final.loglike(data = data, pred = res.0$pred, info = info, type = "EM")
.loglike <- f.final.loglike(data = data, pred = res$pred, info = info, type = "EM")
#
## LOG RATIO OF THE TWO
.loglike.ratio <- 2*(.loglike - .loglike.0)
names(.loglike.ratio) <- paste(names(.loglike.ratio), ".ratio", sep = "")
#
## COMPUTE TEST P-VALUE
.lratio.test <- pchisq(.loglike.ratio, df = .df, lower.tail = F)
#
## COLLECT OUTPUT VALUES
.lratio.ut <- c(loglike.0 = unname(.loglike.0), loglike = unname(.loglike), df = .df, p.value.overall = unname(.lratio.test))
#
##
return(.lratio.ut)
}
