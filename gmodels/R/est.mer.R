# est.mer.R
# generate estimable output for mer objects using mcmcsamp()
# Randall Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created April 25, 2006
# Updated 2012-04-19 for S4 version of lmer object

## est.mer <- function(obj, cm, beta0, conf.int, show.beta0, n.sim)
## {

##   samp <- lme4:::mcmcsamp(obj, n.sim)
##   ##  samp.summ <- summary(samp)

##   samp.cm <- t(cm %*% samp@fixef)

##   # calculate requested statistics
##   est <- apply(samp.cm, 2, mean)
##   stderr <- apply(samp.cm, 2, sd)

##   pval <- sapply(1:length(beta0),
##                  function(i){percentile(beta0[i], samp.cm[,i])})
##   pval <- ifelse(pval <= .5, 2*pval, 2*(1-pval))

##   if(is.null(conf.int))
##     {
##       lower.ci <- NULL
##       upper.ci <- NULL
##     }
##   else
##     {
##       alpha <- 1-conf.int
##       samp.ci <- sapply(1:length(beta0),
##                         function(i)
##                         {
##                           quantile(samp.cm[,i], probs=c(alpha/2, 1-alpha/2))
##                         }
##                         )

##       lower.ci <- samp.ci[1,]
##       upper.ci <- samp.ci[2,]
##     }

##   # return results
##   if(!show.beta0)
##     beta0 <- NULL

##   samp.stats <- cbind('beta0' = beta0,
##                       'Estimate' = est,
##                       'Std. Error' = stderr,
##                       'p value' = pval,
##                       'Lower.CI' = lower.ci,
##                       'Upper.CI' = upper.ci)

##   row.names(samp.stats) <- paste('(', apply(cm, 1, paste, collapse=" "),
##                                  ')', sep='')

##   return(samp.stats)
## }

