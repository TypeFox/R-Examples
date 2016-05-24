#' A common problem faced by journal reviewers and authors is the question of
#' whether the results of a replication study are consistent with the original
#' published study. One solution to this problem is to examine the effect size
#' from the original study and generate the range of effect sizes that could
#' reasonably be obtained (due to random sampling) in a replication attempt
#' (i.e., calculate a replication interval). If a replication effect size falls
#' outside the replication interval then that value could not have occurred to
#' due the effects of sampling error alone. Alternatively, if a replication
#' effect size falls within the replication interval then the replication
#' results could have reasonably occurred due to the effects of sampling error
#' alone. This package calculates the replication interval for two different types of effect sizes (i.e., correlation \emph{r}, standardized mean difference \emph{d} ).
#'\tabular{ll}{
#'Package: \tab replicationInterval\cr
#'Type: \tab Package\cr
#'Version: \tab 1.0.0\cr
#'Date: \tab 2015-07-07\cr
#'License: \tab Unlimited\cr
#'}
#'\code{\link{ri.d}} creates a replication interval for a standardized mean difference (i.e., \emph{d}-value) \cr
#'\code{\link{ri.r}} creates a replication interval for a correlation (i.e., \emph{r} )\cr
#'
#'@name replicationInterval-package
#'@aliases replicationInterval
#'@docType package
#'@title Replication Interval Functions
#'@author 
#'\tabular{ll}{
#'Author: \tab David J. Stanley \email{dstanley@@uoguelph.ca}\cr
#'Maintainer: \tab David J. Stanley \email{dstanley@@uoguelph.ca}
#'}
#'@references 
#'Stanley, D.J., & Spence, J.R.(2014). Expectations for replications: Are yours realistic? \emph{Perspectives on Psychological Science, 9}, 305-318.\cr\cr
#'Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and calculation of confidence intervals that are based on central and noncentral distributions. \emph{Educational and Psychological Measurement, 61(4)}, 532-574.
#'
#'@keywords package
#'@examples
#' ri.d(d=.65,n1=20,n2=20)
#' ri.d(d=.65,n1=20,n2=20,explain=TRUE)
#' ri.d(d=.65,n1=20,n2=20,rep.n1=40,rep.n2=40,explain=TRUE)


#' ri.r(r=.35,n=40)
#' ri.r(r=.35,n=40,explain=TRUE)
#' ri.r(r=.35,n=40,rep.n=80,explain=TRUE)

NULL



quantile.d <- function (probability.in.tail,ncd,n1,n2,is.lower.tail) {
     nct.value=ncd/sqrt((1/n1)+(1/n2))
     df.t=n1+n2-2
     t.tail=stats::qt(p=probability.in.tail, ncp=nct.value,df=df.t,lower.tail=is.lower.tail)
     d.tail=t.tail*sqrt((1/n1)+(1/n2))
     return(d.tail)
}


ri.d.output <- function(interval.output) {
     conf.level.percent=round(interval.output$conf.level*100,0)
     
     d <- interval.output$d
     lower.conf.limit <- interval.output$lower.confidence.interval.d 
     upper.conf.limit <- interval.output$upper.confidence.interval.d
     n1 <- interval.output$confidence.interval.n1
     n2 <- interval.output$confidence.interval.n2
     
     lower.population.lower.bound <- interval.output$lower.population.lower.bound.d
     lower.population.upper.bound <- interval.output$lower.population.upper.bound.d 
     
     upper.population.lower.bound <- interval.output$upper.population.lower.bound.d
     upper.population.upper.bound <- interval.output$upper.population.upper.bound.d 
     
     lower.replication.interval <- interval.output$lower.replication.interval.d
     upper.replication.interval <- interval.output$upper.replication.interval.d
     
     n1.replication <- interval.output$replication.interval.n1
     n2.replication <- interval.output$replication.interval.n2
     
#      text.ci <- sprintf("Due to sampling error, the published/sample d = %1.2f (N1 = %d and N2 = %d) could have been created by a population d-value as low as %1.2f (lower-population d-value) or as high as %1.2f (upper-population d-value).",d,n1,n2,lower.conf.limit,upper.conf.limit)
     
     text.ci <- sprintf("The %d%% confidence interval for the original sample d = %1.2f (N1 = %d and N2 = %d) is [%1.2f, %1.2f]. This confidence interval indicates a plausible range of population d-values that could have created the original sample d-value. Thus, due to the effects of random sampling error, the original sample d = %1.2f (N1 = %d and N2 = %d) could have been created by a population d-value as low as %1.2f (lower-population d-value) or as high as %1.2f (upper-population d-value).",conf.level.percent, d,n1,n2,lower.conf.limit,upper.conf.limit,d,n1,n2,lower.conf.limit,upper.conf.limit)
     
     text.lower.population <- sprintf("If the lower-population d-value of %1.2f created the original d = %1.2f (due to sampling error) then %d%% of replication d-values (using N1 = %d and N2 = % d) will fall between %1.2f and %1.2f.",lower.conf.limit,d,conf.level.percent,n1.replication,n2.replication,lower.population.lower.bound,lower.population.upper.bound)
     
     
     text.upper.population <- sprintf("If the upper-population d-value of %1.2f created the original d = %1.2f (due to sampling error) then %d%% of replication d-values (using N1 = %d and N2 = % d) will fall between %1.2f and %1.2f.",upper.conf.limit,d,conf.level.percent,n1.replication,n2.replication,upper.population.lower.bound,upper.population.upper.bound)
     
     
     text.summary <- sprintf("That is, if the study was replicated (using N1 = %d and N2 = %d) the researcher could obtain any d-value in the %1.2f to %1.2f interval due to sampling error.",n1.replication,n2.replication,lower.replication.interval,upper.replication.interval)
     
     text.replication.interval.brief <- sprintf("Replication Interval for d = %1.2f is {%1.2f, %1.2f; N1 = %d, N2 = %d} based on the %d%% Confidence Interval [%1.2f, %1.2f].", d,lower.replication.interval,upper.replication.interval,n1.replication,n2.replication,conf.level.percent,lower.conf.limit,upper.conf.limit)
     
     text.replication.interval <- sprintf("Thus, for d = %1.2f (when replications are conducted with N1 = %d and N2 = %d) the replication interval is {%1.2f, %1.2f}.", d,n1.replication,n2.replication,lower.replication.interval,upper.replication.interval)
     
     
     text.explain <- paste(text.ci,text.lower.population,text.upper.population, text.summary,text.replication.interval)
     
     text.method.section <- sprintf("If our replication study obtains an effect between d = %1.2f and d = %1.2f this will be evidence that our findings are not inconsistent with the original study effect of d = %1.2f, %d%% CI [%1.2f,%1.2f]{%1.2f, %1.2f; N1 = %d, N2 = %d}.",lower.replication.interval,upper.replication.interval,d,conf.level.percent,lower.conf.limit,upper.conf.limit,lower.replication.interval, upper.replication.interval,n1.replication,n2.replication)
     
     
     

     
     text.out <- list()
     text.out$replication.interval <- text.replication.interval.brief
     text.out$replication.interval.explanation <- text.explain
     text.out$method.section.text <- text.method.section
     
     return(text.out)
}


#' Creates a replication interval based for a \emph{d}-value (i.e., standardized mean difference) 
#' @param d Original study: Sample \emph{d}-value (standardized mean difference) created with pooled variance denominator. See formulas 4.18 and 4.19 (p.26) in Borenstein, Hedges, Higgins, & Rothstein (2009).
#' @param n1 Original study: Sample size for group 1
#' @param n2 Original study: Sample size for group 2
#' @param rep.n1 (optional) Replication study: Sample size for group 1. If not specified, n1 is used.
#' @param rep.n2 (optional) Replication study: Sample size for group 2. If not specified, n2 is used.
#' @param conf.level (optional 0 to 1 value) Confidence level desired (0 to 1). If not specified .95 (i.e., 95 percent) will be used.
#' @param explain (optional boolean) Default is FALSE. If TRUE, text output explaining the interval is provided. 
#' @param extended.output (optional boolean) Default is FALSE. If TRUE, additional details (e.g., confidence interval) provided in numeric return output.
#' @param manuscript.text (optional boolean) Default is TRUE. If TRUE, present text for Method section.
#' @return A list of values (\code{lower.replication.interval.d, upper.replication.interval.d}) containing the replication interval (and related statistics if requested with the \code{extended.output} argument).
#' @references
#' Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2009). \emph{Introduction to meta-analysis}. John Wiley & Sons.\cr\cr
#'Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and calculation of confidence intervals that are based on central and noncentral distributions. \emph{Educational and Psychological Measurement, 61(4)}, 532-574.
#' @examples
#' ri.d(d=.65,n1=20,n2=20)
#' ri.d(d=.65,n1=20,n2=20,explain=TRUE)
#' @export 
ri.d <- function (d,n1,n2,rep.n1=NA, rep.n2=NA,conf.level=.95,explain=FALSE,extended.output=FALSE,manuscript.text = TRUE) {
     
      if (is.na(rep.n1)) {
        rep.n1=n1
      }
      if (is.na(rep.n2)) {
        rep.n2=n2
      }
  
     d.value.ci=MBESS::ci.smd(smd=d,n.1=n1,n.2=n2,conf.level=conf.level)
     d <- d.value.ci$smd
     
     conf.lower.bound <- d.value.ci$Lower.Conf.Limit.smd
     conf.upper.bound <- d.value.ci$Upper.Conf.Limit.smd
     probability.in.tail=(1 - conf.level)/2
          
     lower.population.lower.bound <- quantile.d(probability.in.tail=probability.in.tail, ncd=conf.lower.bound,n1=rep.n1,n2=rep.n2,is.lower.tail=TRUE)
     lower.population.upper.bound <- quantile.d(probability.in.tail=probability.in.tail, ncd=conf.lower.bound,n1=rep.n1,n2=rep.n2,is.lower.tail=FALSE)
     upper.population.lower.bound <- quantile.d(probability.in.tail=probability.in.tail, ncd=conf.upper.bound,n1=rep.n1,n2=rep.n2,is.lower.tail=TRUE)
     upper.population.upper.bound <- quantile.d(probability.in.tail=probability.in.tail, ncd=conf.upper.bound,n1=rep.n1,n2=rep.n2,is.lower.tail=FALSE)
     
     interval.output <- list()
     interval.output$conf.level <- conf.level
     interval.output$d <- d
     interval.output$lower.confidence.interval.d <- conf.lower.bound
     interval.output$upper.confidence.interval.d <- conf.upper.bound
     interval.output$confidence.interval.n1 <- n1
     interval.output$confidence.interval.n2 <- n2
     interval.output$lower.population.lower.bound.d <- lower.population.lower.bound
     interval.output$lower.population.upper.bound.d <- lower.population.upper.bound
     interval.output$upper.population.lower.bound.d <- upper.population.lower.bound
     interval.output$upper.population.upper.bound.d <- upper.population.upper.bound
     interval.output$replication.interval.n1 <- rep.n1
     interval.output$replication.interval.n2 <- rep.n2
     interval.output$lower.replication.interval.d <- lower.population.lower.bound
     interval.output$upper.replication.interval.d <- upper.population.upper.bound

     percent.ri.due.to.ci <- ((conf.upper.bound - conf.lower.bound) /(interval.output$upper.replication.interval.d - interval.output$lower.replication.interval.d))*100
     interval.output$percent.replication.interval.due.to.confidence.interval <- round(percent.ri.due.to.ci)
     
     
     
     
     interval.output.brief <- list()
     interval.output.brief$lower.confidence.interval.d <- conf.lower.bound
     interval.output.brief$upper.confidence.interval.d <- conf.upper.bound
     interval.output.brief$lower.replication.interval.d <- lower.population.lower.bound
     interval.output.brief$upper.replication.interval.d <- upper.population.upper.bound
     
     if (extended.output==TRUE) {
          interval.results <- interval.output
     } else {
          interval.results <- interval.output.brief
     }
     
     
     if (manuscript.text==TRUE | explain == TRUE) {
          rText<-ri.d.output(interval.output)
          cat("\n")
          cat(rText$replication.interval)
          cat("\n")
     }      
     
     if (manuscript.text == TRUE) {
          cat("\n")
          cat("Method Section Snippet:")
          cat("\n")
          cat(rText$method.section.text)
          cat("\n")
     }
     
     if (explain==TRUE) {
          cat("\nExplanation:\n")
          cat(rText$replication.interval.explanation)
          cat("\n\n")
     }
     
     if (manuscript.text==TRUE | explain == TRUE) {
          cat("\n\n")     
     }      
     
#      if (explain==TRUE) {
#           cat("\n")
#           rText<-ri.d.output(interval.output)
#           cat(rText[[1]])
#           cat("\n\n")
#           cat(rText[[2]])
#           cat("\n\n\n")
#           
#      }

     
     return(interval.results)
}


















quantile.r <- function (probability.in.tail,ncr,n,is.lower.tail) {
     ncz <- atanh(ncr)
     ncz.se <- 1 / sqrt(n-3)
     z.tail <- stats::qnorm(p=probability.in.tail, mean=ncz, sd=ncz.se,lower.tail=is.lower.tail)
     r.tail <- tanh(z.tail)
     return(r.tail)
}

ci.r <- function(r,n,conf.level) {
     probability.in.tail <- (1- conf.level)/2
     
     obs.z <- atanh(r)
     obs.z.se <- 1/sqrt(n-3)
     
     high.z <- stats::qnorm(p=probability.in.tail, mean=obs.z,sd=obs.z.se,lower.tail=FALSE)
     high.r <- tanh(high.z)
     
     low.z <- stats::qnorm(p=probability.in.tail, mean=obs.z,sd=obs.z.se,lower.tail=TRUE)
     low.r <- tanh(low.z)
     
     
     conf.interval.output <- list()
     
     conf.interval.output$r <- r
     conf.interval.output$lower.conf.limit.r <- low.r
     conf.interval.output$upper.conf.limit.r <- high.r
     return(conf.interval.output)
}

ri.r.output <- function(interval.output) {
     conf.level.percent=round(interval.output$conf.level*100,0)
     
     r <- interval.output$r
     lower.conf.limit <- interval.output$lower.confidence.interval.r 
     upper.conf.limit <- interval.output$upper.confidence.interval.r
     n <- interval.output$upper.confidence.interval.n
     
     lower.population.lower.bound <- interval.output$lower.population.lower.bound.r
     lower.population.upper.bound <- interval.output$lower.population.upper.bound.r 
     
     upper.population.lower.bound <- interval.output$upper.population.lower.bound.r
     upper.population.upper.bound <- interval.output$upper.population.upper.bound.r 
     
     lower.replication.interval <- interval.output$lower.replication.interval.r
     upper.replication.interval <- interval.output$upper.replication.interval.r
     
     n.replication <- interval.output$replication.interval.n
     
#      text.ci <- sprintf("Due to sampling error, the published/sample r = %1.2f (N = %d) could have been created by a population correlation as low as %1.2f (lower-population correlation) or as high as %1.2f (upper-population correlation).",r,n,lower.conf.limit,upper.conf.limit)
     
     text.ci <- sprintf("The %d%% confidence interval for the original sample r = %1.2f (N = %d) is [%1.2f, %1.2f]. This confidence interval indicates a plausible range of population correlations that could have created the original sample correlation. Thus, due to the effects of random sampling error, the original sample r = %1.2f (N = %d) could have been created by a population correlation as low as %1.2f (lower-population correlation) or as high as %1.2f (upper-population correlation).",conf.level.percent,r,n,lower.conf.limit,upper.conf.limit,r,n,lower.conf.limit,upper.conf.limit)     
     
     
     
     text.lower.population <- sprintf("If the lower-population correlation of %1.2f created the original r = %1.2f (due to sampling error) then %d%% of replication correlations (using N = %d) will fall between %1.2f and %1.2f.",lower.conf.limit,r,conf.level.percent,n.replication,lower.population.lower.bound,lower.population.upper.bound)
     
     
     text.upper.population <- sprintf("If the upper-population correlation of %1.2f created the original r = %1.2f (due to sampling error) then %d%% of replication correlations (using N = %d) will fall between %1.2f and %1.2f.",upper.conf.limit,r,conf.level.percent,n.replication,upper.population.lower.bound,upper.population.upper.bound)
     
     
     text.summary <- sprintf("That is, if the study was replicated (using N = %d) the researcher could obtain any correlation in the %1.2f to %1.2f interval due to sampling error.",n.replication,lower.replication.interval,upper.replication.interval)
     
     text.replication.interval.brief <- sprintf("Replication Interval for r = %1.2f is {%1.2f, %1.2f; N = %d} based on the %d%% Confidence Interval [%1.2f, %1.2f].", r,lower.replication.interval,upper.replication.interval,n.replication,conf.level.percent,lower.conf.limit,upper.conf.limit)
     
     
     text.replication.interval <- sprintf("Thus, for r = %1.2f (when replications are conducted with N = %d) the replication interval is {%1.2f, %1.2f}.", r,n.replication,lower.replication.interval,upper.replication.interval)
     
     
     text.explain <- paste(text.ci,text.lower.population,text.upper.population, text.replication.interval)
     
     text.method.section <- sprintf("If our replication study obtains a correlation between r = %1.2f and r = %1.2f this will be evidence that our findings are not inconsistent with the original study correlation of r = %1.2f, %d%% CI [%1.2f,%1.2f]{%1.2f, %1.2f; N = %d}.",lower.replication.interval,upper.replication.interval,r,conf.level.percent,lower.conf.limit,upper.conf.limit,lower.replication.interval, upper.replication.interval,n.replication)
     
     
     
     text.out <- list()
     text.out$replication.interval <-text.replication.interval.brief
     text.out$replication.interval.explanation <- text.explain
     text.out$method.section.text <- text.method.section
     return(text.out)
}


#' Creates a replication interval based on a published/sample correlation. 
#' @param r Original study: Correlation
#' @param n Original study: Sample size 
#' @param rep.n (optional) Replication study: Sample size. If not specified, n is used.
#' @param conf.level (optional 0 to 1 value) Confidence level desired (0 to 1). If not specified .95 (i.e., 95 percent) will be used.
#' @param explain (optional boolean) Default is FALSE. If TRUE, text output explaining the interval is provided. 
#' @param extended.output (optional boolean) Default is FALSE. If TRUE, additional details (e.g., confidence interval) provided in numeric return output.
#' @param manuscript.text (optional boolean) Default is TRUE. If TRUE, present text for Method section.
#' @return A list of values (\code{lower.replication.interval.r, upper.replication.interval.r}) containing the replication interval (and related statistics if requested with the \code{extended.output} argument).
#' @examples
#' ri.r(r=.35,n=40)
#' ri.r(r=.35,n=40,explain=TRUE)
#' @export
ri.r <- function (r,n,rep.n=NA,conf.level = .95,explain=FALSE,extended.output=FALSE, manuscript.text = TRUE) {

     if (is.na(rep.n)) {
       rep.n <- n
     }
     
     r.value.ci=ci.r(r=r,n=n,conf.level=conf.level)
     r <- r.value.ci$r
     
     
     conf.lower.bound <- r.value.ci$lower.conf.limit.r
     conf.upper.bound <- r.value.ci$upper.conf.limit.r
     
     probability.in.tail=(1 - conf.level)/2
     
     lower.population.lower.bound <- quantile.r(probability.in.tail=probability.in.tail, ncr=conf.lower.bound,n=rep.n,is.lower.tail=TRUE)
     lower.population.upper.bound <- quantile.r(probability.in.tail=probability.in.tail, ncr=conf.lower.bound,n=rep.n,is.lower.tail=FALSE)
     
     
     upper.population.lower.bound <- quantile.r(probability.in.tail=probability.in.tail, ncr=conf.upper.bound,n=rep.n,is.lower.tail=TRUE)
     upper.population.upper.bound <- quantile.r(probability.in.tail=probability.in.tail, ncr=conf.upper.bound,n=rep.n,is.lower.tail=FALSE)
     
     
     
     
     
     interval.output <- list()
     interval.output$conf.level <- conf.level
     
     interval.output$r <- r
     interval.output$lower.confidence.interval.r <- conf.lower.bound
     interval.output$upper.confidence.interval.r <- conf.upper.bound
     interval.output$upper.confidence.interval.n <- n
     
     interval.output$lower.population.lower.bound.r <- lower.population.lower.bound
     interval.output$lower.population.upper.bound.r <- lower.population.upper.bound
     
     interval.output$upper.population.lower.bound.r <- upper.population.lower.bound
     interval.output$upper.population.upper.bound.r <- upper.population.upper.bound
     
     interval.output$replication.interval.n <- rep.n
     
     interval.output$lower.replication.interval.r <- lower.population.lower.bound
     interval.output$upper.replication.interval.r <- upper.population.upper.bound
         
     percent.ri.due.to.ci <- ((conf.upper.bound - conf.lower.bound) /(interval.output$upper.replication.interval.r - interval.output$lower.replication.interval.r))*100
     interval.output$percent.replication.interval.due.to.confidence.interval <- round(percent.ri.due.to.ci)
     
     interval.output.brief <- list()
     interval.output.brief$lower.confidence.interval.r <- conf.lower.bound
     interval.output.brief$upper.confidence.interval.r <- conf.upper.bound
     interval.output.brief$lower.replication.interval.r <- lower.population.lower.bound
     interval.output.brief$upper.replication.interval.r <- upper.population.upper.bound
    
      if (extended.output==TRUE) {
          interval.results <- interval.output
     } else {
          interval.results <- interval.output.brief
     }
     
     
     if (manuscript.text==TRUE | explain == TRUE) {
          rText<-ri.r.output(interval.output)
          cat("\n")
          cat(rText$replication.interval)
          cat("\n")
     }      
     
     if (manuscript.text == TRUE) {
          cat("\n")
          cat("Method Section Snippet:")
          cat("\n")
          cat(rText$method.section.text)
          cat("\n")
     }

     if (explain==TRUE) {
          cat("\nExplanation:\n")
          cat(rText$replication.interval.explanation)
          cat("\n\n")
     }
     
     if (manuscript.text==TRUE | explain == TRUE) {
          cat("\n\n")     
     }      
     
     
     return(interval.results)
}









