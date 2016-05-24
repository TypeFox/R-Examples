
diseq <- function(x, ...)
{
	UseMethod("diseq")
}

diseq.genotype <- function(x, ...)
  {
    if (nallele(x) < 2)
      {
        warning("Only 1 Marker allele. Returning NA")
        return(NA)
      }

    observed.no <- table( factor(allele(x,1), levels=allele.names(x)),
                          factor(allele(x,2), levels=allele.names(x)) )
    observed <- prop.table(observed.no)
    observed <- 1/2 * (observed + t(observed) )

    retval <- diseq.table(observed)
    retval$observed.no <- observed.no
    retval$call <- match.call()
    retval
  }

diseq.table <- function(x, ...)
{
  observed <- x
  allele.freq <- apply(observed,1,sum)
  # equal to: allele.freq <- apply(observed,2,sum)

  expected <- outer(allele.freq, allele.freq, "*")

  oeTab <- cbind(Obs=c(observed),
                 Exp=c(expected),
                 "Obs-Exp"=c(observed - expected))
  rownames(oeTab) <- outer(rownames(observed),
                           rownames(observed), paste, sep="/")

  diseq <- observed - expected
  diag(diseq) <- NA

  dmax.positive <- expected
  # equals: max( p(i)p(j), p(j)p(i) )

  dmax.negative <- outer(allele.freq, allele.freq, pmin ) - expected
  # equals: min( p(i) * (1 - p(j)), p(j)( 1 - (1-p(i) ) ) )

  dprime <- diseq / ifelse( diseq > 0, dmax.positive, dmax.negative )

  # r gives the pairwise correlation coefficient for pairs containing at lease
  # one allele from the specified pair.
  # For two alleles:
  #    corr coefficient = diseq / sqrt( p(a) * (1-p(a) ) * p(b) * (1-p(b)) )
  #p.1.minus.p <- allele.freq * (1-allele.freq)
  #r <-  -diseq  / sqrt( outer( p.1.minus.p, p.1.minus.p, "*") )
  r.denom <- sqrt(    allele.freq  %*% t(  allele.freq)) *
             sqrt( (1-allele.freq) %*% t(1-allele.freq))
  r <- -diseq  / r.denom

  # above formula works unchanged for 2 alleles, but requires adjustment
  # for multiple alleles.
  # r <- r * (length(allele.freq) - 1)

  offdiag.expected <- expected
  diag(offdiag.expected) <- NA
  sum.expected <- sum(offdiag.expected, na.rm=TRUE)

  if(all(dim(x)==2)) # 2 allele case
    {
      diseq.overall <- diseq[1,2]
      dprime.overall <- dprime[1,2]
      r.overall <- r[1,2]
      R2.overall <- r.overall^2
    }
  else
    {
      diseq.overall <- sum( abs(diseq) * expected , na.rm=TRUE ) / sum.expected
      dprime.overall <- sum( abs(dprime) * expected , na.rm=TRUE ) / sum.expected
      r.overall <- sum( abs(r) * expected , na.rm=TRUE ) / sum.expected
      R2.overall <- r.overall^2
    }

  #diag(r) <- 1.0

  retval <- list(
                 call = match.call(),
                 observed=observed,
                 expected=expected,
                 table=oeTab,
                 allele.freq=allele.freq,
                 D=diseq,
                 Dprime=dprime,
                 r=r,
                 R2=r^2,
                 D.overall=diseq.overall,
                 Dprime.overall=dprime.overall,
                 r.overall = r.overall,
                 R2.overall = R2.overall
                 )

  class(retval) <- "diseq"
  retval
}

print.diseq  <-  function(x, show=c("D","D'","r","R^2","table"), ...)
  {

    cat("\n")
    if(!is.null(x$locus))
      {
        cat("\n")
        print( x$locus )
      }
    cat("\n")
    cat("Call: \n")
    print(x$call)
    cat("\n")
    if("D" %in% show)
      {
        cat("Disequlibrium for each allele pair (D)\n")
        cat("\n")
        print(x$D)
        cat("\n")
      }
    if("D'" %in% show)
      {
        cat("Disequlibrium for each allele pair (D')\n")
        cat("\n")
        print(x$Dprime)
        cat("\n")
      }
    if("r" %in% show)
      {
        cat("Correlation coefficient for each allele pair (r)\n")
        cat("\n")
        print(x$r)
        cat("\n")
      }
    if("R^2" %in% show)
      {
        cat("R^2 for each allele pair\n")
        cat("\n")
        print(x$R2)
        cat("\n")
      }
    if("table" %in% show)
      {
        cat("Observed vs Expected frequency table\n")
        cat("\n")
        print(x$table)
        cat("\n")
      }

    if( any(c("D","D'","r") %in% show))
      {
        if( ncol(x$r) <= 2 )
          cat("Overall Values\n")
        else
          cat("Overall Values (mean absolute-value weighted by expected allele frequency)\n")
        cat("\n")

        if("D" %in% show)
          cat("  D  :  ", x$D.overall, "\n", sep="")
        if("D'" %in% show)
          cat("  D' :  ", x$Dprime.overall, "\n", sep="")
        if("r" %in% show)
          cat("  r  :  ", x$r.overall, "\n", sep="")
        if("R^2" %in% show)
          cat("  R^2:  ", x$R2.overall, "\n", sep="")
        cat("\n")
      }

    cat("\n")
  }

diseq.ci <- function(x, R=1000, conf=0.95, correct=TRUE, na.rm=TRUE, ...)
{
  if (!("genotype") %in% class(x) )
    stop("x must inherit from class 'genotype'.")

  if( any(is.na(x) ) )
    {
      if( na.rm)
        x <- na.omit(x)
      else
        stop("Missing values and NaN's not allowed if `na.rm' is FALSE.")
    }

  # step 1 - generate summary table
  observed.no <- table( factor(allele(x,1), levels=allele.names(x)),
                        factor(allele(x,2), levels=allele.names(x)) )
  observed <- prop.table(observed.no)
  observed <- 1/2 * (observed + t(observed) )

  # step 2 - make table into a probability vector for calling rmultinom
  n <- sum(observed.no)
  prob.vector <- c(observed)

  # step 3 - sample R multinomials with the specified frequenceis
  # (include observed data to avoid bias)
  resample.data <- cbind(c(observed.no),
                         rmultz2( n, prob.vector, R ) )

  bootfun <- function(x) {
    observed[,] <- x/n
    observed <- 1/2 * (observed + t(observed) )
    d <-  diseq(observed)
    c( "Overall D  "=d$D.overall,
       "Overall D' "=d$Dprime.overall,
       "Overall r  "=d$r.overall,
       "Overall R^2"=d$R2.overall)
  }

  results <- apply( resample.data, 2, bootfun )

  alpha.2 <- (1-conf)/2

#  ci <- t(apply(results, 1,
#              quantile, c( alpha.2 , 1-alpha.2), na.rm=TRUE ))

  if(length(allele.names(x))<=2)
    {
      ci <- t(apply(results, 1, function(x) quantile(x, c(0.025, 0.975),
                                                     na.rm=na.rm ) ) )

      warning.text <- paste("The R^2 disequlibrium statistics is bounded",
                            "between [0,1].  The confidence ",
                            "intervals for R^2 values near 0 and 1 are",
                            "ill-behaved.", sep=" ")

      if(correct)
        {
          warning.text <- paste(warning.text, "A rough correction has",
                                "been applied, but the intervals still",
                                "may not be correct for R^2 values near",
                                "0 or 1.",
                                sep=" ")

          X <- results["Overall R^2",]
          ci["Overall R^2",] <- ci.balance(X,X[1],confidence=conf,
                                           minval=0,maxval=1)$ci
        }
    }
  else
    {
      warning.text <- paste("For more than two alleles, overall",
                            "disequlibrium statistics are bounded",
                            "between [0,1].  Because of this, confidence",
                            "intervals for values near 0 and 1 are",
                            "ill-behaved.", sep=" ")

      if(correct)
        {
          warning.text <- paste(warning.text, "A rough correction has been applied, but",
                                "the intervals still may not be correct for values near 0 or 1.",
                                sep=" ")

          ci <- t(apply(results, 1,
                        function(x)
                        ci.balance(x,x[1],confidence=conf,
                                   minval=0,maxval=1)$ci ))
        }
      else
        ci <- t(apply(results, 1, function(x) quantile(x, c(0.025, 0.975) ) ) )

      warning(paste(strwrap(c(warning.text,"\n"),prefix="  "),collapse="\n") )
    }

  na.count <-  function(x) sum(is.na(x))
  nas <- apply( results, 1, na.count)

  zero.in.range <- (ci[,1] <= 0) & (ci[,2] >= 0)

  ci <- cbind( "Observed"=results[,1], ci, "NAs"=nas,
               "Zero in Range"=zero.in.range )

  outside.ci <- (ci[,1] < ci[,2]) | (ci[,1] > ci[,3])

  if( any(outside.ci) )
    warning("One or more observed value outide of confidence interval. Check results.")

  if(any(nas>0))
    warning("NAs returned from diseq call")


  retval <- list(
         call=match.call(),
         R=R,
         conf=conf,
         ci=ci,
         warning.text=warning.text
         )

  retval
}
