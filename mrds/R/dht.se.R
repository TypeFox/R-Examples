#' Variance and confidence intervals for density and abundance estimates
#'
#' Computes standard error, cv, and log-normal confidence intervals for
#' abundance and density within each region (if any) and for the total of all
#' the regions.  It also produces the correlation matrix for regional and total
#' estimates.
#'
#' The variance has two components: 1) variation due to uncertanity from
#' estimation of the detection function and 2) variation in abundance due to
#' random sample selection.  The first component is computed using a delta
#' method estimate of variance (\code{\link{DeltaMethod}} (Huggins 1989, 1991,
#' Borchers et al. 1998) in which the first derivatives of the abundance
#' estimator with respect to the parameters in the detection function are
#' computed numerically.  The second component can be computed in one of three
#' ways as set by the option \code{varflag} with values 0,1,2.
#'
#' A value of 0 is to use a binomial variance for the number of observations
#' and it is only useful if the sampled region is the survey region and the
#' objects are not clustered which will not occur very often. If covered region
#' is less than the survey region the variance estimator is scaled up but it
#' will be a poor estimator and the confidence interval will likely not achieve
#' the nominal level.
#'
#' A value of 1 uses the variance for the encounter rate of (Fewster et al.
#' (2009), estimator R2 (which has been shown to have better properties than
#' the previous default of Buckland et al. 2001 pg 78-79)).  If
#' \code{group=FALSE} the variance of the mean group size is also included.
#' This variance estimator is not appropriate if \code{size} or a derivative of
#' \code{size} is used in the any of the detection function models.
#'
#' In general if any covariates are used in the models, the default option 2 is
#' preferable.  It uses a variance estimator based on that suggested by Innes
#' et al. (2002) which used the formula for the variance ecounter rate but
#' replaces the number of observations per sample with the estimated abundance
#' per sample.  The difference between the version used here and that in Innes
#' et al. (2002) is that Innes et al. use an estimator with form similar to
#' that of Buckland et al. (2001), while the estimator here uses a form based
#' on Fewster et al. (2009, estimator R2).
#'
#' For more on encounter rate variance estimation, see \link{varn}.
#'
#' Exceptions to the above occur if there is only one sample in a stratum. In
#' that case it uses Poisson assumption (var(x)=x) and it assumes a known
#' variance so z=1.96 is used for critical value. In all other cases the
#' degrees of freedom for the t-distribution assumed for the log(abundance) or
#' log(density) is based on the Satterthwaite approximation (Buckland et al.
#' 2001 pg 90) for the degrees of freedom (df).  The df are weighted by the
#' squared cv in combining the two sources of variation because of the assumed
#' log-normal distribution because the components are multiplicative.  For
#' combining df for the sampling variance across regions they are weighted by
#' the variance because it is a sum across regions.
#'
#' A non-zero correlation between regional estimates can occur from using a
#' common detection function across regions.  This is reflected in the
#' correlation matrix of the regional and total estimates which is given in the
#' value list.  It is only needed if subtotals of regional estimates are
#' needed.
#'
#' @param model ddf model object
#' @param region.table table of region values
#' @param samples table of samples(replicates)
#' @param obs table of observations
#' @param options list of options that can be set (see \code{\link{dht}})
#' @param numRegions number of regions
#' @param estimate.table table of estimate values
#' @param Nhat.by.sample estimated abundances by sample
#' @export
#' @return List with 2 elements: \item{estimate.table}{completed table with se,
#'   cv and confidence limits} \item{vc }{correlation matrix of estimates}
#' @note This function is called by \code{dht} and it is not expected that the
#'   user will call this function directly but it is documented here for
#'   completeness and for anyone expanding the code or using this function in
#'   their own code
#' @author Jeff Laake
#' @seealso \code{\link{dht}}, \code{\link{print.dht}}
#' @references see \code{\link{dht}}
#' @keywords utility
#' @importFrom stats qnorm qt var
dht.se <- function(model, region.table, samples, obs, options, numRegions,
                   estimate.table, Nhat.by.sample){
  #  Functions Used:  DeltaMethod, dht.deriv (in DeltaMethod), varn

  # Define function: compute.df
  compute.df<- function(k,type){
    if(type=="O1" | type=="O2"| type=="O3"){
      H.O <- k - 1
      k.h.O <- rep(2, H.O)
      df <- sum(k.h.O - 1)
    }else{
      if(type=="S1" | type=="S2"){
        H.S <- floor(k/2)
        k.h.S <- rep(2, H.S)
        if(k %% 2 > 0) k.h.S[H.S] <- 3
        df <- sum(k.h.S - 1)
      }else{
        df <- k-1
      }
    }
    return(df)
  }

  # First compute variance component due to estimation of detection function
  # parameters. This uses the delta method and produces a v-c matrix if more
  # than one strata
  if(!is.null(model$par)){
    vcov <- solvecov(model$hessian)$inv
    vc1.list <- DeltaMethod(model$par, dht.deriv, vcov, options$pdelta,
        model = model, samples = samples, obs = obs, options = options)
    vc1 <- vc1.list$variance
  }else{
    vc1.list <- list(variance=0)
    vc1 <- 0
  }

  # Next compute the component due to sampling of both lines and of the
  # detection process itself
  # There are 3 different options here:
  #  1) varflag=0; Binomial variance of detection process - only applicable if
  #   survey region=covered region although it will scale up but it would be
  #   a poor estimator
  #  2) varflag=1; delta method, with varn based on Fewster et al (2009)
  #   estimator R2 (var(n/L))
  #  3) varflag=2; Innes et al variance estimator (var(N/L), except changed to
  #   resemble the form of estimator R2 of Fewster et al (2009))
  # Exceptions to the above occur if there is only one sample in a stratum.
  #   In that case it uses Poisson approximation.

  scale <- region.table$Area/region.table$CoveredArea

  # If no areas were given or varflag=0 use approach #1 (varflag=0)
  # Note: vc2 is of proper dimension because Region.Label for obs is setup
  # with all levels of the Region.Label from the region.table.
  if(sum(region.table$Area) == 0 | options$varflag == 0){
    if(options$group){
      vc2 <- by((1 - obs$pdot)/obs$pdot^2, obs$Region.Label,sum)
    }else{
      vc2 <- by(obs$size^2 * (1 - obs$pdot)/obs$pdot^2, obs$Region.Label, sum)
    }
    # set missing value to 0
    vc2[is.na(vc2)] <- 0

    if(sum(region.table$Area) != 0){
      vc2 <- vc2 * scale[1:numRegions]^2
    }
  }else{
    # Otherwise compute variance for varflag=1 or 2
    vc2 <- rep(0, numRegions)
    # 26 jan 06 jll; changed to use object rather than distance; also
    # overwrites existing n because that can be sum(size) rather than count
    nobs <- tapply(obs$object, obs$Label, length)
    nobs <- data.frame(Label = names(nobs),
                       n = as.vector(nobs)[!is.na(nobs)])
    Nhat.by.sample$n <- NULL
    # when there are no sighings
    if(nrow(nobs) > 0){
      Nhat.by.sample <- merge(Nhat.by.sample, nobs, by.x = "Label",
          by.y = "Label", all.x = TRUE)
      Nhat.by.sample$n[is.na(Nhat.by.sample$n)] <- 0
    }else{
      Nhat.by.sample <- cbind(Nhat.by.sample, n = rep(0,nrow(Nhat.by.sample)))
    }

    # Compute number of lines per region for df calculation
    if(numRegions > 1){
      estimate.table$k <- c(as.vector(tapply(samples$Effort,
                                             samples$Region.Label, length)), 0)
      estimate.table$k[numRegions + 1] <- sum(estimate.table$k)
    }else{
      estimate.table$k <- as.vector(tapply(samples$Effort,
                                           samples$Region.Label, length))
    }

    # If individual abundance being computed, calculate mean and variance
    # of mean for group size.
    if(!options$group){
      if(length(obs$size) > 0){
        vars <- by(obs$size, obs$Region.Label, var)/
                  by(obs$size, obs$Region.Label, length)
        sbar <- by(obs$size, obs$Region.Label, mean)
        sobs <- data.frame(Region.Label = names(sbar),
                           vars         = as.vector(vars),
                           sbar         = as.vector(sbar))
      }else{
        sobs = data.frame(Region.Label=levels(obs$Region.Label),
                          vars=rep(NA,length(levels(obs$Region.Label))),
                          sbar=rep(NA,length(levels(obs$Region.Label))))
      }
      Nhat.by.sample <- merge(Nhat.by.sample, sobs, by.x = "Region.Label",
                              by.y = "Region.Label", all.x = TRUE)
      Nhat.by.sample$sbar[is.na(Nhat.by.sample$sbar)] <- 0
      Nhat.by.sample$vars[is.na(Nhat.by.sample$vars)] <- 0
    }else{
      # If group abundance is being estimated, set mean=1, var=0
      Nhat.by.sample$sbar <- rep(1, dim(Nhat.by.sample)[1])
      Nhat.by.sample$vars <- rep(0, dim(Nhat.by.sample)[1])
    }

    # sort Nhat.by.sample by Region.Label and Sample.Label
    Nhat.by.sample <- Nhat.by.sample[order(Nhat.by.sample$Region.Label,Nhat.by.sample$Sample.Label),]

    # Loop over each region and compute each variance;
    # jll 11/11/04 - changes made in following code using
    # Effort.x (effort per line) rather than previous errant code
    # that used Effort.y (effort per region)
    for(i in 1:numRegions){
      stratum.data <- Nhat.by.sample[as.character(Nhat.by.sample$Region.Label)==
                                     as.character(region.table$Region[i]), ]
      Ni <- sum(stratum.data$Nhat)
      Li <- sum(stratum.data$Effort.x)
      sbar <- stratum.data$sbar[1]
      vars <- stratum.data$vars[1]

      if (options$group) vars <- 0

      if(length(stratum.data$Effort.y) == 1){
        if (options$varflag == 1){
          vc2[i] <- Ni^2 * (1/stratum.data$n + vars/sbar^2)
        }else{
          vc2[i] <- Ni^2 * (1/Ni + vars/sbar^2)
        }
      }else if (options$varflag == 1){
        vc2[i] <- (Ni * Li)^2 * varn(stratum.data$Effort.x,
                                     stratum.data$n,type=options$ervar)/
                                sum(stratum.data$n)^2 + Ni^2 * vars/sbar^2
      }else{
        if(options$varflag==2){
          vc2[i] <- varn(stratum.data$Effort.x/(scale[i] * Li),
                         stratum.data$Nhat/scale[i],type=options$ervar)
        }else{
          vc2[i] <- varn(stratum.data$Effort.x/(scale[i] * Li),
                        stratum.data$Nhat/scale[i], type=options$ervar)
        }
      }
    }
  }

  vc2[is.nan(vc2)] <- 0

  # Construct v-c matrix for encounter rate variance given computed ps
  # The cov between regional estimate and total estimate is simply var for
  #  regional estimate.
  # Assumes no cov between regions due to independent sample selection.
  if(numRegions > 1){
    v2 <- vc2
    vc2 <- diag(c(vc2, sum(vc2)))
    vc2[1:numRegions, (numRegions + 1)] <- v2
    vc2[(numRegions + 1), 1:numRegions] <- v2
  }else if (length(vc2) > 1){
    vc2 <- diag(vc2)
  }else{
    vc2 <- as.matrix(vc2)
  }

  vc <- vc1 + vc2

  # deal with missing values and 0 estimates.
  estimate.table$se <- sqrt(diag(vc))
  estimate.table$se[is.nan(estimate.table$se)] <- 0
  estimate.table$cv <- estimate.table$se/estimate.table$Estimate
  estimate.table$cv[is.nan(estimate.table$cv)] <- 0

  # work out the confidence intervals
  # if the options$ci.width is set, then use that, else default to
  # 95% CI
  if(is.null(options$ci.width)){
    ci.width <- 0.025
  }else{
    ci.width <- (1-options$ci.width)/2
  }

  # Use satterthwaite approx for df and log-normal distribution for
  # 95% intervals
  if(options$varflag != 0){
    # set df from replicate lines to a min of 1 which avoids divide by zero
    # Following 2 lines added and references to estimate.table$k changed to df
    df <- estimate.table$k
    df <- sapply(df,compute.df,type=options$ervar)
    df[df<1] <- 1

    if(is.na(vc1) || all(vc1==0)){
      estimate.table$df <- df
    }else{
      estimate.table$df <- estimate.table$cv^4/((diag(vc1)/
                            estimate.table$Estimate^2)^2/(length(model$fitted) -
                            length(model$par)) +
                            (diag(vc2)/estimate.table$Estimate^2)^2/df)
    }

    # compute proper satterthwaite
    # df for total estimate assuming sum of indep region estimates; uses
    # variances instead of cv's because it is a sum of means for encounter
    # rate portion of variance (df.total)
    if(numRegions>1){
      df.total <- (diag(vc2)[numRegions+1])^2/
                   sum((diag(vc2)^2/df)[1:numRegions])
      if(all(vc1==0)){
          estimate.table$df[numRegions+1] <- df.total
      }else{
        estimate.table$df[numRegions+1] <- estimate.table$cv[numRegions+1]^4 /
                    ((diag(vc1)[numRegions+1]/
                    estimate.table$Estimate[numRegions+1]^2)^2/
                    (length(model$fitted)-length(model$par))
                    + (diag(vc2)[numRegions+1]/
                        estimate.table$Estimate[numRegions+1]^2)^2/df.total)
      }
    }

    estimate.table$df[estimate.table$df < 1 &estimate.table$df >0] <- 1
    cvalue <- exp((abs(qt(ci.width, estimate.table$df)) *
                  sqrt(log(1 + estimate.table$cv^2))))
  }else{
    # intervals for varflag=0; sets df=0
    # and uses normal approximation
    cvalue <- exp((abs(qnorm(ci.width)) * sqrt(log(1 + estimate.table$cv^2))))
    estimate.table$df <- rep(0,dim(estimate.table)[1])
  }

  # deal with missing values and divide by 0 issues
  estimate.table$df[is.nan(estimate.table$df)] <- 0
  estimate.table$lcl  <-  estimate.table$Estimate/cvalue
  estimate.table$lcl[is.nan(estimate.table$lcl)] <- 0
  estimate.table$ucl  <-  estimate.table$Estimate * cvalue
  estimate.table$ucl[is.nan(estimate.table$ucl)] <- 0
  estimate.table$k  <-  NULL

  return(list(estimate.table = estimate.table,
              vc             = vc,
              vc1            = vc1.list,
              vc2            = vc2 ))
}
