cdf.est <- function(z, wgt, x=NULL, y=NULL, stratum=NULL, cluster=NULL,
   wgt1=NULL, x1=NULL, y1=NULL, popsize=NULL, popcorrect=FALSE, pcfsize=NULL,
   N.cluster=NULL, stage1size=NULL, support=NULL, sizeweight=FALSE, swgt=NULL,
   swgt1=NULL, vartype="Local", conf=95, cdfval=NULL, pctval=c(5,10,25,50,75,90,95),
   check.ind=TRUE, warn.ind=NULL, warn.df=NULL, warn.vec=NULL) {

################################################################################
# Function: cdf.est
# Programmer: Tom Kincaid
# Date: July 12, 2000
# Last Revised: October 10, 2012
# Description:
#   This function calculates an estimate of the cumulative distribution function
#   (CDF) for the proportion (expressed as percent) and the total of a response 
#   variable, where the response variable may be defined for either a finite
#   or an extensive resource.  Optionally, for a finite resource, the size-
#   weighted CDF can be calculated.  In addition the standard error of the
#   estimated CDF and confidence bounds are calculated.  The user can supply the
#   set of values at which the CDF is estimated.  For the CDF of a proportion,
#   the Horvitz-Thompson ratio estimator, i.e., the ratio of two Horvitz-Thompson
#   estimators, is used to calculate the CDF estimate.  For the CDF of a total,
#   the user can supply the known size of the resource or the known sum of the
#   size-weights of the resource, as appropriate.  For the CDF of a total when
#   either the size of the resource or the sum of the size-weights of the
#   resource is provided, the classic ratio estimator is used to calculate the
#   CDF estimate, where that estimator is the product of the known value and the
#   Horvitz-Thompson ratio estimator.   For the CDF of a total when neither the
#   size of the resource nor the sum of the size-weights of the resource is
#   provided, the Horvitz-Thompson estimator is used to calculate the CDF
#   estimate.  Variance estimates for the estimated CDF are calculated using
#   either the local mean variance estimator or the simple random sampling (SRS)
#   variance estimator.  The choice of variance estimator is subject to user
#   control. The local mean variance estimator requires the x-coordinate and the
#   y-coordinate of each site.  The SRS variance estimator uses the independent
#   random sample approximation to calculate joint inclusion probabilities.
#   Confidence bounds are calculated using a Normal distribution multiplier.  
#   In addition the function uses the estimated CDF to calculate percentile
#   estimates.  Estimated confidence bounds for the percentile estimates are
#   calculated.  The user can supply the set of values for which percentiles
#   estimates are desired.  Optionally, the user can use the default set of
#   percentiles.  The function can accommodate a stratified sample.  For a
#   stratified sample, separate estimates and standard errors are calculated
#   for each stratum, which are used to produce estimates and standard errors 
#   for all strata combined.  Strata that contain a single value are removed.
#   For a stratified sample, when either the size of the resource or the sum of
#   the size-weights of the resource is provided for each stratum, those values
#   are used as stratum weights for calculating the estimates and standard 
#   errors for all strata combined.  For a stratified sample when neither the
#   size of the resource nor the sum of the size-weights of the resource is
#   provided for each stratum, estimated values are used as stratum weights for
#   calculating the estimates and standard errors for all strata combined.  The
#   function can accommodate single-stage and two-stage samples for both
#   stratified and unstratified sampling designs.  Finite population and
#   continuous population correction factors can be utilized in variance
#   estimation.  The function checks for compatibility of input values and
#   removes missing values.
# Arguments:
#   z = the response value for each site.
#     for each site.
#   wgt = the final adjusted weight (inverse of the sample inclusion
#     probability) for each site, which is either the weight for a single-stage
#     sample or the stage two weight for a two-stage sample.
#   x = x-coordinate for location for each site, which is either the
#     x-coordinate for a single-stage sample or the stage two x-coordinate for a
#     two-stage sample.  The default is NULL.
#   y = y-coordinate for location for each site, which is either the
#     y-coordinate for a single-stage sample or the stage two y-coordinate for a
#     two-stage sample.  The default is NULL.
#   stratum = the stratum for each site.  The default is NULL.
#   cluster = the stage one sampling unit (primary sampling unit or cluster)
#     code for each site.  The default is NULL.
#   wgt1 = the final adjusted stage one weight for each site.  The default is
#     NULL.
#   x1 = the stage one x-coordinate for location for each site.  The default is
#     NULL.
#   y1 = the stage one y-coordinate for location for each site.  The default is
#     NULL.
#   popsize = known size of the resource, which is used to perform ratio
#     adjustment to estimators expressed using measurement units for the
#     resource and to calculate strata proportions for calculating estimates for
#     a stratified sample.  For a finite resource, this argument is either the
#     total number of sampling units or the known sum of size-weights.  For an
#     extensive resource, this argument is the measure of the resource, i.e.,
#     either known total length for a linear resource or known total area for an
#     areal resource.  For a stratified sample this variable must be a vector
#     containing a value for each stratum and must have the names attribute set
#     to identify the stratum codes.  The default is NULL.
#   popcorrect = a logical value that indicates whether finite or continuous
#     population correction factors should be employed during variance
#     estimation, where TRUE = use the correction factors and FALSE = do not use
#     the correction factors.  The default is FALSE.
#   pcfsize = size of the resource, which is required for calculation of finite
#     and continuous population correction factors for a single-stage sample.
#     For a stratified sample this argument must be a vector containing a value
#     for each stratum and must have the names attribute set to identify the
#     stratum codes.  The default is NULL.
#   N.cluster = the number of stage one sampling units in the resource, which is
#     required for calculation of finite and continuous population correction
#     factors for a two-stage sample.  For a stratified sample this argument
#     must be a vector containing a value for each stratum and must have the
#     names attribute set to identify the stratum codes.  The default is NULL.
#   stage1size = size of the stage one sampling units of a two-stage sample,
#     which is required for calculation of finite and continuous population
#     correction factors for a two-stage sample and must have the names
#     attribute set to identify the stage one sampling unit codes.  For a
#     stratified sample, the names attribute must be set to identify both
#     stratum codes and stage one sampling unit codes using a convention where
#     the two codes are separated by the & symbol, e.g., "Stratum 1&Cluster 1".
#     The default is NULL.
#   support = the support value for each site - the value one (1) for a site
#     from a finite resource or the measure of the sampling unit associated with
#     a site from an extensive resource, which is required for calculation of
#     finite and continuous population correction factors.  The default is NULL.
#   sizeweight = a logical value that indicates whether size-weights should be
#     used in the analysis, where TRUE = use the size-weights and FALSE = do not
#     use the size-weights.  The default is FALSE.
#   swgt = the size-weight for each site, which is the stage two size-weight for
#     a two-stage sample.  The default is NULL.
#   swgt1 = the stage one size-weight for each site.  The default is NULL.
#   vartype = the choice of variance estimator, where "Local" = local mean
#     estimator and "SRS" = SRS estimator.  The default is "Local".
#   conf = the confidence level.  The default is 95%.
#   cdfval = the set of values at which the CDF is estimated.  If a set of
#    values is not provided, then the sorted set of unique values of the
#    response variable is used.  The default is NULL.
#   pctval = the set of values at which percentiles are estimated.  The default
#    set is: {5, 10, 25, 50, 75, 90, 95}.
#   check.ind = a logical value that indicates whether compatability checking of
#     the input values is conducted, where TRUE = conduct compatibility checking
#     and FALSE = do not conduct compatibility checking.  The default is TRUE.
#   warn.ind = a logical value that indicates whether warning messages were
#     generated, where TRUE = warning messages were generated and FALSE =
#     warning messages were not generated.  The default is NULL.
#   warn.df = a data frame for storing warning messages.  The default is NULL.
#   warn.vec = a vector that contains names of the population type, the
#     subpopulation, and an indicator.  The default is NULL.
# Output:
#   If the function was called by the cont.analysis function, then output is an
#   object in list format composed of a list named Results, which contains
#   estimates and confidence bounds, and a data frame named warn.df, which
#   contains warning messages.  The Results list is composed of two data frames:
#   one data frame named CDF, which contains the CDF estimates, and a second
#   data frame named Pct, which contains the percentile estimates.  If the
#   function was called directly, then output is the Results list.
# Other Functions Required:
#   input.check - check input values for errors, consistency, and compatibility
#     with analytical functions
#   wnas - remove missing values
#   vecprint - takes an input vector and outputs a character string with line
#     breaks inserted
#   cdf.nresp - calculate the number of response values less than or equal to
#     each of the set of values at which the CDF is estimated
#   cdf.prop - calculate the CDF for the proportion
#   cdf.total - calculate the CDF for the total
#   cdf.size.prop - calculate the size-weighted CDF for the proportion
#   cdf.size.total - calculate the size-weighted CDF for the total
#   cdfvar.prop - calculate variance of the CDF for the proportion
#   cdfvar.total - calculate variance of the CDF for the total
#   cdfvar.size.prop - calculate variance of the size-weighted CDF for the
#     proportion
#   cdfvar.size.total - calculate variance of the size-weighted CDF for the
#     total
# Examples:
#   z <- rnorm(100, 10, 1)
#   wgt <- runif(100, 10, 100)
#   cdfval <- seq(min(z), max(z), length=20)
#   cdf.est(z, wgt, vartype="SRS", cdfval=cdfval)
#
#   x <- runif(100)
#   y <- runif(100)
#   cdf.est(z, wgt, x, y, cdfval=cdfval)
################################################################################

# As necessary, create a data frame for warning messages

   if(is.null(warn.ind)) {
      warn.ind <- FALSE
      warn.df <- NULL
      warn.vec <- rep(NA, 3)
   }
   fname <- "cdf.est"

# Check for existence of the response variable and determine the number of
# values

   if(is.null(z))
      stop("\nValues must be provided for the response variable.")
   if(!is.numeric(z))
      stop("\nValues for the response variable must be numeric.")
   nresp <- length(z)

# Assign a logical value to the indicator variable for a stratified sample

   stratum.ind <- length(unique(stratum)) > 1

# If the sample is stratified, convert stratum to a factor, determine stratum 
# levels, and calculate number of strata,

   if(stratum.ind) {
      stratum <- factor(stratum)
      stratum.levels <- levels(stratum)
      nstrata <- length(stratum.levels)
   } else {
      stratum.levels <- NULL
      nstrata <- NULL
   }

# Assign a logical value to the indicator variable for a two stage sample

   cluster.ind <- length(unique(cluster)) > 1

# Assign the value of popcorrect to the indicator variable for use of the
# population correction factor

   pcfactor.ind <- popcorrect

# Assign the value of sizeweight to the indicator variable for use of size
# weights

   swgt.ind <- sizeweight

# Begin the section that checks for compatibility of input values

   if(check.ind) {

# If the sample has two stages, convert cluster to a factor, determine cluster 
# levels, and calculate number of clusters

   if(cluster.ind) {
      if(stratum.ind) {
         cluster.in <- cluster
         cluster <- tapply(cluster, stratum, factor)
         cluster.levels <- sapply(cluster, levels, simplify=FALSE)
         ncluster <- sapply(cluster.levels, length)
      } else {
         cluster <- factor(cluster)
         cluster.levels <- levels(cluster)
         ncluster <- length(cluster.levels)
      }
   }

# Check for compatibility of input values

      temp <- input.check(nresp, wgt, NULL, NULL, x, y, stratum.ind, stratum,
         stratum.levels, nstrata, cluster.ind, cluster, cluster.levels,
         ncluster, wgt1, x1, y1, popsize, pcfactor.ind, pcfsize, N.cluster,
         stage1size, support, swgt.ind, swgt, swgt1, vartype, conf, cdfval,
         pctval)

      popsize <- temp$popsize
      pcfsize <- temp$pcfsize
      N.cluster <- temp$N.cluster
      stage1size <- temp$stage1size

# If the sample was stratified and had two stages, then reset cluster to its 
# input value

   if(stratum.ind && cluster.ind)
      cluster <- cluster.in

# End the section that checks for compatibility of input values

   }

# Remove missing values

   if(vartype == "Local") {
      if(swgt.ind) {
         if(stratum.ind) {
            if(cluster.ind)
               temp <- wnas(list(z=z, wgt=wgt, x=x, y=y, stratum=stratum,
                   cluster=cluster, wgt1=wgt1, x1=x1, y1=y1, swgt=swgt,
                   swgt1=swgt1))
            else
               temp <- wnas(list(z=z, wgt=wgt, x=x, y=y, stratum=stratum,
                  swgt=swgt))
         } else {
            if(cluster.ind)
               temp <- wnas(list(z=z, wgt=wgt, x=x, y=y, cluster=cluster,
                  wgt1=wgt1, x1=x1, y1=y1, swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(z=z, wgt=wgt, x=x, y=y, swgt=swgt))
         }
         z <- temp$z
         wgt <- temp$wgt
         x <- temp$x
         y <- temp$y
         if(stratum.ind)
            stratum <- temp$stratum
         if(cluster.ind) {
            cluster <- temp$cluster
            wgt1 <- temp$wgt1
            x1 <- temp$x1
            y1 <- temp$y1
            swgt1 <- temp$swgt1
         }
         swgt <- temp$swgt
      } else {
         if(stratum.ind) {
            if(cluster.ind)
               temp <- wnas(list(z=z, wgt=wgt, x=x, y=y, stratum=stratum,
                  cluster=cluster, wgt1=wgt1, x1=x1, y1=y1))
            else
               temp <- wnas(list(z=z, wgt=wgt, x=x, y=y, stratum=stratum))
         } else {
            if(cluster.ind)
               temp <- wnas(list(z=z, wgt=wgt, x=x, y=y, cluster=cluster,
                  wgt1=wgt1, x1=x1, y1=y1))
            else
               temp <- wnas(list(z=z, wgt=wgt, x=x, y=y))
         }
         z <- temp$z
         wgt <- temp$wgt
         x <- temp$x
         y <- temp$y
         if(stratum.ind)
            stratum <- temp$stratum
         if(cluster.ind) {
            cluster <- temp$cluster
            wgt1 <- temp$wgt1
            x1 <- temp$x1
            y1 <- temp$y1
         }
      }
   } else {
      if(swgt.ind) {
         if(stratum.ind) {
            if(cluster.ind)
               temp <- wnas(list(z=z, wgt=wgt, stratum=stratum, cluster=cluster,
                  wgt1=wgt1, swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(z=z, wgt=wgt, stratum=stratum, swgt=swgt))
         } else {
            if(cluster.ind)
               temp <- wnas(list(z=z, wgt=wgt, cluster=cluster, wgt1=wgt1,
                  swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(z=z, wgt=wgt, swgt=swgt))
         }
         z <- temp$z
         wgt <- temp$wgt
         if(stratum.ind)
            stratum <- temp$stratum
         if(cluster.ind) {
            cluster <- temp$cluster
            wgt1 <- temp$wgt1
            swgt1 <- temp$swgt1
         }
         swgt <- temp$swgt
      } else {
         if(stratum.ind) {
            if(cluster.ind)
               temp <- wnas(list(z=z, wgt=wgt, stratum=stratum, cluster=cluster,
                  wgt1=wgt1))
            else
               temp <- wnas(list(z=z, wgt=wgt, stratum=stratum))
         } else {
            if(cluster.ind)
               temp <- wnas(list(z=z, wgt=wgt, cluster=cluster, wgt1=wgt1))
            else
               temp <- wnas(list(z=z, wgt=wgt))
         }
         z <- temp$z
         wgt <- temp$wgt
         if(stratum.ind)
            stratum <- temp$stratum
         if(cluster.ind) {
            cluster <- temp$cluster
            wgt1 <- temp$wgt1
         }
      }
   }

# For a stratified sample, check for strata that no longer contain any values,
# as necesssary adjust popsize, remove strata that contain a single value, and
# output a warning message

   if(stratum.ind) {
      stratum <- factor(stratum)
      stratum.levels.old <- stratum.levels
      stratum.levels <- levels(stratum)
      nstrata.old <- nstrata
      nstrata <- length(stratum.levels)
      if(nstrata < nstrata.old) {
         warn.ind <- TRUE
         temp <- match(stratum.levels, stratum.levels.old)
         temp.str <- vecprint(stratum.levels.old[-temp])
         warn <- paste("The following strata no longer contain any values and were removed from the \nanalysis:\n", temp.str, sep="")
         act <- "Strata were removed from the analysis.\n"
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2], indicator=warn.vec[3],
            stratum=NA, warning=I(warn), action=I(act)))
         if(!is.null(popsize))
            popsize <- popsize[temp]
      }

      ind <- FALSE
      for(i in 1:nstrata) {
         stratum.i <- stratum == stratum.levels[i]
         if(sum(stratum.i) == 1) {
            warn.ind <- TRUE
            warn <- paste("The stratum named", stratum.levels[i], "contains a single value and was removed from the analysis.\n")
            act <- "Stratum was removed from the analysis.\n"
            warn.df <- rbind(warn.df, data.frame(func=I(fname),
               subpoptype=warn.vec[1], subpop=warn.vec[2],
               indicator=warn.vec[3], stratum=NA, warning=I(warn),
               action=I(act)))
            z <- z[!stratum.i]
            wgt <- wgt[!stratum.i]
            if(vartype == "Local") {
               x <- x[!stratum.i]
               y <- y[!stratum.i]
            }
            stratum <- stratum[!stratum.i]
            if(cluster.ind) {
               cluster <- cluster[!stratum.i]
               wgt1 <- wgt1[!stratum.i]
               if(vartype == "Local") {
                  x1 <- x1[!stratum.i]
                  y1 <- y1[!stratum.i]
               }
            }
            if(swgt.ind) {
               swgt <- swgt[!stratum.i]
               if(cluster.ind)
                  swgt1 <- swgt1[!stratum.i]
            }
            if(!is.null(popsize))
               popsize <- popsize[names(popsize) != stratum.levels[i]]
            ind <- TRUE
         }
      }
      if(ind) {
         stratum <- factor(stratum)
         stratum.levels <- levels(stratum)
         nstrata <- length(stratum.levels)
      }

# Check whether the number of strata is one

      if(nstrata == 1) {
         warn.ind <- TRUE
         warn <- "Only a single stratum was available for the analysis.\n"
         act <- "An unstratified data analysis was used.\n"
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2], indicator=warn.vec[3],
            stratum=NA, warning=I(warn), action=I(act)))
         stratum.ind <- FALSE
      }
   }

# Check whether the vector of response values is empty

   if(length(z) == 0)
      stop("\nEstimates cannot be calculated since the vector of response values is empty.")

# If the sample has two stages, determine whether there are any stage one
# sampling units with a sufficient number of sites to allow variance calculation

   if(cluster.ind) {
      temp <- sapply(split(cluster, cluster), length) == 1
      if(all(temp)) {
         stop("\nA variance estimate cannot be calculated since all of the stage one sampling \nunit(s) contain a single stage two sampling unit.")
      }
      if(any(temp)) {
         temp.str <- vecprint(names(temp)[temp])
         warn <- paste("Since the following stage one sampling units contain a single stage two \nsampling unit, a variance estimate cannot be calculated and the mean of the \nvariance estimates for stage one sampling units with two or more sites will \nbe used:\n", temp.str, sep="")
         act <- "The mean of the variance estimates will be used.\n"
         warn.df <- rbind(warn.df, data.frame(func=I(fname), subpoptype=NA,
            subpop=NA, indicator=NA, stratum=NA, warning=I(warn),
            action=I(act)))
      }
   }

# Calculate the confidence bound multiplier

   mult <- qnorm(0.5 + (conf/100)/2)

# Calculate additional required values

   if(is.null(cdfval))
      cdfval <- sort(unique(z))
   ncdfval <- length(cdfval)
   nvec <- 1:ncdfval
   npctval <- length(pctval)
   if(swgt.ind) {
      if(!is.null(popsize)) {
         sum.popsize <- sum(popsize)
      } else {
         if(stratum.ind) {
            if(cluster.ind) {
               popsize.hat <- tapply(wgt*swgt*wgt1*swgt1, stratum, sum)
               sum.popsize.hat <- sum(wgt*swgt*wgt1*swgt1)
            } else {
               popsize.hat <- tapply(wgt*swgt, stratum, sum)
               sum.popsize.hat <- sum(wgt*swgt)
            }
         } else {
            if(cluster.ind)
               popsize.hat <- sum(wgt*swgt*wgt1*swgt1)
            else
               popsize.hat <- sum(wgt*swgt)
         }
      }
   } else {
      if(!is.null(popsize)) {
         sum.popsize <- sum(popsize)
      } else {
         if(stratum.ind) {
            if(cluster.ind) {
               popsize.hat <- tapply(wgt*wgt1, stratum, sum)
               sum.popsize.hat <- sum(wgt*wgt1)
            } else {
               popsize.hat <- tapply(wgt, stratum, sum)
               sum.popsize.hat <- sum(wgt)
            }
         } else {
            if(cluster.ind)
               popsize.hat <- sum(wgt*wgt1)
            else
               popsize.hat <- sum(wgt)
         }
      }
   }

# Branch to handle stratified and unstratified data

   if(stratum.ind) {

# Begin the section for stratified data

# Create the object for the estimates

   Results <- vector("list", 2)
   names(Results) <- c("CDF", "Pct")

# Begin subsection for CDF estimates

# Create the data frame for CDF estimates

   rslt <- data.frame(array(0, c(ncdfval, 10)))
   dimnames(rslt) <- list(1:ncdfval, c("Value", "NResp", "Estimate.P",
      "StdError.P", paste("LCB", conf, "Pct.P", sep=""), paste("UCB", conf,
      "Pct.P", sep=""), "Estimate.U", "StdError.U", paste("LCB", conf, "Pct.U",
      sep=""), paste("UCB", conf, "Pct.U", sep="")))
   rslt[,1] <- cdfval

# Begin the subsection for individual strata

   for(i in 1:nstrata) {

# Calculate the CDF estimates and variance estimates

   stratum.i <- stratum == stratum.levels[i]
   if(swgt.ind) {
      cdfest.p <- cdf.size.prop(z[stratum.i], wgt[stratum.i], cdfval,
         cluster.ind, cluster[stratum.i], wgt1[stratum.i], swgt[stratum.i],
         swgt1[stratum.i])
      temp <- cdfvar.size.prop(z[stratum.i], wgt[stratum.i], x[stratum.i],
         y[stratum.i], cdfval, cdfest.p, stratum.ind, stratum.levels[i],
         cluster.ind, cluster[stratum.i], wgt1[stratum.i], x1[stratum.i],
         y1[stratum.i], pcfactor.ind, pcfsize[i], N.cluster[i], stage1size[[i]],
         support[stratum.i], swgt[stratum.i], swgt1[stratum.i], vartype,
         warn.ind, warn.df, warn.vec)
      varest.p <- temp$varest
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df

      cdfest.u <- cdf.size.total(z[stratum.i], wgt[stratum.i], cdfval,
         cluster.ind, cluster[stratum.i], wgt1[stratum.i], popsize[i],
         swgt[stratum.i], swgt1[stratum.i])
      temp <- cdfvar.size.total(z[stratum.i], wgt[stratum.i], x[stratum.i],
         y[stratum.i], cdfval, cdfest.u, stratum.ind, stratum.levels[i],
         cluster.ind, cluster[stratum.i], wgt1[stratum.i], x1[stratum.i],
         y1[stratum.i], popsize[i], pcfactor.ind, pcfsize[i], N.cluster[i],
         stage1size[[i]], support[stratum.i], swgt[stratum.i], swgt1[stratum.i],
         vartype, warn.ind, warn.df, warn.vec)
      varest.u <- temp$varest
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df
   } else {
      cdfest.p <- cdf.prop(z[stratum.i], wgt[stratum.i], cdfval, cluster.ind,
         cluster[stratum.i], wgt1[stratum.i])
      temp <- cdfvar.prop(z[stratum.i], wgt[stratum.i], x[stratum.i],
         y[stratum.i], cdfval, cdfest.p, stratum.ind, stratum.levels[i],
         cluster.ind, cluster[stratum.i], wgt1[stratum.i], x1[stratum.i],
         y1[stratum.i], pcfactor.ind, pcfsize[i], N.cluster[i], stage1size[[i]],
         support[stratum.i], vartype, warn.ind, warn.df, warn.vec)
      varest.p <- temp$varest
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df

      cdfest.u <- cdf.total(z[stratum.i], wgt[stratum.i], cdfval, cluster.ind,
         cluster[stratum.i], wgt1[stratum.i], popsize[i])
      temp <- cdfvar.total(z[stratum.i], wgt[stratum.i], x[stratum.i],
         y[stratum.i], cdfval, cdfest.u, stratum.ind, stratum.levels[i],
         cluster.ind, cluster[stratum.i], wgt1[stratum.i], x1[stratum.i],
         y1[stratum.i], popsize[i], pcfactor.ind, pcfsize[i], N.cluster[i],
         stage1size[[i]], support[stratum.i], vartype, warn.ind, warn.df,
         warn.vec)
      varest.u <- temp$varest
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df
   }

# Add estimates to the data frame for all strata combined

   if(!is.null(popsize)) {
      rslt[,3] <- rslt[,3] + (popsize[i]/sum.popsize)*cdfest.p
      rslt[,4] <- rslt[,4] + ((popsize[i]/sum.popsize)^2)*varest.p
   } else {
      rslt[,3] <- rslt[,3] + (popsize.hat[i]/sum.popsize.hat)*cdfest.p
      rslt[,4] <- rslt[,4] + ((popsize.hat[i]/sum.popsize.hat)^2)*varest.p
   }
   rslt[,7] <- rslt[,7] + cdfest.u
   rslt[,8] <- rslt[,8] + varest.u

# End the subsection for individual strata

   }

# Begin the subsection for all strata combined

# Calculate the number of response values, standard errors, and confidence
# bounds

   rslt[,2] <- cdf.nresp(z, cdfval)

   rslt[,4] <- sqrt(rslt[,4])

   rslt[,8] <- sqrt(rslt[,8])

   rslt[,5] <- 100*pmax(rslt[,3] - mult*rslt[,4], 0)
   rslt[,6] <- 100*pmin(rslt[,3] + mult*rslt[,4], 1)
   rslt[,3] <- 100*rslt[,3]
   rslt[,4] <- 100*rslt[,4]

   rslt[,9] <- pmax(rslt[,7] - mult*rslt[,8], 0)
   if(!is.null(popsize))
      rslt[,10] <- pmin(rslt[,7] + mult*rslt[,8], sum.popsize)
   else
      rslt[,10] <- rslt[,7] + mult*rslt[,8]

# Assign results to the data frame for estimates

   Results$CDF <- data.frame(rslt)

# End the subsection for all strata combined

# End subsection for CDF estimates

# Begin subsection for percentile estimates

# Create the data frame for percentile estimates

   rslt <- data.frame(array(0, c(npctval, 10)))
   dimnames(rslt) <- list(1:npctval, c("Statistic", "NResp", "Estimate.P",
      "StdError.P", paste("LCB", conf, "Pct.P", sep=""), paste("UCB", conf,
      "Pct.P", sep=""), "Estimate.U", "StdError.U", paste("LCB", conf, "Pct.U",
      sep=""), paste("UCB", conf, "Pct.U", sep="")))
   rslt[,1] <- paste(pctval, "Pct", sep="")
   rslt[,4] <- I(character(npctval))
   rslt[,8] <- I(character(npctval))

# Convert the input percentile values to proportions and to percentage of total,
# as appropriate

   pctval.p <- pctval/100
   if(!is.null(popsize))
      pctval.u <- (pctval/100)*sum.popsize
   else
      pctval.u <- (pctval/100)*sum.popsize.hat

# Determine whether all response values are equal and assign estimates when true

   if(min(z) == max(z)) {
      rslt[,2] <- length(z)
      rslt[,c(3,5,6,7,9,10)] <- max(z)

   } else {

# Calculate percentile estimates

      cdfest.p <- Results$CDF$Estimate.P/100
      cdfest.u <- Results$CDF$Estimate.U
      for(j in 1:npctval) {
         high <- ifelse(length(nvec[cdfest.p >= pctval.p[j]]) > 0,
            min(nvec[cdfest.p >= pctval.p[j]]), NA)
         low <- ifelse(length(nvec[cdfest.p <= pctval.p[j]]) > 0,
            max(nvec[cdfest.p <= pctval.p[j]]), NA)
         if(is.na(high)) {
            rslt[j,3] <- NA
         } else if(is.na(low)) {
            rslt[j,3] <- cdfval[high]
         } else {
            if(high > low)
               ival <- (pctval.p[j] - cdfest.p[low]) / (cdfest.p[high] -
                  cdfest.p[low])
            else
               ival <- 1
            rslt[j,3] <- ival*cdfval[high] + (1-ival)*cdfval[low]
         }

         high <- ifelse(length(nvec[cdfest.u >= pctval.u[j]]) > 0,
            min(nvec[cdfest.u >= pctval.u[j]]), NA)
         low <- ifelse(length(nvec[cdfest.u <= pctval.u[j]]) > 0,
            max(nvec[cdfest.u <= pctval.u[j]]), NA)
         if(is.na(high)) {
            rslt[j,7] <- NA
         } else if(is.na(low)) {
            rslt[j,7] <- cdfval[high]
         } else {
            if(high > low)
               ival <- (pctval.u[j] - cdfest.u[low]) / (cdfest.u[high] -
                  cdfest.u[low])
            else
               ival <- 1
            rslt[j,7] <- ival*cdfval[high] + (1-ival)*cdfval[low]
         }
      }

# Determine number of response values

      rslt[,2] <- cdf.nresp(z, rslt[,3])

# Calculate variance of the inverse percentile estimates

      temp.p <- !is.na(rslt[,3])
      varest.p <- ifelse(temp.p, 0, NA)
      lbound.p <- rep(NA, npctval)
      ubound.p <- rep(NA, npctval)
      temp.u <- !is.na(rslt[,7])
      varest.u <- ifelse(temp.u, 0, NA)
      lbound.u <- rep(NA, npctval)
      ubound.u <- rep(NA, npctval)
      for(i in 1:nstrata) {
         stratum.i <- stratum == stratum.levels[i]
         if(swgt.ind) {
            if(!is.null(popsize)) {
               temp <- cdfvar.size.prop(z[stratum.i], wgt[stratum.i],
                  x[stratum.i], y[stratum.i], rslt[temp.p,3], pctval.p[temp.p],
                  stratum.ind, stratum.levels[i], cluster.ind,
                  cluster[stratum.i], wgt1[stratum.i], x1[stratum.i],
                  y1[stratum.i], pcfactor.ind, pcfsize[i], N.cluster[i],
                  stage1size[[i]], support[stratum.i], swgt[stratum.i],
                  swgt1[stratum.i], vartype, warn.ind, warn.df, warn.vec)
               varest.p[temp.p] <- varest.p[temp.p] +
                  ((popsize[i]/sum.popsize)^2)*temp$varest
               warn.ind <- temp$warn.ind
               warn.df <- temp$warn.df
            } else {
               temp <- cdfvar.size.prop(z[stratum.i], wgt[stratum.i],
                  x[stratum.i], y[stratum.i], rslt[temp.p,3], pctval.p[temp.p],
                  stratum.ind, stratum.levels[i], cluster.ind,
                  cluster[stratum.i], wgt1[stratum.i], x1[stratum.i],
                  y1[stratum.i], pcfactor.ind, pcfsize[i], N.cluster[i],
                  stage1size[[i]], support[stratum.i], swgt[stratum.i],
                  swgt1[stratum.i], vartype, warn.ind, warn.df, warn.vec)
               varest.p[temp.p] <- varest.p[temp.p] + ((popsize.hat[i]/
                  sum.popsize.hat)^2)*temp$varest
               warn.ind <- temp$warn.ind
               warn.df <- temp$warn.df
            }
            temp <- cdfvar.size.total(z[stratum.i], wgt[stratum.i],
               x[stratum.i], y[stratum.i], rslt[temp.u,7], pctval.u[temp.u],
               stratum.ind, stratum.levels[i], cluster.ind, cluster[stratum.i],
               wgt1[stratum.i], x1[stratum.i], y1[stratum.i], popsize[i],
               pcfactor.ind, pcfsize[i], N.cluster[i], stage1size[[i]],
               support[stratum.i], swgt[stratum.i], swgt1[stratum.i], vartype,
               warn.ind, warn.df, warn.vec)
            varest.u[temp.u] <- varest.u[temp.u] + temp$varest
            warn.ind <- temp$warn.ind
            warn.df <- temp$warn.df
         } else {
            if(!is.null(popsize)) {
               temp <- cdfvar.prop(z[stratum.i], wgt[stratum.i], x[stratum.i],
                  y[stratum.i], rslt[temp.p,3], pctval.p[temp.p], stratum.ind,
                  stratum.levels[i], cluster.ind, cluster[stratum.i],
                  wgt1[stratum.i], x1[stratum.i], y1[stratum.i], pcfactor.ind,
                  pcfsize[i], N.cluster[i], stage1size[[i]], support[stratum.i],
                  vartype, warn.ind, warn.df, warn.vec)
               varest.p[temp.p] <- varest.p[temp.p] +
                  ((popsize[i]/sum.popsize)^2)*temp$varest
               warn.ind <- temp$warn.ind
               warn.df <- temp$warn.df
            } else {
               temp <- cdfvar.prop(z[stratum.i], wgt[stratum.i], x[stratum.i],
                  y[stratum.i], rslt[temp.p,3], pctval.p[temp.p], stratum.ind,
                  stratum.levels[i], cluster.ind, cluster[stratum.i],
                  wgt1[stratum.i], x1[stratum.i], y1[stratum.i], pcfactor.ind,
                  pcfsize[i], N.cluster[i], stage1size[[i]], support[stratum.i],
                  vartype, warn.ind, warn.df, warn.vec)
               varest.p[temp.p] <- varest.p[temp.p] + ((popsize.hat[i]/
                  sum.popsize.hat)^2)*temp$varest
               warn.ind <- temp$warn.ind
               warn.df <- temp$warn.df
            }
            temp <- cdfvar.total(z[stratum.i], wgt[stratum.i], x[stratum.i],
               y[stratum.i], rslt[temp.u,7], pctval.u[temp.u], stratum.ind,
               stratum.levels[i], cluster.ind, cluster[stratum.i],
               wgt1[stratum.i], x1[stratum.i], y1[stratum.i], popsize[i],
               pcfactor.ind, pcfsize[i], N.cluster[i], stage1size[[i]],
               support[stratum.i], vartype, warn.ind, warn.df, warn.vec)
            varest.u[temp.u] <- varest.u[temp.u] + temp$varest
            warn.ind <- temp$warn.ind
            warn.df <- temp$warn.df
         }
      }

# Calculate confidence bounds of the inverse percentile estimates

      lbound.p[temp.p] <- pmax(pctval.p[temp.p] - mult*sqrt(varest.p[temp.p]),
         0)
      ubound.p[temp.p] <- pmin(pctval.p[temp.p] + mult*sqrt(varest.p[temp.p]),
         1)
      lbound.u[temp.u] <- pmax(pctval.u[temp.u] - mult*sqrt(varest.u[temp.u]),
         0)
      if(!is.null(popsize))
         ubound.u[temp.u] <- pmin(pctval.u[temp.u] +
            mult*sqrt(varest.u[temp.u]), sum.popsize)
      else
         ubound.u[temp.u] <- pctval.u[temp.u] + mult*sqrt(varest.u[temp.u])

# Calculate confidence bounds of the percentile estimates

      for(j in 1:npctval) {
         if(temp.p[j]) {
            high <- ifelse(length(nvec[cdfest.p >= lbound.p[j]]) > 0,
               min(nvec[cdfest.p >= lbound.p[j]]), NA)
            low <- ifelse(length(nvec[cdfest.p <= lbound.p[j]]) > 0,
               max(nvec[cdfest.p <= lbound.p[j]]), NA)
            if(is.na(high)) {
               rslt[j,5] <- NA
            } else if(is.na(low)) {
               rslt[j,5] <- cdfval[high]
            } else {
               if(high > low)
                  ival <- (lbound.p[j] - cdfest.p[low]) / (cdfest.p[high] -
                     cdfest.p[low])
               else
                  ival <- 1
               rslt[j,5] <- ival*cdfval[high] + (1-ival)*cdfval[low]
            }

            high <- ifelse(length(nvec[cdfest.p >= ubound.p[j]]) > 0,
               min(nvec[cdfest.p >= ubound.p[j]]), NA)
            low <- ifelse(length(nvec[cdfest.p <= ubound.p[j]]) > 0,
               max(nvec[cdfest.p <= ubound.p[j]]), NA)
            if(is.na(high)) {
               rslt[j,6] <- NA
            } else if(is.na(low)) {
               rslt[j,6] <- cdfval[high]
            } else {
               if(high > low)
                  ival <- (ubound.p[j] - cdfest.p[low]) / (cdfest.p[high] -
                     cdfest.p[low])
               else
                  ival <- 1
               rslt[j,6] <- ival*cdfval[high] + (1-ival)*cdfval[low]
            }
         } else {
            rslt[j,5] <- NA
            rslt[j,6] <- NA
         }

         if(temp.u[j]) {
            high <- ifelse(length(nvec[cdfest.u >= lbound.u[j]]) > 0,
               min(nvec[cdfest.u >= lbound.u[j]]), NA)
            low <- ifelse(length(nvec[cdfest.u <= lbound.u[j]]) > 0,
               max(nvec[cdfest.u <= lbound.u[j]]), NA)
            if(is.na(high)) {
               rslt[j,9] <- NA
            } else if(is.na(low)) {
               rslt[j,9] <- cdfval[high]
            } else {
               if(high > low)
                  ival <- (lbound.u[j] - cdfest.u[low]) / (cdfest.u[high] -
                     cdfest.u[low])
               else
                  ival <- 1
               rslt[j,9] <- ival*cdfval[high] + (1-ival)*cdfval[low]
            }

            high <- ifelse(length(nvec[cdfest.u >= ubound.u[j]]) > 0,
               min(nvec[cdfest.u >= ubound.u[j]]), NA)
            low <- ifelse(length(nvec[cdfest.u <= ubound.u[j]]) > 0,
               max(nvec[cdfest.u <= ubound.u[j]]), NA)
            if(is.na(high)) {
               rslt[j,10] <- NA
            } else if(is.na(low)) {
               rslt[j,10] <- cdfval[high]
            } else {
               if(high > low)
                  ival <- (ubound.u[j] - cdfest.u[low]) / (cdfest.u[high] -
                     cdfest.u[low])
               else
                  ival <- 1
               rslt[j,10] <- ival*cdfval[high] + (1-ival)*cdfval[low]
            }
         } else {
            rslt[j,9] <- NA
            rslt[j,10] <- NA
         }
      }
   }

# Assign results to the data frame for estimates

   Results$Pct <- rslt

# End the subsection for percentile estimates

# End the section for stratified data

   } else {

# Begin the section for unstratified data

# Check whether the vector of response values contains a single element

   if(length(z) == 1)
      stop("\nEstimates cannot be calculated since the vector of response values contains a \nsingle element.")

# Begin subsection for CDF estimates

# Create the object for the estimates

   Results <- vector("list", 2)
   names(Results) <- c("CDF", "Pct")

# Calculate the number of response values, CDF estimates, standard error
# estimates, and confidence bounds

   nresp <- cdf.nresp(z, cdfval)

   if(swgt.ind) {
      cdfest.p <- cdf.size.prop(z, wgt, cdfval, cluster.ind, cluster, wgt1,
         swgt, swgt1)
      temp <- cdfvar.size.prop(z, wgt, x, y, cdfval, cdfest.p, stratum.ind,
         NULL, cluster.ind, cluster, wgt1, x1, y1, pcfactor.ind, pcfsize,
         N.cluster, stage1size, support, swgt, swgt1, vartype, warn.ind,
         warn.df, warn.vec)
      sdest.p <- sqrt(temp$varest)
      lbound.p <- pmax(cdfest.p - mult*sdest.p, 0)
      ubound.p <- pmin(cdfest.p + mult*sdest.p, 1)
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df

      cdfest.u <- cdf.size.total(z, wgt, cdfval, cluster.ind, cluster, wgt1,
         popsize, swgt, swgt1)
      temp <- cdfvar.size.total(z, wgt, x, y, cdfval, cdfest.u, stratum.ind,
         NULL, cluster.ind, cluster, wgt1, x1, y1, popsize, pcfactor.ind,
         pcfsize, N.cluster, stage1size, support, swgt, swgt1, vartype,
         warn.ind, warn.df, warn.vec)
      sdest.u <- sqrt(temp$varest)
      lbound.u <- pmax(cdfest.u - mult*sdest.u, 0)
      if(!is.null(popsize))
         ubound.u <- pmin(cdfest.u + mult*sdest.u, popsize)
      else
         ubound.u <- cdfest.u + mult*sdest.u
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df
   } else {
      cdfest.p <- cdf.prop(z, wgt, cdfval, cluster.ind, cluster, wgt1)
      temp <- cdfvar.prop(z, wgt, x, y, cdfval, cdfest.p, stratum.ind, NULL,
         cluster.ind, cluster, wgt1, x1, y1, pcfactor.ind, pcfsize, N.cluster,
         stage1size, support, vartype, warn.ind, warn.df, warn.vec)
      sdest.p <- sqrt(temp$varest)
      lbound.p <- pmax(cdfest.p - mult*sdest.p, 0)
      ubound.p <- pmin(cdfest.p + mult*sdest.p, 1)
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df

      cdfest.u <- cdf.total(z, wgt, cdfval, cluster.ind, cluster, wgt1, popsize)
      temp <- cdfvar.total(z, wgt, x, y, cdfval, cdfest.u, stratum.ind, NULL,
         cluster.ind, cluster, wgt1, x1, y1, popsize, pcfactor.ind, pcfsize,
         N.cluster, stage1size, support, vartype, warn.ind, warn.df, warn.vec)
      sdest.u <- sqrt(temp$varest)
      lbound.u <- pmax(cdfest.u - mult*sdest.u, 0)
      if(!is.null(popsize))
         ubound.u <- pmin(cdfest.u + mult*sdest.u, popsize)
      else
         ubound.u <- cdfest.u + mult*sdest.u
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df
   }

   rslt <- cbind(cdfval, nresp, 100*cdfest.p, 100*sdest.p, 100*lbound.p,
      100*ubound.p, cdfest.u, sdest.u, lbound.u, ubound.u)
   dimnames(rslt) <- list(1:ncdfval, c("Value", "NResp", "Estimate.P",
      "StdError.P", paste("LCB", conf, "Pct.P", sep=""), paste("UCB", conf,
      "Pct.P", sep=""), "Estimate.U", "StdError.U", paste("LCB", conf, "Pct.U",
      sep=""), paste("UCB", conf, "Pct.U", sep="")))

# Assign results to the data frame for estimates

   Results$CDF <- data.frame(rslt)

# End subsection for CDF estimates

# Begin subsection for percentile estimates

# Create the data frame for percentile estimates

   rslt <- data.frame(array(0, c(npctval, 10)))
   dimnames(rslt) <- list(1:npctval, c("Statistic", "NResp", "Estimate.P",
      "StdError.P", paste("LCB", conf, "Pct.P", sep=""), paste("UCB", conf, "Pct.P",
      sep=""), "Estimate.U", "StdError.U", paste("LCB", conf, "Pct.U", sep=""),
      paste("UCB", conf, "Pct.U", sep="")))
   rslt[,1] <- paste(pctval, "Pct", sep="")
   rslt[,4] <- I(character(npctval))
   rslt[,8] <- I(character(npctval))


# Convert the input percentile values to proportions or to percentage of total,
# as appropriate

   pctval.p <- pctval/100
   if(!is.null(popsize))
      pctval.u <- (pctval/100)*popsize
   else
      pctval.u <- (pctval/100)*popsize.hat

# Determine whether all response values are equal and assign estimates when true

   if(min(z) == max(z)) {
      rslt[,2] <- length(z)
      rslt[,c(3,5,6,7,9,10)] <- max(z)

   } else {

# Calculate percentile estimates

      cdfest.p <- Results$CDF$Estimate.P/100
      cdfest.u <- Results$CDF$Estimate.U
      for(j in 1:npctval) {
         high <- ifelse(length(nvec[cdfest.p >= pctval.p[j]]) > 0,
            min(nvec[cdfest.p >= pctval.p[j]]), NA)
         low <- ifelse(length(nvec[cdfest.p <= pctval.p[j]]) > 0,
            max(nvec[cdfest.p <= pctval.p[j]]), NA)
         if(is.na(high)) {
            rslt[j,3] <- NA
         } else if(is.na(low)) {
            rslt[j,3] <- cdfval[high]
         } else {
            if(high > low)
               ival <- (pctval.p[j] - cdfest.p[low]) / (cdfest.p[high] -
                  cdfest.p[low])
            else
               ival <- 1
            rslt[j,3] <- ival*cdfval[high] + (1-ival)*cdfval[low]
         }

         high <- ifelse(length(nvec[cdfest.u >= pctval.u[j]]) > 0,
            min(nvec[cdfest.u >= pctval.u[j]]), NA)
         low <- ifelse(length(nvec[cdfest.u <= pctval.u[j]]) > 0,
            max(nvec[cdfest.u <= pctval.u[j]]), NA)
         if(is.na(high)) {
            rslt[j,7] <- NA
         } else if(is.na(low)) {
            rslt[j,7] <- cdfval[high]
         } else {
            if(high > low)
               ival <- (pctval.u[j] - cdfest.u[low]) / (cdfest.u[high] -
                  cdfest.u[low])
            else
               ival <- 1
            rslt[j,7] <- ival*cdfval[high] + (1-ival)*cdfval[low]
         }
      }

# Determine number of response values

      rslt[,2] <- cdf.nresp(z, rslt[,3])

# Calculate confidence bounds of the inverse percentile estimates

      temp.p <- !is.na(rslt[,3])
      sdest.p <- rep(NA, npctval)
      lbound.p <- rep(NA, npctval)
      ubound.p <- rep(NA, npctval)
      temp.u <- !is.na(rslt[,7])
      sdest.u <- rep(NA, npctval)
      lbound.u <- rep(NA, npctval)
      ubound.u <- rep(NA, npctval)
      if(swgt.ind) {
         temp <- cdfvar.size.prop(z, wgt, x, y, rslt[temp.p,3],
            pctval.p[temp.p], stratum.ind, NULL, cluster.ind, cluster, wgt1, x1,
            y1, pcfactor.ind, pcfsize, N.cluster, stage1size, support, swgt,
            swgt1, vartype, warn.ind, warn.df, warn.vec)
         sdest.p[temp.p] <- sqrt(temp$varest)
         lbound.p[temp.p] <- pmax(pctval.p[temp.p] - mult*sdest.p[temp.p], 0)
         ubound.p[temp.p] <- pmin(pctval.p[temp.p] + mult*sdest.p[temp.p], 1)
         warn.ind <- temp$warn.ind
         warn.df <- temp$warn.df

         if(!is.null(popsize)) {
            temp <- cdfvar.size.total(z, wgt, x, y, rslt[temp.u,7],
               pctval.u[temp.u], stratum.ind, NULL, cluster.ind, cluster, wgt1,
               x1, y1, popsize, pcfactor.ind, pcfsize, N.cluster, stage1size,
               support, swgt, swgt1, vartype, warn.ind, warn.df, warn.vec)
            sdest.u[temp.u] <- sqrt(temp$varest)
            ubound.u[temp.u] <- pmin(pctval.u[temp.u] + mult*sdest.u[temp.u],
               popsize)
            warn.ind <- temp$warn.ind
            warn.df <- temp$warn.df
         } else {
            temp <- cdfvar.size.total(z, wgt, x, y, rslt[temp.u,7],
               pctval.u[temp.u], stratum.ind, NULL, cluster.ind, cluster, wgt1,
               x1, y1, popsize, pcfactor.ind, pcfsize, N.cluster, stage1size,
               support, swgt, swgt1, vartype, warn.ind, warn.df, warn.vec)
            sdest.u[temp.u] <- sqrt(temp$varest)
            ubound.u[temp.u] <- pctval.u[temp.u] + mult*sdest.u[temp.u]
            warn.ind <- temp$warn.ind
            warn.df <- temp$warn.df
         }
         lbound.u[temp.u] <- pmax(pctval.u[temp.u] - mult*sdest.u[temp.u], 0)
      } else {
         temp <- cdfvar.prop(z, wgt, x, y, rslt[temp.p,3], pctval.p[temp.p],
            stratum.ind, NULL, cluster.ind, cluster, wgt1, x1, y1, pcfactor.ind,
            pcfsize, N.cluster, stage1size, support, vartype, warn.ind, warn.df,
            warn.vec)
         sdest.p[temp.p] <- sqrt(temp$varest)
         lbound.p[temp.p] <- pmax(pctval.p[temp.p] - mult*sdest.p[temp.p], 0)
         ubound.p[temp.p] <- pmin(pctval.p[temp.p] + mult*sdest.p[temp.p], 1)
         warn.ind <- temp$warn.ind
         warn.df <- temp$warn.df

         if(!is.null(popsize)) {
            temp <- cdfvar.total(z, wgt, x, y, rslt[temp.u,7], pctval.u[temp.u],
               stratum.ind, NULL, cluster.ind, cluster, wgt1, x1, y1, popsize,
               pcfactor.ind, pcfsize, N.cluster, stage1size, support, vartype,
               warn.ind, warn.df, warn.vec)
            sdest.u[temp.u] <- sqrt(temp$varest)
            ubound.u[temp.u] <- pmin(pctval.u[temp.u] + mult*sdest.u[temp.u],
               popsize)
            warn.ind <- temp$warn.ind
            warn.df <- temp$warn.df
         } else {
            temp <- cdfvar.total(z, wgt, x, y, rslt[temp.u,7], pctval.u[temp.u],
               stratum.ind, NULL, cluster.ind, cluster, wgt1, x1, y1, popsize,
               pcfactor.ind, pcfsize, N.cluster, stage1size, support, vartype,
               warn.ind, warn.df, warn.vec)
            sdest.u[temp.u] <- sqrt(temp$varest)
            ubound.u[temp.u] <- pctval.u[temp.u] + mult*sdest.u[temp.u]
            warn.ind <- temp$warn.ind
            warn.df <- temp$warn.df
         }
         lbound.u[temp.u] <- pmax(pctval.u[temp.u] - mult*sdest.u[temp.u], 0)
      }

# Calculate confidence bounds of the percentile estimates

      for(j in 1:npctval) {
         high <- ifelse(length(nvec[cdfest.p >= lbound.p[j]]) > 0,
            min(nvec[cdfest.p >= lbound.p[j]]), NA)
         low <- ifelse(length(nvec[cdfest.p <= lbound.p[j]]) > 0,
            max(nvec[cdfest.p <= lbound.p[j]]), NA)
         if(is.na(high)) {
            rslt[j,5] <- NA
         } else if(is.na(low)) {
            rslt[j,5] <- cdfval[high]
         } else {
            if(high > low)
               ival <- (lbound.p[j] - cdfest.p[low]) / (cdfest.p[high] -
                  cdfest.p[low])
            else
               ival <- 1
            rslt[j,5] <- ival*cdfval[high] + (1-ival)*cdfval[low]
         }

         high <- ifelse(length(nvec[cdfest.p >= ubound.p[j]]) > 0,
            min(nvec[cdfest.p >= ubound.p[j]]), NA)
         low <- ifelse(length(nvec[cdfest.p <= ubound.p[j]]) > 0,
            max(nvec[cdfest.p <= ubound.p[j]]), NA)
         if(is.na(high)) {
            rslt[j,6] <- NA
         } else if(is.na(low)) {
            rslt[j,6] <- cdfval[high]
         } else {
            if(high > low)
               ival <- (ubound.p[j] - cdfest.p[low]) / (cdfest.p[high] -
                  cdfest.p[low])
            else
               ival <- 1
            rslt[j,6] <- ival*cdfval[high] + (1-ival)*cdfval[low]
         }

         high <- ifelse(length(nvec[cdfest.u >= lbound.u[j]]) > 0,
            min(nvec[cdfest.u >= lbound.u[j]]), NA)
         low <- ifelse(length(nvec[cdfest.u <= lbound.u[j]]) > 0,
            max(nvec[cdfest.u <= lbound.u[j]]), NA)
         if(is.na(high)) {
            rslt[j,9] <- NA
         } else if(is.na(low)) {
            rslt[j,9] <- cdfval[high]
         } else {
            if(high > low)
               ival <- (lbound.u[j] - cdfest.u[low]) / (cdfest.u[high] -
                  cdfest.u[low])
            else
               ival <- 1
            rslt[j,9] <- ival*cdfval[high] + (1-ival)*cdfval[low]
         }

         high <- ifelse(length(nvec[cdfest.u >= ubound.u[j]]) > 0,
            min(nvec[cdfest.u >= ubound.u[j]]), NA)
         low <- ifelse(length(nvec[cdfest.u <= ubound.u[j]]) > 0,
            max(nvec[cdfest.u <= ubound.u[j]]), NA)
         if(is.na(high)) {
            rslt[j,10] <- NA
         } else if(is.na(low)) {
            rslt[j,10] <- cdfval[high]
         } else {
            if(high > low)
               ival <- (ubound.u[j] - cdfest.u[low]) / (cdfest.u[high] -
                  cdfest.u[low])
            else
               ival <- 1
            rslt[j,10] <- ival*cdfval[high] + (1-ival)*cdfval[low]
         }
      }
   }

# Assign results to the data frame for estimates

   Results$Pct <- rslt

# End subsection for percentile estimates

# End section for unstratified data

   }

# Depending on whether the function was called directly or was called by
# cont.analysis, return appropriate results

   if(is.na(warn.vec[1])) {

# As necessary, output a message indicating that warning messages were generated
# during execution of the program

      if(warn.ind) {
         warn.df <<- warn.df
         if(nrow(warn.df) == 1)
            cat("During execution of the program, a warning message was generated.  The warning \nmessage is stored in a data frame named 'warn.df'.  Enter the following command \nto view the warning message: warnprnt()\n")
         else
            cat(paste("During execution of the program,", nrow(warn.df), "warning messages were generated.  The warning \nmessages are stored in a data frame named 'warn.df'.  Enter the following \ncommand to view the warning messages: warnprnt() \nTo view a subset of the warning messages (say, messages number 1, 3, and 5), \nenter the following command: warnprnt(m=c(1,3,5))\n"))
      }

# Return the Results data frame

      Results

   } else {

# Return the Results data frame, the warn.ind logical value, and the warn.df
# data frame

      names(Results$Pct)[3:6] <- substr(names(Results$Pct)[3:6], 1,
         nchar(names(Results$Pct)[3:6])-2)
      list(Results=Results, warn.ind=warn.ind, warn.df=warn.df)
   }
}
