total.est <- function(z, wgt, x=NULL, y=NULL, stratum=NULL, cluster=NULL,
   wgt1=NULL, x1=NULL, y1=NULL, popsize=NULL, popcorrect=FALSE, pcfsize=NULL,
   N.cluster=NULL, stage1size=NULL, support=NULL, sizeweight=FALSE, swgt=NULL,
   swgt1=NULL, vartype="Local", conf=95, check.ind=TRUE, warn.ind=NULL,
   warn.df=NULL, warn.vec=NULL) {

################################################################################
# Function: total.est
# Programmer: Tom Kincaid
# Date: December 18, 2000
# Last Revised: April 6, 2011
# Description:
#   This function calculates estimates of the population total, mean, variance, 
#   and standard deviation of a response variable, where the response variable
#   may be defined for either a finite or an extensive resource.  In addition 
#   the standard error of the population estimates and confidence bounds are 
#   calculated.  The Horvitz-Thompson estimator is used to calculate the 
#   total, variance, and standard deviation estimates.  The Horvitz-Thompson 
#   ratio estimator, i.e., the ratio of two Horvitz-Thompson estimators, is used 
#   to calculate the mean estimate.  Variance estimates are calculated using 
#   either the local mean variance estimator or the simple random sampling (SRS) 
#   variance estimator.  The choice of variance estimator is subject to user 
#   control.  The local mean variance estimator requires the x-coordinate and the
#   y-coordinate of each site.  The SRS variance estimator uses the independent 
#   random sample approximation to calculate joint inclusion probabilities.  
#   Confidence bounds are calculated using a Normal distribution multiplier.
#   The function can accommodate a stratified sample.  For a stratified sample, 
#   separate estimates and standard errors are calculated for each stratum, which
#   are used to produce estimates and standard errors for all strata combined.  
#   Strata that contain a single value are removed.  For a stratified sample, 
#   when either the size of the resource or the sum of the size-weights of the 
#   resource is provided for each stratum, those values are used as stratum 
#   weights for calculating the estimates and standard errors for all strata 
#   combined.  For a stratified sample when neither the size of the resource nor 
#   the sum of the size-weights of the resource is provided for each stratum, 
#   estimated values are used as stratum weights for calculating the estimates 
#   and standard errors for all strata combined.  The function can accommodate 
#   single-stage and two-stage samples for both stratified and unstratified 
#   sampling designs.  Finite population and continuous population correction 
#   factors can be utilized in variance estimation.  The function checks for 
#   compatibility of input values and removes missing values.
# Arguments:
#   z = the response value for each site.
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
#   popsize = known size of the resource, which is used to calculate strata
#     proportions for calculating estimates for a stratified sample.  For a
#     finite resource, this argument is either the total number of sampling
#     units or the known sum of size-weights.  For an extensive resource, this
#     argument is the measure of the resource, i.e., either known total length
#     for a linear resource or known total area for an areal resource.  For a
#     stratified sample this variable must be a vector containing a value for
#     each stratum and must have the names attribute set to identify the stratum
#     codes.  The default is NULL.
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
#   object in list format composed of the Results data frame, which contains
#   estimates and confidence bounds, and the warn.df data frame, which contains
#   warning messages.  If the function was called directly, then output is the
#   Results data frame.
# Other Functions Required:
#   input.check - check input values for errors, consistency, and compatibility
#     with analytical functions
#   wnas - remove missing values
#   vecprint - takes an input vector and outputs a character string with line
#     breaks inserted
#   total.var - calculate variance of the total, mean, variance, and standard
#     deviation estimates
# Examples:
#   z <- rnorm(100, 10, 1)
#   wgt <- runif(100, 10, 100)
#   total.est(z, wgt, vartype="SRS")
#
#   x <- runif(100)
#   y <- runif(100)
#   total.est(z, wgt, x, y)
################################################################################

# As necessary, create a data frame for warning messages

   if(is.null(warn.ind)) {
      warn.ind <- FALSE
      warn.df <- NULL
      warn.vec <- rep(NA, 3)
   }
   fname <- "total.est"

# Check for existence of the response variable and determine the number of values

   if(is.null(z))
      stop("\nValues for the response variable must be provided.")
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
         stage1size, support, swgt.ind, swgt, swgt1, vartype, conf)

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

# Calculate confidence bound multiplier

   mult <- qnorm(0.5 + (conf/100)/2)

# Calculate additional required values

   if(swgt.ind) {
      if(!is.null(popsize))
         sum.popsize <- sum(popsize)
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
   } else {
      if(!is.null(popsize))
         sum.popsize <- sum(popsize)
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

# Branch to handle stratified and unstratified data

   if(stratum.ind) {

# Begin the section for stratified data

# Create the data frame for estimates for all strata combined

   Results <- data.frame(array(0, c(4, 6)))
   dimnames(Results) <- list(1:4, c("Statistic", "NResp", "Estimate",
      "StdError", paste("LCB", conf, "Pct", sep=""), paste("UCB", conf, "Pct",
      sep="")))
   Results[,1] <- c("Total", "Mean", "Variance", "Std. Deviation")
   Results[,2] <- rep(length(z), 4)

# Begin the subsection for individual strata

   for(i in 1:nstrata) {

# Calculate required values

      stratum.i <- stratum == stratum.levels[i]
      z.st <- z[stratum.i]
      if(swgt.ind) {
         if(cluster.ind) {
            w2 <- wgt[stratum.i]*swgt[stratum.i]
            w1 <- wgt1[stratum.i]*swgt1[stratum.i]
         } else {
            w <- wgt[stratum.i]*swgt[stratum.i]
         }
      } else {
         if(cluster.ind) {
            w2 <- wgt[stratum.i]
            w1 <- wgt1[stratum.i]
         } else {
            w <- wgt[stratum.i]
         }
      }

# Calculate the total estimate

      if(cluster.ind)
         total.est <- sum(w2*w1*z.st)
      else
         total.est <- sum(w*z.st)

# Calculate the mean estimate

      mean.est <- total.est / popsize.hat[i]

# Calculate the variance and standard deviation estimates

      if(cluster.ind)
         var.est <- sum(w2*w1*((z.st - mean.est)^2)) / (popsize.hat[i])
      else
         var.est <- sum(w*((z.st - mean.est)^2)) / (popsize.hat[i])
      sd.est <- sqrt(var.est)

# Calculate variance estimates for the estimated quantities

      if(cluster.ind)
         temp <- total.var(z.st, w2, x[stratum.i], y[stratum.i], mean.est,
            var.est, sd.est, stratum.ind, stratum.levels[i], cluster.ind,
            cluster[stratum.i], w1, x1[stratum.i], y1[stratum.i], pcfactor.ind,
            pcfsize[i], N.cluster[i], stage1size[[i]], support[stratum.i],
            vartype, warn.ind, warn.df, warn.vec)
      else
         temp <- total.var(z.st, w, x[stratum.i], y[stratum.i], mean.est,
            var.est, sd.est, stratum.ind, stratum.levels[i], cluster.ind,
            pcfactor.ind=pcfactor.ind, pcfsize=pcfsize[i],
            support=support[stratum.i], vartype=vartype, warn.ind=warn.ind,
            warn.df=warn.df, warn.vec=warn.vec)
      total.var <- temp$varest[1]
      mean.var <- temp$varest[2]
      var.var <- temp$varest[3]
      sd.var <- temp$varest[4]
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df

# Add estimates to the data frame for all strata combined

      Results[1, 3] <- Results[1, 3] + total.est
      Results[1, 4] <- Results[1, 4] + total.var
      if(!is.null(popsize)) {
         Results[2:4, 3] <- Results[2:4, 3] + (popsize[i]/sum.popsize) *
            c(mean.est, var.est, sd.est)
         Results[2:4, 4] <- Results[2:4, 4] + ((popsize[i]/sum.popsize)^2) *
            c(mean.var, var.var, sd.var)
      } else {
         Results[2:4, 3] <- Results[2:4, 3] + (popsize.hat[i]/sum.popsize.hat) *
            c(mean.est, var.est, sd.est)
         Results[2:4, 4] <- Results[2:4, 4] +
            ((popsize.hat[i]/sum.popsize.hat)^2) * c(mean.var, var.var, sd.var)
      }

# End the subsection for individual strata

   }

# Begin the subsection for all strata combined

# Calculate standard errors and confidence bounds

   Results[, 4] <- sqrt(Results[, 4])
   Results[, 5] <- Results[, 3] - mult*Results[, 4]
   Results[3:4, 5] <- pmax(Results[3:4, 5], 0)
   Results[, 6] <- Results[, 3] + mult*Results[, 4]

# End the subsection for all strata combined

# End the section for stratified data

   } else {

# Begin the section for unstratified data

# Check whether the vector of response values contains a single element

   if(length(z) == 1)
      stop("\nEstimates cannot be calculated since the vector of response values contains a \nsingle element.")

# Create the data frame for estimates

   Results <- data.frame(array(0, c(4, 6)))
   dimnames(Results) <- list(1:4, c("Statistic", "NResp", "Estimate",
      "StdError", paste("LCB", conf, "Pct", sep=""), paste("UCB", conf, "Pct",
      sep="")))
   Results[,1] <- c("Total", "Mean", "Variance", "Std. Deviation")
   Results[,2] <- rep(length(z), 4)

# Calculate required values

   if(swgt.ind) {
      if(cluster.ind) {
         w2 <- wgt*swgt
         w1 <- wgt1*swgt1
      } else {
         w <- wgt*swgt
      }
   } else {
      if(cluster.ind) {
         w2 <- wgt
         w1 <- wgt1
      } else {
         w <- wgt
      }
   }

# Calculate the total estimate

   if(cluster.ind)
      total.est <- sum(w2*w1*z)
   else
      total.est <- sum(w*z)

# Calculate the mean estimate

   mean.est <- total.est / popsize.hat

# Calculate the variance and standard deviation estimates

   if(cluster.ind)
      var.est <- sum(w2*w1*((z - mean.est)^2)) / (popsize.hat)
    else
      var.est <- sum(w*((z - mean.est)^2)) / (popsize.hat)
   sd.est <- sqrt(var.est)

# Calculate standard error estimates for the estimated quantities

      if(cluster.ind)
         temp <- total.var(z, w2, x, y, mean.est, var.est, sd.est, stratum.ind,
            NULL, cluster.ind, cluster, w1, x1, y1, pcfactor.ind, pcfsize,
            N.cluster, stage1size, support, vartype, warn.ind, warn.df,
            warn.vec)
      else
         temp <- total.var(z, w, x, y, mean.est, var.est, sd.est, stratum.ind,
             NULL, cluster.ind, pcfactor.ind=pcfactor.ind, pcfsize=pcfsize,
             support=support, vartype=vartype, warn.ind=warn.ind,
             warn.df=warn.df, warn.vec=warn.vec)
      total.sd <- sqrt(temp$varest[1])
      mean.sd <- sqrt(temp$varest[2])
      var.sd <- sqrt(temp$varest[3])
      sd.sd <- sqrt(temp$varest[4])
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df

# Add estimates to the data frame for results

   Results[, 3] <- c(total.est, mean.est, var.est, sd.est)
   Results[, 4] <- c(total.sd, mean.sd, var.sd, sd.sd)

# Calculate confidence bounds

   Results[, 5] <- Results[, 3] - mult*Results[, 4]
   Results[3:4, 5] <- pmax(Results[3:4, 5], 0)
   Results[, 6] <- Results[, 3] + mult*Results[, 4]

# End the section for unstratified data

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

      list(Results=Results, warn.ind=warn.ind, warn.df=warn.df)
   }
}
