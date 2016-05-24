category.est <- function(catvar, wgt, x=NULL, y=NULL, stratum=NULL,
   cluster=NULL, wgt1=NULL, x1=NULL, y1=NULL, popsize=NULL, popcorrect=FALSE,
   pcfsize=NULL, N.cluster=NULL, stage1size=NULL, support=NULL,
   sizeweight=FALSE, swgt=NULL, swgt1=NULL, vartype="Local", conf=95,
   check.ind=TRUE, warn.ind=NULL, warn.df=NULL, warn.vec=NULL) {

################################################################################
# Function: category.est
# Programmer: Tom Kincaid
# Date: August 4, 2000
# Last Revised: April 6, 2011
# Description:
#   This function estimates proportion (expressed as percent) and size of a
#   resource in each of a set of categories and can also be used to estimate
#   proportion and size for site status categories.  Upper and lower confidence
#   bounds also are estimated.  Proportion estimates are calculated using the
#   Horvitz-Thompson ratio estimator, i.e., the ratio of two Horvitz-Thompson
#   estimators.  The numerator of the ratio estimates the size of the category.
#   The denominator of the ratio estimates the size of the resource.  
#   Variance estimates for the proportion estimates are calculated using either 
#   the local mean variance estimator or the simple random sampling (SRS) 
#   variance estimator.  The choice of variance estimator is subject to user
#   control.  The local mean variance estimator requires the x-coordinate and
#   the y-coordinate of each site.  The SRS variance estimator uses the
#   independent random sample approximation to calculate joint inclusion
#   probabilities.  Confidence bounds are calculated using a Normal distribution
#   multiplier.  For a finite resource size is the number of units in the
#   resource.  For an extensive resource size is the measure (extent) of the
#   resource, i.e., length, area, or volume.  Size estimates are calculated
#   using the Horvitz-Thompson estimator.  Variance estimates for the size
#   estimates are calculated using either the local mean variance estimator or
#   the SRS variance estimator.  The function can accommodate a stratified
#   sample.  For a stratified sample, separate estimates and standard errors are
#   calculated for each stratum, which are used to produce estimates and
#   standard errors for all strata combined.  Strata that contain a single value
#   are removed.  For a stratified sample, when either the size of the resource
#   or the sum of the size-weights for the resource is provided for each
#   stratum, those values are used as stratum weights for calculating the
#   estimates and standard errors for all strata combined.  In addition, when
#   either of those known values is provided for each stratum, size estimates
#   are obtained by multiplying the proportion estimate, i.e., the Horvitz-
#   Thompson ratio estimator, by the known value for the stratum.  For a
#   stratified sample when neither the size of the resource nor the sum of the
#   size-weights of the resource is provided for each stratum, estimated values
#   are used as stratum weights for calculating the estimates and standard
#   errors for all strata combined.  The function can accommodate single-stage
#   and two-stage samples for both stratified and unstratified sampling designs.
#   Finite population and continuous population correction factors can be
#   utilized in variance estimation.  The function checks for compatibility of
#   input values and removes missing values.
# Arguments:
#   catvar = the value of the categorical response variable or the site status
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
#   If the function was called by the cat.analysis function, then output is an
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
#   catvar.prop - calculate variance of the proportion estimates
#   catvar.size - calculate variance of the size estimates
# Examples:
#   catvar <- rep(c("north", "south", "east", "west"), rep(25, 4))
#   wgt <- runif(100, 10, 100)
#   category.est(catvar, wgt, vartype="SRS")
#
#   x <- runif(100)
#   y <- runif(100)
#   category.est(catvar, wgt, x, y)
################################################################################

# As necessary, create a data frame for warning messages

   if(is.null(warn.ind)) {
      warn.ind <- FALSE
      warn.df <- NULL
      warn.vec <- rep(NA, 3)
   }
   fname <- "category.est"

# Check for existence of the category variable and determine the number of values

   if(is.null(catvar))
      stop("\nValues for the categorical variable must be provided.")
   nresp <- length(catvar)

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
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y, stratum=stratum, cluster=cluster, wgt1=wgt1, x1=x1, y1=y1, swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y, stratum=stratum, swgt=swgt))
         } else {
            if(cluster.ind)
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y, cluster=cluster, wgt1=wgt1, x1=x1, y1=y1, swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y, swgt=swgt))
         }
         catvar <- temp$catvar
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
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y, stratum=stratum, cluster=cluster, wgt1=wgt1, x1=x1, y1=y1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y, stratum=stratum))
         } else {
            if(cluster.ind)
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y, cluster=cluster, wgt1=wgt1, x1=x1, y1=y1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y))
         }
         catvar <- temp$catvar
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
               temp <- wnas(list(catvar=catvar, wgt=wgt, stratum=stratum, cluster=cluster, wgt1=wgt1, swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt, stratum=stratum, swgt=swgt))
         } else {
            if(cluster.ind)
               temp <- wnas(list(catvar=catvar, wgt=wgt, cluster=cluster, wgt1=wgt1, swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt, swgt=swgt))
         }
         catvar <- temp$catvar
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
               temp <- wnas(list(catvar=catvar, wgt=wgt, stratum=stratum, cluster=cluster, wgt1=wgt1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt, stratum=stratum))
         } else {
            if(cluster.ind)
               temp <- wnas(list(catvar=catvar, wgt=wgt, cluster=cluster, wgt1=wgt1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt))
         }
         catvar <- temp$catvar
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
            catvar <- catvar[!stratum.i]
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

# Check whether the vector of categorical variable values is empty

   if(length(catvar) == 0)
      stop("\nEstimates cannot be calculated since the vector of categorical variable values \nis empty.")

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

# Determine levels of the category variable

   catvar.levels <- levels(factor(catvar))
   nlevels <- length(catvar.levels)

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

# Begin the subsection for proportion estimates

# Create the data frame for estimates for all strata combined

   Results <- data.frame(array(0, c(nlevels+1, 10)))
   dimnames(Results) <- list(1:(nlevels+1), c("Category", "NResp", "Estimate.P",
      "StdError.P", paste("LCB", conf, "Pct.P", sep=""), paste("UCB", conf,
      "Pct.P", sep=""), "Estimate.U", "StdError.U", paste("LCB", conf, "Pct.U",
      sep=""), paste("UCB", conf, "Pct.U", sep="")))
   Results[,1] <- c(catvar.levels, "Total")

# For known popsize, create matrices for proportion and variance estimates

   if(!is.null(popsize)) {
      prop.popsize <- matrix(0, nlevels, nstrata)
      varest.popsize <- matrix(0, nlevels, nstrata)
   }

# Begin the subsection for individual strata

   for(i in 1:nstrata) {

# Calculate the proportion estimates

   stratum.i <- stratum == stratum.levels[i]
   z <- factor(catvar[stratum.i])
   z.levels <- levels(z)
   m <- length(z.levels)
   if(swgt.ind) {
      if(cluster.ind) {
         w2 <- wgt[stratum.i]*swgt[stratum.i]
         w1 <- wgt1[stratum.i]*swgt1[stratum.i]
         prop <- tapply(w2*w1, z, sum) / popsize.hat[i]
         nval <- tapply(w2, z, length)
      } else {
         w <- wgt[stratum.i]*swgt[stratum.i]
         prop <- tapply(w, z, sum) / popsize.hat[i]
         nval <- tapply(w, z, length)
      }
   } else {
      if(cluster.ind) {
         w2 <- wgt[stratum.i]
         w1 <- wgt1[stratum.i]
         prop <- tapply(w2*w1, z, sum) / popsize.hat[i]
         nval <- tapply(w2, z, length)
      } else {
         w <- wgt[stratum.i]
         prop <- tapply(w, z, sum) / popsize.hat[i]
         nval <- tapply(w, z, length)
      }
   }

# Calculate the variance estimates

   if(cluster.ind) {
      temp <- catvar.prop(z, w2, x[stratum.i], y[stratum.i], prop, stratum.ind,
         stratum.levels[i], cluster.ind, cluster[stratum.i], w1, x1[stratum.i],
         y1[stratum.i], pcfactor.ind, NULL, N.cluster[i], stage1size[[i]],
         support[stratum.i], vartype, warn.ind, warn.df, warn.vec)
   } else {
      temp <- catvar.prop(z, w, x[stratum.i], y[stratum.i], prop, stratum.ind,
         stratum.levels[i], cluster.ind, pcfactor.ind=pcfactor.ind,
         pcfsize=pcfsize[i], support=support[stratum.i], vartype=vartype,
         warn.ind=warn.ind, warn.df=warn.df, warn.vec=warn.vec)
   }
   varest <- temp$varest
   warn.ind <- temp$warn.ind
   warn.df <- temp$warn.df

# Combine estimates in a matrix

   rslt <- cbind(nval, prop, varest)
   if(m < nlevels) {
      temp <- array(0, c(nlevels, 3))
      k <- 1
      for(j in 1:nlevels) {
         if(k <= m) {
            if(z.levels[k] != catvar.levels[j]) {
               temp[j,] <- rep(NA, 3)
            } else {
               temp[j,] <- rslt[k,]
               k <- k+1
            }
         } else {
            temp[j,] <- rep(NA, 3)
         }
      }
      rslt <- temp
   }

# For known popsize, add estimates to the  matrices for proportion and variance
# estimates

   if(!is.null(popsize)) {
      prop.popsize[,i] <- rslt[,2]
      varest.popsize[,i] <- rslt[,3]
   }

# Add estimates to the data frame for all strata combined

   Results[1:nlevels, 2][!is.na(rslt[,1])] <- Results[1:nlevels, 2][!is.na(rslt[,1])] + rslt[,1][!is.na(rslt[,1])]
   if(!is.null(popsize)) {
      Results[1:nlevels, 3][!is.na(rslt[,2])] <- Results[1:nlevels, 3][!is.na(rslt[,2])] + (popsize[i]/sum.popsize)*rslt[,2][!is.na(rslt[,2])]
      Results[1:nlevels, 4][!is.na(rslt[,3])] <- Results[1:nlevels, 4][!is.na(rslt[,3])] + ((popsize[i]/sum.popsize)^2)*rslt[,3][!is.na(rslt[,3])]
   } else {
      Results[1:nlevels, 3][!is.na(rslt[,2])] <- Results[1:nlevels, 3][!is.na(rslt[,2])] + (popsize.hat[i]/sum.popsize.hat)*rslt[,2][!is.na(rslt[,2])]
      Results[1:nlevels, 4][!is.na(rslt[,3])] <- Results[1:nlevels, 4][!is.na(rslt[,3])] + ((popsize.hat[i]/sum.popsize.hat)^2)*rslt[,3][!is.na(rslt[,3])]
   }

# End the subsection for individual strata

   }

# Begin the subsection for all strata combined

# Calculate standard errors and confidence bounds

   Results[1:nlevels, 4] <- sqrt(Results[1:nlevels, 4])
   Results[1:nlevels, 5] <- 100*pmax(Results[1:nlevels, 3] - mult*Results[1:nlevels, 4], 0)
   Results[1:nlevels, 6] <- 100*pmin(Results[1:nlevels, 3] + mult*Results[1:nlevels, 4], 1)
   Results[1:nlevels, 3] <- 100*Results[1:nlevels, 3]
   Results[1:nlevels, 4] <- 100*Results[1:nlevels, 4]

# Provide appropriate values for the "Total" row

   if(!is.null(popsize))
      Results[nlevels+1, 2:6] <- c(sum(Results[1:nlevels, 2]), 100, NA, NA, NA)
   else
      Results[nlevels+1, 2:6] <- c(sum(Results[1:nlevels, 2]), 100, 0, 100, 100)

# End the subsection for all strata combined

# End the subsection for proportion estimates

# Begin the subsection for size estimates

# Begin the subsection for individual strata

   for(i in 1:nstrata) {

# Estimate the size of each category

   stratum.i <- stratum == stratum.levels[i]
   z <- factor(catvar[stratum.i])
   z.levels <- levels(z)
   m <- length(z.levels)
   if(swgt.ind) {
      if(!is.null(popsize)) {
         size <- c(popsize[i]*prop.popsize[,i], popsize[i])
      } else {
         if(cluster.ind) {
            w2 <- wgt[stratum.i]*swgt[stratum.i]
            w1 <- wgt1[stratum.i]*swgt1[stratum.i]
            size <- c(tapply(w2*w1, z, sum), popsize.hat[i])
         } else {
            w <- wgt[stratum.i]*swgt[stratum.i]
            size <- c(tapply(w, z, sum), popsize.hat[i])
         }
      }
   } else {
      if(!is.null(popsize)) {
         size <- c(popsize[i]*prop.popsize[,i], popsize[i])
      } else {
         if(cluster.ind) {
            w2 <- wgt[stratum.i]
            w1 <- wgt1[stratum.i]
            size <- c(tapply(w2*w1, z, sum), popsize.hat[i])
         } else {
            w <- wgt[stratum.i]
            size <- c(tapply(w, z, sum), popsize.hat[i])
         }
      }
   }

# Calculate the variance estimates

   if(!is.null(popsize)) {
      varest <- c((popsize[i]^2)*varest.popsize[,i], NA)
   } else {
      if(cluster.ind) {
         temp <- catvar.size(z, w2, x[stratum.i], y[stratum.i], size,
            stratum.ind, stratum.levels[i], cluster.ind, cluster[stratum.i], w1,
            x1[stratum.i], y1[stratum.i], pcfactor.ind, NULL, N.cluster[i],
            stage1size[[i]], support[stratum.i], vartype, warn.ind, warn.df, warn.vec)
      } else {
         temp <- catvar.size(z, w, x[stratum.i], y[stratum.i], size,
            stratum.ind, stratum.levels[i], cluster.ind,
            pcfactor.ind=pcfactor.ind, pcfsize=pcfsize[i],
            support=support[stratum.i], vartype=vartype, warn.ind=warn.ind,
            warn.df=warn.df, warn.vec=warn.vec)
      }
      varest <- temp$varest
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df
   }

# Combine estimates in a matrix

   rslt <- cbind(size, varest)
   if(is.null(popsize) && m < nlevels) {
      temp <- array(0, c(nlevels+1, 2))
      k <- 1
      for(j in 1:nlevels) {
         if(k <= m) {
            if(z.levels[k] != catvar.levels[j]) {
               temp[j,] <- rep(NA, 2)
            } else {
                temp[j,] <- rslt[k,]
                k <- k+1
            }
         } else {
            temp[j,] <- rep(NA, 2)
         }
      }
      temp[nlevels+1,] <- rslt[m+1,]
      rslt <- temp
   }

# Add estimates to the data frame for all strata combined

   Results[,7][!is.na(rslt[,1])] <- Results[,7][!is.na(rslt[,1])] +
      rslt[,1][!is.na(rslt[,1])]
   Results[,8][!is.na(rslt[,2])] <- Results[,8][!is.na(rslt[,2])] +
      rslt[,2][!is.na(rslt[,2])]

# End the subsection for individual strata

   }

# Begin the subsection for all strata combined

# Calculate standard errors and confidence bounds

   Results[,8] <- sqrt(Results[,8])
   Results[,9] <- pmax(Results[,7] - mult*Results[,8], 0)
   if(!is.null(popsize)) {
      Results[,10] <- pmin(Results[,7] + mult*Results[,8], sum.popsize)
   } else {
      Results[,10] <- Results[,7] + mult*Results[,8]
   }

# For known popsize, adjust values for the "Total" row

   if(!is.null(popsize))
      Results[nlevels+1, 8:10] <- rep(NA, 3)

# End the subsection for all strata combined

# End the subsection for size estimates

# End the section for stratified data

   } else {

# Begin the section for unstratified data

# Check whether the vector of category variable values contains a single element

   if(length(catvar) == 1)
      stop("\nEstimates cannot be calculated since the vector of categorical variable values \ncontains a single element.")

# Create the data frame for estimates

   Results <- data.frame(array(0, c(nlevels+1, 10)))
   dimnames(Results) <- list(1:(nlevels+1), c("Category", "NResp", "Estimate.P",
      "StdError.P", paste("LCB", conf, "Pct.P", sep=""), paste("UCB", conf,
      "Pct.P", sep=""), "Estimate.U", "StdError.U", paste("LCB", conf, "Pct.U",
      sep=""), paste("UCB", conf, "Pct.U", sep="")))
   Results[,1] <- c(catvar.levels, "Total")

# Begin the subsection for proportion estimates

# Calculate the number of values and the proportion of each category

   z <- factor(catvar)
   if(swgt.ind) {
      if(cluster.ind) {
         w2 <- wgt*swgt
         w1 <- wgt1*swgt1
         prop <- tapply(w2*w1, z, sum) / popsize.hat
         nval <- tapply(w2, z, length)
      } else {
         w <- wgt*swgt
         prop <- tapply(w, z, sum) / popsize.hat
         nval <- tapply(w, z, length)
      }
   } else {
      if(cluster.ind) {
         w2 <- wgt
         w1 <- wgt1
         prop <- tapply(w2*w1, z, sum) / popsize.hat
         nval <- tapply(w2, z, length)
      } else {
         w <- wgt
         prop <- tapply(w, z, sum) / popsize.hat
         nval <- tapply(w, z, length)
      }
   }

# Calculate the standard error estimates

   if(cluster.ind) {
      temp <- catvar.prop(z, w2, x, y, prop, stratum.ind, NULL, cluster.ind,
         cluster, w1, x1, y1, pcfactor.ind, NULL, N.cluster, stage1size,
         support, vartype, warn.ind, warn.df, warn.vec)
   } else {
      temp <- catvar.prop(z, w, x, y, prop, stratum.ind, NULL, cluster.ind,
         pcfactor.ind=pcfactor.ind, pcfsize=pcfsize, support=support,
         vartype=vartype, warn.ind=warn.ind, warn.df=warn.df, warn.vec=warn.vec)
   }
   sdest <- sqrt(temp$varest)
   warn.ind <- temp$warn.ind
   warn.df <- temp$warn.df

# Calculate confidence bounds

   lbound <- 100*pmax(prop - mult*sdest, 0)
   ubound <- 100*pmin(prop + mult*sdest, 1)
   prop <- 100*prop
   sdest <- 100*sdest

# Add estimates to the data frame for results

   Results[1:nlevels,2:6] <- cbind(nval, prop, sdest, lbound, ubound)

# Provide appropriate values for the "Total" row

   if(!is.null(popsize))
      Results[nlevels+1, 2:6] <- c(sum(Results[1:nlevels, 2]), 100, NA, NA, NA)
   else
      Results[nlevels+1, 2:6] <- c(sum(Results[1:nlevels, 2]), 100, 0, 100, 100)

# End the subsection for proportion estimates

# Begin the subsection for size estimates

# Estimate the size of each category

   if(!is.null(popsize))
      size <- c(popsize*Results[1:nlevels, 3]/100, popsize)
   else
      size <- c(popsize.hat*(prop/100), popsize.hat)

# Calculate the standard error estimates

   if(!is.null(popsize)) {
      sdest <- popsize*Results[,4]/100
   } else {
      if(cluster.ind) {
         temp <- catvar.size(z, w2, x, y, size, stratum.ind, NULL, cluster.ind,
            cluster, w1, x1, y1, pcfactor.ind, NULL, N.cluster, stage1size,
            support, vartype, warn.ind, warn.df, warn.vec)
      } else {
         temp <- catvar.size(z, w, x, y, size, stratum.ind, NULL, cluster.ind,
            pcfactor.ind=pcfactor.ind, pcfsize=pcfsize, support=support,
            vartype=vartype, warn.ind=warn.ind, warn.df=warn.df,
            warn.vec=warn.vec)
      }
      sdest <- sqrt(temp$varest)
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df
   }

# Calculate confidence bounds

   lbound <- pmax(size - mult*sdest, 0)
   if(!is.null(popsize)) {
      ubound <- pmin(size + mult*sdest, popsize)
   } else {
      ubound <- size + mult*sdest
   }

# Add estimates to the data frame for results

   Results[,7:10] <- cbind(size, sdest, lbound, ubound)

# End the subsection for size estimates

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
