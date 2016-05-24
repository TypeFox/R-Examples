relrisk.est <- function(response, stressor, response.levels=c("Poor", "Good"),
   stressor.levels=c("Poor", "Good"), wgt, xcoord=NULL, ycoord=NULL,
   stratum=NULL, cluster=NULL, wgt1=NULL, xcoord1=NULL, ycoord1=NULL,
   popcorrect=FALSE, pcfsize=NULL, N.cluster=NULL, stage1size=NULL,
   support=NULL, sizeweight=FALSE, swgt=NULL, swgt1=NULL, vartype="Local",
   conf=95, check.ind=TRUE, warn.ind=NULL, warn.df=NULL, warn.vec=NULL) {

################################################################################
# Function: relrisk.est
# Purpose: Compute the relative risk estimate
# Programmers: Tom Kincaid, Tony Olsen, John Vansickle
# Date: May 4, 2004
# Last Revised: April 6, 2011
# Description:
#   This function calculates the relative risk estimate for a 2x2 table of cell
#   counts defined by a categorical response variable and a categorical
#   explanatory (stressor) variable for an unequal probability design.  Relative
#   risk is the ratio of two probabilities: the numerator is the probability
#   that the first level of the response variable is observed given occurrence
#   of the first level of the stressor variable, and the denominator is the
#   probability that the first level of the response variable is observed given
#   occurrence of the second level of the stressor variable.  The numerator
#   probability and denominator probability are estimated using cell and
#   marginal totals from a 2x2 table of cell counts defined by a categorical
#   response variable and a categorical stressor variable. An estimate of the
#   numerator probability is provided by the ratio of the cell total defined by
#   the first level of response variable and the first level of the stressor
#   variable to the marginal total for the first level of the stressor variable.
#   An estimate of the denominator probability is provided by the ratio of the
#   cell total defined by the first level of response variable and the second
#   level of the stressor variable to the marginal total for the second level of
#   the stressor variable.  Cell and marginal totals are estimated using the
#   Horvitz-Thompson estimator.  The standard error of the log of the relative
#   risk estimate and confidence limits for the estimate also are calculated.
#   The standard error is calculated using a first-order Taylor series
#   linearization (Sarndal et al., 1992).
# Arguments:
#   response = the categorical response variable values.
#   stressor = the categorical explanatory (stressor) variable values.
#   response.levels = category values (levels) for the categorical response 
#     variable, where the first level is used for calculating the numerator and 
#     the denominator of the relative risk estimate.  If response.levels is not 
#     supplied, then values "Poor" and "Good" are used for the first level and 
#     second level of the response variable, respectively.  The default is 
#     c("Poor", "Good").
#   stressor.levels = category values (levels) for the categorical stressor 
#     variable, where the first level is used for calculating the numerator of 
#     the relative risk estimate and the second level is used for calculating 
#     the denominator of the estimate.  If stressor.levels is not supplied, then 
#     values "Poor" and "Good" are used for the first level and second level of 
#     the stressor variable, respectively.  The default is c("Poor", "Good").
#   wgt = the final adjusted weight (inverse of the sample inclusion
#     probability) for each site, which is either the weight for a single-stage
#     sample or the stage two weight for a two-stage sample.
#   xcoord = x-coordinate for location for each site, which is either the x-
#     coordinate for a single-stage sample or the stage two x-coordinate for a
#     two-stage sample.  The default is NULL.
#   ycoord = y-coordinate for location for each site, which is either the y-
#     coordinate for a single-stage sample or the stage two y-coordinate for a
#     two-stage sample.  The default is NULL.
#   stratum = the stratum for each site.  The default is NULL.
#   cluster = the stage one sampling unit (primary sampling unit or cluster) 
#     code for each site.  The default is NULL.
#   wgt1 = the final adjusted stage one weight for each site.  The default is 
#     NULL.
#   xcoord1 = the stage one x-coordinate for location for each site.  The 
#     default is NULL.
#   ycoord1 = the stage one y-coordinate for location for each site.  The 
#     default is NULL.
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
#     two-stage sample.  The default is NULL.
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
# Results:
#   If the function was called by the relrisk.analysis function, then output is
#   an object in list format composed of the Results list, which contains
#   estimates and confidence bounds, the warn.ind logical value, which indicates
#   whether warning messages were generated, and the warn.df data frame, which
#   contains warning messages.  If the function was called directly, then output
#   is the Results list, which contains the following components:
#     RelRisk - the relative risk estimate
#     RRnum - numerator ("elevated" risk) of the relative risk estimate
#     RRdenom - denominator ("baseline" risk) of the relative risk estimate
#     RRlog.se - standard error for the log of the relative risk estimate
#     ConfLimits - confidence limits for the relative risk estimate
#     WeightTotal - sum of the final adjusted weights
#     CellCounts - cell and margin counts for the 2x2 table
#     CellProportions - estimated cell proportions for the 2x2 table
# Other Functions Required:
#   vecprint - takes an input vector and outputs a character string with line 
#     breaks inserted
#   input.check - check input values for errors, consistency, and compatibility 
#     with analytical functions
#   wnas - remove missing values
#   relrisk.var - calculate values required for estimating variance of the
#     relative risk estimate
# Examples:
#   response <- sample(c("Poor", "Good"), 100, replace=TRUE)
#   stressor <- sample(c("Poor", "Good"), 100, replace=TRUE)
#   wgt <- runif(100, 10, 100)
#   relrisk.est(response, stressor, wgt=wgt, vartype="SRS")
#
#   xcoord <- runif(100)
#   ycoord <- runif(100)
#   relrisk.est(response, stressor, wgt=wgt, xcoord=xcoord, ycoord=ycoord)
################################################################################

# As necessary, create a data frame for warning messages
   if(is.null(warn.ind)) {
      warn.ind <- FALSE
      warn.df <- NULL
      warn.vec <- rep(NA, 3)
   }
   fname <- "relrisk.est"

# Check for existence of the response and stressor variables and determine the
# number of values
   if(is.null(response))
      stop("\nValues for the response variable must be provided.")
   if(is.null(stressor))
      stop("\nValues for the stressor variable must be provided.")
   nresp <- length(response)

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

#
# Begin the section that checks for compatibility of input values
#

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
   temp <- input.check(nresp, wgt, NULL, NULL, xcoord, ycoord, stratum.ind,
      stratum, stratum.levels, nstrata, cluster.ind, cluster, cluster.levels,
      ncluster, wgt1, xcoord1, ycoord1, NULL, pcfactor.ind, pcfsize,
      N.cluster, stage1size, support, swgt.ind, swgt, swgt1, vartype, conf)
   pcfsize <- temp$pcfsize
   N.cluster <- temp$N.cluster
   stage1size <- temp$stage1size

# If the sample was stratified and had two stages, then reset cluster to its 
# input value
   if(stratum.ind && cluster.ind)
      cluster <- cluster.in

#
# End the section that checks for compatibility of input values
#

   }

# Remove missing values
   if(vartype == "Local") {
      if(swgt.ind) {
         if(stratum.ind) {
            if(cluster.ind)
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, xcoord=xcoord, ycoord=ycoord, stratum=stratum,
                  cluster=cluster, wgt1=wgt1, xcoord1=xcoord1, ycoord1=ycoord1,
                  swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, xcoord=xcoord, ycoord=ycoord, stratum=stratum,
                  swgt=swgt))
         } else {
            if(cluster.ind)
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, xcoord=xcoord, ycoord=ycoord, cluster=cluster,
                  wgt1=wgt1, xcoord1=xcoord1, ycoord1=ycoord1, swgt=swgt,
                  swgt1=swgt1))
            else
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, xcoord=xcoord, ycoord=ycoord, swgt=swgt))
         }
         response <- temp$response
         stressor <- temp$stressor
         wgt <- temp$wgt
         xcoord <- temp$xcoord
         ycoord <- temp$ycoord
         if(stratum.ind)
            stratum <- temp$stratum
         if(cluster.ind) {
            cluster <- temp$cluster
            wgt1 <- temp$wgt1
            xcoord1 <- temp$xcoord1
            ycoord1 <- temp$ycoord1
            swgt1 <- temp$swgt1
         }
         swgt <- temp$swgt
      } else {
         if(stratum.ind) {
            if(cluster.ind)
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, xcoord=xcoord, ycoord=ycoord, stratum=stratum,
                  cluster=cluster, wgt1=wgt1, xcoord1=xcoord1, ycoord1=ycoord1))
            else
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, xcoord=xcoord, ycoord=ycoord, stratum=stratum))
         } else {
            if(cluster.ind)
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, xcoord=xcoord, ycoord=ycoord, cluster=cluster,
                  wgt1=wgt1, xcoord1=xcoord1, ycoord1=ycoord1))
            else
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, xcoord=xcoord, ycoord=ycoord))
         }
         response <- temp$response
         stressor <- temp$stressor
         wgt <- temp$wgt
         xcoord <- temp$xcoord
         ycoord <- temp$ycoord
         if(stratum.ind)
            stratum <- temp$stratum
         if(cluster.ind) {
            cluster <- temp$cluster
            wgt1 <- temp$wgt1
            xcoord1 <- temp$xcoord1
            ycoord1 <- temp$ycoord1
         }
      }
   } else {
      if(swgt.ind) {
         if(stratum.ind) {
            if(cluster.ind)
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, stratum=stratum, cluster=cluster, wgt1=wgt1,
                  swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, stratum=stratum, swgt=swgt))
         } else {
            if(cluster.ind)
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, cluster=cluster, wgt1=wgt1, swgt=swgt,
                  swgt1=swgt1))
            else
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, swgt=swgt))
         }
         response <- temp$response
         stressor <- temp$stressor
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
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, stratum=stratum, cluster=cluster, wgt1=wgt1))
            else
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, stratum=stratum))
         } else {
            if(cluster.ind)
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt, cluster=cluster, wgt1=wgt1))
            else
               temp <- wnas(list(response=response, stressor=stressor,
                  wgt=wgt))
         }
         response <- temp$response
         stressor <- temp$stressor
         wgt <- temp$wgt
         if(stratum.ind)
            stratum <- temp$stratum
         if(cluster.ind) {
            cluster <- temp$cluster
            wgt1 <- temp$wgt1
         }
      }
   }

# Determine levels of the response and stressor variables
   response <- factor(response)
   temp <- match(levels(response), response.levels)
      if(any(is.na(temp)))
         stop("\nThe category values in response.levels must match the categorical response variable \ncodes.")
   response <- factor(as.vector(response), levels=response.levels)

   stressor <- factor(stressor)
   temp <- match(levels(stressor), stressor.levels)
      if(any(is.na(temp)))
         stop("\nThe category values in stressor.levels must match the categorical stressor variable \ncodes.")
   stressor <- factor(as.vector(stressor), levels=stressor.levels)

# For a stratified sample, check for strata that no longer contain any values,
# as necesssary remove strata that contain a single value and output a warning
# message
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
            response <- response[!stratum.i]
            stressor <- stressor[!stratum.i]
            wgt <- wgt[!stratum.i]
            if(vartype == "Local") {
               xcoord <- xcoord[!stratum.i]
               ycoord <- ycoord[!stratum.i]
            }
            stratum <- stratum[!stratum.i]
            if(cluster.ind) {
               cluster <- cluster[!stratum.i]
               wgt1 <- wgt1[!stratum.i]
               if(vartype == "Local") {
                  xcoord1 <- xcoord1[!stratum.i]
                  ycoord1 <- ycoord1[!stratum.i]
               }
            }
            if(swgt.ind) {
               swgt <- swgt[!stratum.i]
               if(cluster.ind)
                  swgt1 <- swgt1[!stratum.i]
            }
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
   nresp <- length(response)
   if(nresp == 0)
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

# Compute the sum of the weights
   if(swgt.ind) {
      if(cluster.ind) {
         popsize.hat <- sum(wgt*swgt*wgt1*swgt1)
      } else {
         popsize.hat <- sum(wgt*swgt)
      }
   } else {
      if(cluster.ind) {
         popsize.hat <- sum(wgt*wgt1)
      } else {
         popsize.hat <- sum(wgt)
      }
   }

#
# Branch to handle stratified and unstratified data
#

   if(stratum.ind) {

#
# Begin the section for stratified data
#

# Initialize variables for all strata combined
      wgt.total <- 0
      varest <- 0

#
# Begin the subsection for individual strata
#

      for(i in 1:nstrata) {

# Calculate required values
         stratum.i <- stratum == stratum.levels[i]
         response.st <- response[stratum.i]
         stressor.st <- stressor[stratum.i]
         if(swgt.ind) {
            if(cluster.ind) {
               wgt.st <- wgt[stratum.i]*swgt[stratum.i]
               wgt1.st <- wgt1[stratum.i]*swgt1[stratum.i]
            } else {
               wgt.st <- wgt[stratum.i]*swgt[stratum.i]
            }
         } else {
            if(cluster.ind) {
               wgt.st <- wgt[stratum.i]
               wgt1.st <- wgt1[stratum.i]
            } else {
               wgt.st <- wgt[stratum.i]
            }
         }

# Compute the 2x2 table of weight totals
         if(cluster.ind) {
            wgt.total.st <- tapply(wgt.st*wgt1.st, list(response=response.st,
               stressor=stressor.st), sum)
         } else {
            wgt.total.st <- tapply(wgt.st, list(response=response.st,
               stressor=stressor.st), sum)
         }
         wgt.total.st[is.na(wgt.total.st)] <- 0

# Calculate the variance-covariance estimate for the cell and marginal totals
         if(cluster.ind) {
            temp <- relrisk.var(response.st, stressor.st, response.levels,
               stressor.levels, wgt.st, xcoord[stratum.i], ycoord[stratum.i],
               stratum.ind, stratum.levels[i], cluster.ind,
               cluster[stratum.i], wgt1.st, xcoord1[stratum.i],
               ycoord1[stratum.i], pcfactor.ind, NULL, N.cluster[i],
               stage1size[[i]], support[stratum.i], vartype, warn.ind,
               warn.df, warn.vec)
         } else {
            temp <- relrisk.var(response.st, stressor.st, response.levels,
               stressor.levels, wgt.st, xcoord[stratum.i], ycoord[stratum.i],
               stratum.ind, stratum.levels[i], cluster.ind,
               pcfactor.ind=pcfactor.ind, pcfsize=pcfsize[i],
               support=support[stratum.i], vartype=vartype,
               warn.ind=warn.ind, warn.df=warn.df,
               warn.vec=warn.vec)
         }
         varest.st <- temp$varest
         warn.ind <- temp$warn.ind
         warn.df <- temp$warn.df

# Add estimates to the variables for all strata combined
         wgt.total <- wgt.total + wgt.total.st
         varest <- varest + varest.st

#
# End the subsection for individual strata
#

      }

# Assign required cell and marginal weight totals
      total1 <- wgt.total[response.levels[1], stressor.levels[1]]
      total2 <- sum(wgt.total[,stressor.levels[1]])
      total3 <- wgt.total[response.levels[1], stressor.levels[2]]
      total4 <- sum(wgt.total[,stressor.levels[2]])
   
# Calculate the estimate of relative risk for all strata combined
      if(total2 == 0 || total4 == 0) {
         rr <- NA
         rr.num <- NA
         rr.denom <- NA
         warn.ind <- TRUE
         temp <- ifelse(total2 == 0, stressor.levels[1], stressor.levels[2])
         warn <- paste("Since there are no observations for level \"", temp, "\" of the stressor \nvariable, the relative risk estimate and its standard error cannot be \ncalculated for stratum \"", stratum.levels[i], "\".  Also, the stratum \nwas removed from the analysis.\n", sep="")
         act <- paste("The relative risk estimate and its standard error were not calculated for \nstratum \"", stratum.levels[i], "\".  Also, the stratum was removed from \nthe analysis.\n", sep="")
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2],
            indicator=warn.vec[3], stratum=NA, warning=I(warn),
            action=I(act)))
      } else if(total1 == 0 && total3 != 0) {
         rr <- 0
         rr.num <- 0
         rr.denom <- total3/total4
         warn.ind <- TRUE
         warn <- paste("Since there are no observations for the cell defined by level \"", response.levels[1], "\" \nof the response variable and level \"", stressor.levels[1], "\" of the stressor \nvariable, the relative risk estimate is zero and standard error of the relative \nrisk estimate cannot be calculated for stratum \"", stratum.levels[i], "\".  \nAlso, the stratum was removed from the analysis.\n", sep="")
         act <- paste("Standard error of the relative risk estimate was not calculated for stratum \n\"", stratum.levels[i], "\".  Also, the stratum was removed from the \nanalysis.\n", sep="")
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2],
            indicator=warn.vec[3], stratum=NA, warning=I(warn),
            action=I(act)))
      } else if(total1 == 0 && total3 == 0) {
         rr <- NA
         rr.num <- total1/total2
         rr.denom <- total3/total4
         warn.ind <- TRUE
         warn <- paste("Since there are no observations for the cell defined by level \"", response.levels[1], "\" \nof the response variable and level \"", stressor.levels[1], "\" of the stressor \nvariable and for the cell defined by level \"", response.levels[1], "\" of the \nresponse variable and level \"", stressor.levels[2], "\" of the stressor variable, \nthe relative risk estimate and its standard error cannot be calculated for \nstratum \"", stratum.levels[i], "\".  Also, the stratum was removed from \nthe analysis.\n", sep="")
         act <- paste("The relative risk estimate and its standard error were not calculated for \nstratum \"", stratum.levels[i], "\".  Also, the stratum was removed from \nthe analysis.\n", sep="")
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2],
            indicator=warn.vec[3], stratum=NA, warning=I(warn),
            action=I(act)))
      } else if(total3 == 0) {
         rr <- NA
         rr.num <- total1/total2
         rr.denom <- total3/total4
         warn.ind <- TRUE
         warn <- paste("Since there are no observations for the cell defined by level \"", response.levels[1], "\" \nof the response variable and level \"", stressor.levels[2], "\" of the stressor \nvariable, the relative risk estimate and its standard error cannot be \ncalculated for stratum \"", stratum.levels[i], "\".  Also, the stratum \nwas removed from the analysis.\n", sep="")
         act <- paste("The relative risk estimate and its standard error were not calculated for \nstratum \"", stratum.levels[i], "\".  Also, the stratum was removed from \nthe analysis.\n", sep="")
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2],
            indicator=warn.vec[3], stratum=NA, warning=I(warn),
            action=I(act)))
      } else {
         rr <- (total1*total4) / (total2*total3)
         rr.num <- total1/total2
         rr.denom <- total3/total4
      }

# Calculate the standard error estimate of the log of relative risk for all
# strata combined
      if(any(c(total1, total2, total3, total4) == 0)) {
         rrlog.se <- NA
      } else {
         pder <- 1/c(total1, -total2, -total3, total4)
         rrlog.se <- sqrt(t(pder) %*% varest %*% pder)
      }

#
# End the section for stratified data
#

   } else {

#
# Begin the section for unstratified data
#

# Check whether the vector of response values contains a single element
      if(nresp == 1)
         stop("\nEstimates cannot be calculated since the vector of response values contains a \nsingle element.")

# If the sample is size-weighted, calculate combined weights
      if(swgt.ind) {
         if(cluster.ind) {
            wgt <- wgt*swgt
            wgt1 <- wgt1*swgt1
         } else {
            wgt <- wgt*swgt
         }
      }

# Compute the 2x2 table of weight totals
      if(cluster.ind) {
         wgt.total <- tapply(wgt*wgt1, list(response=response,
            stressor=stressor), sum)
      } else {
         wgt.total <- tapply(wgt, list(response=response, stressor=stressor),
            sum)
      }
      wgt.total[is.na(wgt.total)] <- 0

# Calculate required cell and marginal weight totals
      total1 <- wgt.total[response.levels[1], stressor.levels[1]]
      total2 <- sum(wgt.total[,stressor.levels[1]])
      total3 <- wgt.total[response.levels[1], stressor.levels[2]]
      total4 <- sum(wgt.total[,stressor.levels[2]])
   
# Calculate the estimate of relative risk
      if(total2 == 0 || total4 == 0) {
         rr <- NA
         rr.num <- NA
         rr.denom <- NA
         warn.ind <- TRUE
         temp <- ifelse(total2 == 0, stressor.levels[1], stressor.levels[2])
         warn <- paste("Since there are no observations for level \"", temp, "\" of the stressor \nvariable, the relative risk estimate and its standard error cannot be \ncalculated.\n", sep="")
         act <- "The relative risk estimate and its standard error were not calculated.\n"
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2], indicator=warn.vec[3],
            stratum=NA, warning=I(warn), action=I(act)))
      } else if(total1 == 0 && total3 != 0) {
         rr <- 0
         rr.num <- 0
         rr.denom <- total3/total4
         warn.ind <- TRUE
         warn <- paste("Since there are no observations for the cell defined by level \"", response.levels[1], "\" \nof the response variable and level \"", stressor.levels[1], "\" of the stressor variable, \nthe relative risk estimate is zero and standard error of the relative risk \nestimate cannot be calculated.\n", sep="")
         act <- "Standard error of the relative risk estimate was not calculated.\n"
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2], indicator=warn.vec[3],
            stratum=NA, warning=I(warn), action=I(act)))
      } else if(total1 == 0 && total3 == 0) {
         rr <- NA
         rr.num <- total1/total2
         rr.denom <- total3/total4
         warn.ind <- TRUE
         warn <- paste("Since there are no observations for the cell defined by level \"", response.levels[1], "\" \nof the response variable and level \"", stressor.levels[1], "\" of the stressor \nvariable and for the cell defined by level \"", response.levels[1], "\" of the \nresponse variable and level \"", stressor.levels[2], "\" of the stressor variable, \nthe relative risk estimate and its standard error cannot be calculated.\n", sep="")
         act <- "The relative risk estimate and its standard error were not calculated.\n"
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2], indicator=warn.vec[3],
            stratum=NA, warning=I(warn), action=I(act)))
      } else if(total3 == 0) {
         rr <- NA
         rr.num <- total1/total2
         rr.denom <- total3/total4
         warn.ind <- TRUE
         warn <- paste("Since there are no observations for the cell defined by level \"", response.levels[1], "\" \nof the response variable and level \"", stressor.levels[2], "\" of the stressor \nvariable, the relative risk estimate and its standard error cannot be \ncalculated.\n", sep="")
         act <- "The relative risk estimate and its standard error were not calculated.\n"
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2], indicator=warn.vec[3],
            stratum=NA, warning=I(warn), action=I(act)))
      } else {
         rr <- (total1*total4) / (total2*total3)
         rr.num <- total1/total2
         rr.denom <- total3/total4
      }

# Determine whether the standard error can be calculated
      if(any(c(total1, total2, total3, total4) == 0)) {
         rrlog.se <- NA
      } else {

# Calculate the variance-covariance estimate for the cell and marginal totals

         if(cluster.ind) {
            temp <- relrisk.var(response, stressor, response.levels,
               stressor.levels, wgt, xcoord, ycoord, stratum.ind, NULL,
               cluster.ind, cluster, wgt1, xcoord1, ycoord1, pcfactor.ind, NULL,
               N.cluster, stage1size, support, vartype, warn.ind, warn.df,
               warn.vec)
         } else {
            temp <- relrisk.var(response, stressor, response.levels,
               stressor.levels, wgt, xcoord, ycoord, stratum.ind, NULL,
               cluster.ind, pcfactor.ind=pcfactor.ind, pcfsize=pcfsize,
               support=support, vartype=vartype, warn.ind=warn.ind,
               warn.df=warn.df, warn.vec=warn.vec)
         }
         varest <- temp$varest
         warn.ind <- temp$warn.ind
         warn.df <- temp$warn.df

# Calculate the standard error estimate of the log of relative risk
         pder <- 1/c(total1, -total2, -total3, total4)
         rrlog.se <- sqrt(t(pder) %*% varest %*% pder)
      }

#
# End section for unstratified data
#

   }

# Calculate confidence limits for the estimate of relative risk
   if(is.na(rrlog.se)) {
      cl <- NA
   } else {
      cl <- c(exp(log(rr) - rrlog.se * mult), exp(log(rr) + rrlog.se * mult))
   }

# Calculate the table of cell and margin counts
   cc <- ftable(addmargins(table(list(response=response,
      stressor=stressor))))

# Calculate the table of cell and margin proportion estimates
   cp <- ftable(addmargins(wgt.total/popsize.hat))

# Create the Results list
   Results <- list(RelRisk=rr, RRnum=rr.num, RRdenom=rr.denom,
      RRlog.se=rrlog.se, ConfLimits=cl, WeightTotal=sum(wgt.total),
      CellCounts=cc, CellProportions=cp )

#
# Depending on whether the function was called directly or was called by
# relrisk.analysis, return appropriate results
#

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

# Return the Results list
      Results

   } else {

# Return the Results list, the warn.ind logical value, and the warn.df
# data frame
      list(Results=Results, warn.ind=warn.ind, warn.df=warn.df)
   }
}
