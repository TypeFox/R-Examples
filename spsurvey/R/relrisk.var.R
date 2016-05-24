relrisk.var <- function(response, stressor, response.levels, stressor.levels,
   wgt, x, y, stratum.ind, stratum.level, cluster.ind, cluster, wgt1, x1, y1,
   pcfactor.ind, pcfsize, N.cluster, stage1size, support, vartype, warn.ind,
   warn.df, warn.vec) {

################################################################################
# Function: relrisk.var
# Purpose: Calculate values required for estimating variance of the relative 
#          risk estimate
# Programmer: Tom Kincaid
# Date: March 9, 2005
# Last Revised: November 2, 2010
# Description:
#   This function calculates the variance-covariance estimate for the cell and
#   marginal totals used to calculate the relative risk estimate.  Either the
#   simple random sampling (SRS) variance estimator or the local mean variance
#   estimator is calculated, which is subject to user control.  The SRS variance
#   estimator uses the independent random sample approximation to calculate
#   joint inclusion probabilities.  The function can  accomodate single-stage
#   and two-stage samples.
# Arguments:
#   response = the categorical response variable.
#   stressor = the categorical stressor variable.
#   response.levels = category values (levels) for the categorical response 
#     variable, where the first level is used for calculating the relative risk 
#     estimate.  If response.levels equals NULL, then values "Poor" and "Good" 
#     are used for the first level and second level of the response variable, 
#     respectively.  The default is NULL.
#   stressor.levels = category values (levels) for the categorical stressor 
#     variable, where the first level is used for calculating the numerator of 
#     the relative risk estimate and the second level is used for calculating 
#     the denominator of the estimate.  If stressor.levels equals NULL, then 
#     values "Poor" and "Good" are used for the first level and second level of 
#     the stressor variable, respectively.  The default is NULL.
#   wgt = the final adjusted weight (inverse of the sample inclusion 
#     probability) for each site, which is either the weight for a single-stage 
#     sample or the stage two weight for a two-stage sample.
#   x = x-coordinate for location for each site, which is either the x-
#     coordinate for a single-stage sample or the stage two x-coordinate for a 
#     two-stage sample.
#   y = y-coordinate for location for each site, which is either the y-
#     coordinate for a single-stage sample or the stage two y-coordinate for a 
#     two-stage sample.
#   stratum.ind = a logical value that indicates whether the sample is 
#     stratified, where TRUE = a stratified sample and FALSE = not a stratified 
#     sample.
#   stratum.level = the stratum level.
#   cluster.ind = a logical value that indicates whether the sample is a two-
#     stage sample, where TRUE = a two-stage sample and FALSE = not a two-stage 
#     sample.
#   cluster = the stage one sampling unit (primary sampling unit or cluster) 
#     code for each site.
#   wgt1 = the final adjusted stage one weight for each site.
#   x1 = the stage one x-coordinate for location for each site.
#   y1 = the stage one y-coordinate for location for each site.
#   pcfactor.ind = a logical value that indicates whether the population
#     correction factor is used during variance estimation, where TRUE = use the
#     population correction factor and FALSE = do not use the factor.
#   pcfsize = size of the resource, which is required for calculation of finite
#     and continuous population correction factors for a single-stage sample.
#     For a stratified sample this argument must be a vector containing a value
#     for each stratum and must have the names attribute set to identify the
#     stratum codes.
#   N.cluster = the number of stage one sampling units in the resource, which is
#     required for calculation of finite and continuous population correction
#     factors for a two-stage sample.  For a stratified sample this variable
#     must be a vector containing a value for each stratum and must have the
#     names attribute set to identify the stratum codes.
#   stage1size = size of the stage one sampling units of a two-stage sample,
#     which is required for calculation of finite and continuous population
#     correction factors for a two-stage sample and must have the names
#     attribute set to identify the stage one sampling unit codes.  For a
#     stratified sample, the names attribute must be set to identify both
#     stratum codes and stage one sampling unit codes using a convention where
#     the two codes are separated by the & symbol, e.g., "Stratum 1&Cluster 1".
#   support = the support value for each site - the value one (1) for a site
#     from a finite resource or the measure of the sampling unit associated with
#     a site from a continuous resource, which is required for calculation of
#     finite and continuous population correction factors.
#   vartype = the choice of variance estimator, where "Local" = local mean 
#     estimator and "SRS" = SRS estimator.
#   warn.ind = a logical value that indicates whether warning messages were
#     generated, where TRUE = warning messages were generated and FALSE =
#     warning messages were not generated.
#   warn.df = a data frame for storing warning messages.
#   warn.vec = a vector that contains names of the population type, the
#     subpopulation, and an indicator.
# Results:
#   An object in list format composed of a vector named varest, which contains
#   the variance-covariance estimate, a logical variable named warn,ind, which
#   is the indicator for warning messges, and a data frame named warn.df, which
#   contains warning messages.
# Other Functions Required:
#   localmean.weight - calculate the weighting matrix for the local mean 
#     variance estimator
#   localmean.cov - calculate the variance/covariance matrix using the local 
#     mean estimator
################################################################################

# Assign the function name

   fname <- "relrisk.var"

# Branch to handle two-stage and single-stage samples

   if(cluster.ind) {

# Begin the section for a two-stage sample

# Calculate additional required values

      cluster <- factor(cluster)
      cluster.levels <- levels(cluster)
      ncluster <- length(cluster.levels)
      response.lst <- split(response, cluster)
      stressor.lst <- split(stressor, cluster)
      if(vartype == "Local") {
         x2.lst <- split(x, cluster)
         y2.lst <- split(y, cluster)
         x1.u <- as.vector(tapply(x1, cluster, unique))
         y1.u <- as.vector(tapply(y1, cluster, unique))
      }
      wgt2.lst <- split(wgt, cluster)
      wgt1.u <- as.vector(tapply(wgt1, cluster, unique))
      if(pcfactor.ind) {
         support.lst <- split(support, cluster)
      } else {
         support.lst <- NULL
      }
      var.ind <- sapply(split(cluster, cluster), length) > 1

# Calculate estimates of the total of the stage two sampling unit response
# values or residuals and the variance of those totals for each stage one
# sampling unit

      total2est <- matrix(0, ncluster, 4)
      var2est <- array(0, c(4, 4, ncluster))
      for(i in 1:ncluster) {

# Calculate the number of response values

         nresp <- length(response.lst[[i]])

# Create indicator variables for required cells and margins

         Ind1 <- (response.lst[[i]] == response.levels[1])*(stressor.lst[[i]] ==
            stressor.levels[1])
         Ind2 <- (stressor.lst[[i]] == stressor.levels[1])
         Ind3 <- (response.lst[[i]] == response.levels[1])*(stressor.lst[[i]] ==
            stressor.levels[2])
         Ind4 <- (stressor.lst[[i]] == stressor.levels[2])
   
# Calculate the matrix of weighted indicator variables

         rm <- cbind(Ind1, Ind2, Ind3, Ind4) * wgt2.lst[[i]]

# Calculate total estimates for the stage one sampling unit

         total2est[i,] <- apply(rm, 2, sum)

# Adjust the variance-covariance estimator for small sample size

         SRSind <- FALSE
         if(vartype == "Local" && nresp < 4) {
            warn.ind <- TRUE
            act <- "The simple random sampling variance estimator was used.\n"
            if(stratum.ind) {
               warn <- paste("There are less than four response values for stage one sampling unit ", cluster.levels[i], "\nin stratum ", stratum.level, ", the simple random sampling variance estimator \nwas used to calculate variance of the category proportion estimates.\n", sep="")
               warn.df <- rbind(warn.df, data.frame(func=I(fname),
                  subpoptype=warn.vec[1], subpop=warn.vec[2],
                  indicator=warn.vec[3], stratum=I(stratum.level),
                  warning=I(warn), action=I(act)))
            } else {
               warn <- paste("There are less than four response values for stage one sampling unit ", cluster.levels[i], ", \nthe simple random sampling variance estimator was used to calculate variance of \nthe category proportion estimates.\n", sep="")
               warn.df <- rbind(warn.df, data.frame(func=I(fname),
                  subpoptype=warn.vec[1], subpop=warn.vec[2],
                  indicator=warn.vec[3], stratum=NA, warning=I(warn),
                  action=I(act)))
            }
            vartype <- "SRS"
            SRSind <- TRUE
         }

# Calculate the population correction factor for the stage two sample

         pcfactor <- ifelse(pcfactor.ind, (stage1size[i] -
            sum(support.lst[[i]]))/stage1size[i], 1)

# Calculate the variance-covariance estimate for the stage one sampling unit

         if(var.ind[i]) {
            if(vartype == "Local") {
               weight.lst <- localmean.weight(x2.lst[[i]], y2.lst[[i]], 1/wgt2.lst[[i]])
               var2est[i,] <- pcfactor*localmean.cov(rm, weight.lst)
            } else {
               var2est[i,] <- pcfactor*nresp*var(rm)
               if(SRSind)
                  vartype <- "Local"
            }
         }
      }

# Assign the mean variance to stage one sampling units with a single stage two
# sampling unit
      for(j in 1:4) {
         ind <- var2est[,j] == 0
         if(sum(ind) > 0) {
            var.mean <- mean(var2est[!ind,j])
            var2est[ind,j] <- var.mean
         }
      }

# Adjust the variance-covariance estimator for small sample size

      if(vartype == "Local" && ncluster < 4) {
         warn.ind <- TRUE
         act <- "The simple random sampling variance estimator was used.\n"
         if(stratum.ind) {
            warn <- paste("There are less than four stage one sampling units in stratum ", stratum.level, ", \nthe simple random sampling variance estimator was used to calculate variance of \nthe category proportion estimates.\n", sep="")
            warn.df <- rbind(warn.df, data.frame(func=I(fname),
               subpoptype=warn.vec[1], subpop=warn.vec[2],
               indicator=warn.vec[3], stratum=I(stratum.level), warning=I(warn),
               action=I(act)))
         } else {
            warn <- paste("There are less than four stage one sampling units, the simple random sampling \nvariance estimator was used to calculate variance of the category proportion \nestimates.\n", sep="")
            warn.df <- rbind(warn.df, data.frame(func=I(fname),
               subpoptype=warn.vec[1], subpop=warn.vec[2],
               indicator=warn.vec[3], stratum=NA, warning=I(warn),
               action=I(act)))
         }
         vartype <- "SRS"
      }

# Calculate the population correction factor for the stage one sample

      pcfactor <- ifelse(pcfactor.ind, (N.cluster - ncluster)/N.cluster, 1)

# Calculate the variance-covariance estimate

      temp <- 0
      for(i in 1:ncluster)
         temp <- temp + var2est[,,i]*wgt1.u[i]
      if(vartype == "Local") {
         weight.lst <- localmean.weight(x1.u, y1.u, 1/wgt1.u)
         varest <- pcfactor*localmean.cov(total2est*matrix(rep(wgt1.u, 4),
            nrow=ncluster), weight.lst) + temp
      } else {
         varest <- pcfactor*ncluster*var(total2est*matrix(rep(wgt1.u, 4),
            nrow=ncluster)) + temp
      }

# End the section for a two-stage sample

   } else {

# Begin the section for a single-stage sample

# Calculate the number of response values

      nresp <- length(response)

# Create indicator variables for required cells and margins

      Ind1 <- (response == response.levels[1])*(stressor == stressor.levels[1])
      Ind2 <- (stressor == stressor.levels[1])
      Ind3 <- (response == response.levels[1])*(stressor == stressor.levels[2])
      Ind4 <- (stressor == stressor.levels[2])
   
# Calculate the matrix of weighted indicator variables

      rm <- cbind(Ind1, Ind2, Ind3, Ind4) * wgt

# Calculate the population correction factor

      pcfactor <- ifelse(pcfactor.ind, (pcfsize - sum(support))/pcfsize, 1)

# Adjust the variance-covariance estimator for small sample size

      if(vartype == "Local" && nresp < 4) {
         warn.ind <- TRUE
         act <- "The simple random sampling variance estimator was used.\n"
         if(stratum.ind) {
            warn <- paste("There are less than four response values in stratum ", stratum.level, ", \nthe simple random sampling variance estimator was used to calculate variance of \nthe category proportion estimates.\n", sep="")
            warn.df <- rbind(warn.df, data.frame(func=I(fname),
               subpoptype=warn.vec[1], subpop=warn.vec[2],
               indicator=warn.vec[3], stratum=I(stratum.level), warning=I(warn),
               action=I(act)))
         } else {
            warn <- "\nThere are less than four response values, the simple random sampling variance \nestimator was used to calculate variance of the category proportion \nestimates.\n"
            warn.df <- rbind(warn.df, data.frame(func=I(fname),
               subpoptype=warn.vec[1], subpop=warn.vec[2],
               indicator=warn.vec[3], stratum=NA, warning=I(warn),
               action=I(act)))
         }
         vartype <- "SRS"
      }

# Calculate the variance-covariance estimate for the cell and marginal totals

      if (vartype == "Local") {
         wgt.lst <- localmean.weight(x=x, y=y, prb=1/wgt)
         varest <- pcfactor*localmean.cov(rm, wgt.lst)
      } else {
         varest <- pcfactor*nresp*var(rm)
      }

# End the section for a single-stage sample

   }

# Return the variance-covariance estimate, the warning message indicator, and
# the warn.df data frame

   list(varest=varest, warn.ind=warn.ind, warn.df=warn.df)
}
