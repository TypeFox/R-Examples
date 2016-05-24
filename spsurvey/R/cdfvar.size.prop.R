cdfvar.size.prop <- function(z, wgt, x, y, val, cdfest, stratum.ind,
   stratum.level, cluster.ind, cluster, wgt1, x1, y1, pcfactor.ind, pcfsize,
   N.cluster, stage1size, support, swgt, swgt1, vartype, warn.ind, warn.df,
   warn.vec) {

################################################################################
# Function: cdfvar.size.prop
# Programmer: Tom Kincaid
# Date: July 12, 2000
# Last Revised: August 18, 2010
# Description:
#   This function calculates variance estimates of the estimated size-weighted
#   cumulative distribution function (CDF) for the proportion of a finite
#   resource.  The set of values at which the CDF is estimated is supplied to
#   the function.  Either the simple random sampling (SRS) variance estimator or
#   the local mean variance estimator is calculated, which is subject to user
#   control.  The SRS variance estimator uses the independent random sample
#   approximation to calculate joint inclusion probabilities.  The function can 
#   accomodate single-stage and two-stage samples.  Finite population and 
#   continuous population correction factors can be utilized in variance 
#   estimation.
# Arguments:
#   z = the response value for each site.
#   wgt = the final adjusted weight (inverse of the sample inclusion
#     probability) for each site, which is either the weight for a single-
#     stage sample or the stage two weight for a two-stage sample.
#   x = x-coordinate for location for each site, which is either the
#     x-coordinate for a single-stage sample or the stage two x-coordinate for a
#     two-stage sample.
#   y = y-coordinate for location for each site, which is either the
#     y-coordinate for a single-stage sample or the stage two y-coordinate for a
#     two-stage sample.
#   val = the set of values at which the CDF is estimated.
#   cdfest = the CDF estimate.
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
#      swgt = the size-weight for each site, which is the stage two size-weight 
#         for a two-stage sample.
#   swgt = the size-weight for each site, which is the stage two size-weight for
#     a two-stage sample.
#   swgt1 = the stage one size-weight for each site.
#   vartype = the choice of variance estimator, where "Local" = local mean
#     estimator and "SRS" = SRS estimator.
#   warn.ind = a logical value that indicates whether warning messages were
#     generated, where TRUE = warning messages were generated and FALSE =
#     warning messages were not generated.
#   warn.df = a data frame for storing warning messages.
#   warn.vec = a vector that contains names of the population type, the
#     subpopulation, and an indicator.
# Output:
#   An object in list format composed of a vector named varest, which contains
#   variance estimates, a logical variable named warn,ind, which is the
#   indicator for warning messges, and a data frame named warn.df, which
#   contains warning messages.
# Other Functions Required:
#   localmean.weight - calculate the weighting matrix for the local mean
#     variance estimator
#   localmean.var - calculate the local mean variance estimator
################################################################################

# Assign the function name

   fname <- "cdfvar.size.prop"

# Branch to handle two-stage and single-stage samples

   if(cluster.ind) {

# Begin the section for a two-stage sample

# Calculate additional required values

      m <- length(val)
      cluster <- factor(cluster)
      cluster.levels <- levels(cluster)
      ncluster <- length(cluster.levels)
      z.lst <- split(z, cluster)
      wgt <- wgt*swgt
      if(vartype == "Local") {
         x2.lst <- split(x, cluster)
         y2.lst <- split(y, cluster)
         x1.u <- as.vector(tapply(x1, cluster, unique))
         y1.u <- as.vector(tapply(y1, cluster, unique))
      }
      wgt2.lst <- split(wgt, cluster)
      wgt1 <- wgt1*swgt1
      wgt1.u <- as.vector(tapply(wgt1, cluster, unique))
      tw2 <- (sum(wgt1*wgt))^2
      if(pcfactor.ind) {
         support.lst <- split(support, cluster)
      } else {
         support.lst <- NULL
      }
      var.ind <- sapply(split(cluster, cluster), length) > 1

# Calculate estimates of the total of the stage two sampling unit residuals 
# and the variance of those totals for each stage one sampling unit

      total2est <- matrix(0, ncluster, m)
      var2est <- matrix(0, ncluster, m)
      for(i in 1:ncluster) {

# Calculate the weighted residuals matrix

         n <- length(z.lst[[i]])
         im <- ifelse(matrix(rep(z.lst[[i]], m), nrow = n) <= matrix(rep(val,
            n), nrow = n, byrow = TRUE), 1, 0)
         rm <- (im - matrix(rep(cdfest, n), nrow = n, byrow = TRUE)) *
            matrix(rep(wgt2.lst[[i]], m), nrow = n)

# Calculate total estimates for the stage one sampling unit

         total2est[i,] <- apply(rm, 2, sum)

# Adjust the variance estimator for small sample size

         SRSind <- FALSE
         if(vartype == "Local" && n < 4) {
            warn.ind <- TRUE
            act <- "The simple random sampling variance estimator was used.\n"
            if(stratum.ind) {
               warn <- paste("There are less than four response values for stage one sampling unit ", cluster.levels[i], "\nin stratum ", stratum.level, ", the simple random sampling variance estimator \nwas used to calculate variance of the CDF estimate.\n", sep="")
               warn.df <- rbind(warn.df, data.frame(func=I(fname),
                  subpoptype=warn.vec[1], subpop=warn.vec[2],
                  indicator=warn.vec[3], stratum=I(stratum.level),
                  warning=I(warn), action=I(act)))
            } else {
               warn <- paste("There are less than four response values for stage one sampling unit ", cluster.levels[i], ", \nthe simple random sampling variance estimator was used to calculate variance of \nthe CDF estimate.\n", sep="")
               warn.df <- rbind(warn.df, data.frame(func=I(fname),
                  subpoptype=warn.vec[1], subpop=warn.vec[2],
                  indicator=warn.vec[3], stratum=NA, warning=I(warn),
                  action=I(act)))
            }
            vartype <- "SRS"
            SRSind <- TRUE
         }

# Calculate the population correction factor for the stage two sample

         pcfactor <- ifelse(pcfactor.ind, (stage1size[i] - sum(support.lst[[i]]))/stage1size[i], 1)

# Calculate variance estimates for the stage one sampling unit

         if(var.ind[i]) {
            if(vartype == "Local") {
               weight.lst <- localmean.weight(x2.lst[[i]], y2.lst[[i]], 1/wgt2.lst[[i]])
               var2est[i,] <- pcfactor*apply(rm, 2, localmean.var, weight.lst)
            } else {
               var2est[i,] <- pcfactor*n*apply(rm, 2, var)
               if(SRSind)
                  vartype <- "Local"
            }
         }
      }

# Assign the mean variance to stage one sampling units with a single stage two
# sampling unit
      for(j in 1:m) {
         ind <- var2est[,j] == 0
         if(sum(ind) > 0) {
            var.mean <- mean(var2est[!ind,j])
            var2est[ind,j] <- var.mean
         }
      }

# Adjust the variance estimator for small sample size

      if(vartype == "Local" && ncluster < 4) {
         warn.ind <- TRUE
         act <- "The simple random sampling variance estimator was used.\n"
         if(stratum.ind) {
            warn <- paste("There are less than four stage one sampling units in stratum ", stratum.level, ", \nthe simple random sampling variance estimator was used to calculate variance of \nthe CDF estimate.\n", sep="")
            warn.df <- rbind(warn.df, data.frame(func=I(fname),
               subpoptype=warn.vec[1], subpop=warn.vec[2],
               indicator=warn.vec[3], stratum=I(stratum.level), warning=I(warn),
               action=I(act)))
         } else {
            warn <- paste("There are less than four stage one sampling units, the simple random sampling \nvariance estimator was used to calculate variance of the CDF estimate.\n", sep="")
            warn.df <- rbind(warn.df, data.frame(func=I(fname),
               subpoptype=warn.vec[1], subpop=warn.vec[2],
               indicator=warn.vec[3], stratum=NA, warning=I(warn),
               action=I(act)))
         }
         vartype <- "SRS"
      }

# Calculate the population correction factor for the stage one sample

      pcfactor <- ifelse(pcfactor.ind, (N.cluster - ncluster)/N.cluster, 1)

# Calculate the variance estimate

      if(vartype == "Local") {
         weight.lst <- localmean.weight(x1.u, y1.u, 1/wgt1.u)
         varest <- (pcfactor*apply(total2est * matrix(rep(wgt1.u, m), nrow =
            ncluster), 2, localmean.var, weight.lst) + apply(var2est *
            matrix(rep(wgt1.u, m), nrow = ncluster), 2, sum)) / tw2
      } else {
         varest <- (pcfactor*ncluster*apply(total2est * matrix(rep(wgt1.u, m),
            nrow = ncluster), 2, var) + apply(var2est * matrix(rep(wgt1.u, m),
            nrow = ncluster), 2, sum))/ tw2
      }

# End of section for a two-stage sample

   } else {

# Begin the section for a single-stage sample

# Calculate additional required values

      n <- length(z)
      m <- length(val)
      wgt <- wgt*swgt
      tw2 <- (sum(wgt))^2

# Calculate the weighted residuals matrix

      im <- ifelse(matrix(rep(z, m), nrow = n) <= matrix(rep(val, n), nrow = n,
         byrow = TRUE), 1, 0)
      rm <- (im - matrix(rep(cdfest, n), nrow = n, byrow = TRUE)) *
         matrix(rep(wgt, m), nrow = n)

# Adjust the variance estimator for small sample size

      if(vartype == "Local" && n < 4) {
         warn.ind <- TRUE
         act <- "The simple random sampling variance estimator was used.\n"
         if(stratum.ind) {
            warn <- paste("There are less than four response values in stratum ", stratum.level, ", \nthe simple random sampling variance estimator was used to calculate variance of \nthe CDF estimate.\n", sep="")
            warn.df <- rbind(warn.df, data.frame(func=I(fname),
               subpoptype=warn.vec[1], subpop=warn.vec[2],
               indicator=warn.vec[3], stratum=I(stratum.level), warning=I(warn),
               action=I(act)))
         } else {
            warn <- "\nThere are less than four response values, the simple random sampling variance \nestimator was used to calculate variance of the CDF estimate.\n"
            warn.df <- rbind(warn.df, data.frame(func=I(fname),
               subpoptype=warn.vec[1], subpop=warn.vec[2],
               indicator=warn.vec[3], stratum=NA, warning=I(warn),
               action=I(act)))
          }
         vartype <- "SRS"
      }

# Calculate the population correction factor

      pcfactor <- ifelse(pcfactor.ind, (pcfsize - sum(support))/pcfsize, 1)

# Calculate the variance estimate

       if(vartype == "Local") {
         weight.lst <- localmean.weight(x, y, 1/wgt)
         varest <- pcfactor*apply(rm, 2, localmean.var, weight.lst) / tw2
      } else {
         varest <- pcfactor*n*apply(rm, 2, var) / tw2
      }

# End section for a single-stage sample

   }

# Return the variance estimate, the warning message indicator, and the warn.df
# data frame

   list(varest=varest, warn.ind=warn.ind, warn.df=warn.df)
}
