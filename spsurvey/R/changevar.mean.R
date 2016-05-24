changevar.mean <- function(z1, z2, wgt, x, y, revisitwgt, mean1, mean2,
   stratum.ind, stratum.level, cluster.ind, cluster, wgt1, x1, y1, pcfactor.ind,
   pcfsize, N.cluster, stage1size, support, vartype, warn.ind, warn.df,
   warn.vec) {

################################################################################
# Function: changevar.mean
# Purpose: Calculate the covariance or correlation estimate of the estimated
#          change in means between two probability surveys
# Programmer: Tom Kincaid
# Date: May 23, 2012
# Description:
#   This function calculates estimates of the variance-covariance matrix
#   of the population proportions in a set of intervals (classes).  The set of
#   values defining upper bounds for the classes is supplied to the function.
#   Either the simple random sampling (SRS) variance estimator or the local mean
#   variance estimator is calculated, which is subject to user control.  The
#   simple random sampling variance estimator uses the independent random sample
#   approximation to calculate joint inclusion probabilities.  The function can 
#   accomodate single-stage and two-stage samples.  Finite population and
#   continuous population correction factors can be utilized in variance
#   estimation.
# Arguments:
#   z1 = the response value for each site for survey one.
#   z2 = the response value for each site for survey two.
#   wgt = the final adjusted weight (inverse of the sample inclusion
#     probability) for each site, which is either the weight for a single-stage
#     sample or the stage two weight for a two-stage sample.
#   x = x-coordinate for location for each site, which is either the x-
#     coordinate for a single-stage sample or the stage two x-coordinate for a
#     two-stage sample.
#   y = y-coordinate for location for each site, which is either the y-
#     coordinate for a single-stage sample or the stage two y-coordinate for a
#     two-stage sample.
#   revisitwgt = a logical value that indicates whether each repeat visit site
#     has the same survey design weight in the two surveys, where TRUE = the
#     weight for each repeat visit site is the same and FALSE = the weight for
#     each repeat visit site is not the same.  When this argument is FALSE, all
#     of the repeat visit sites are assigned equal weights when calculating the
#     covariance component of the change estimate standard error.
#   mean1 = the estimated mean for survey one.
#   mean2 = the estimated mean for survey two.
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
# Output:
#   An object in list format composed of a vector named rslt, which contains the
#   covariance or correlation estimate, a logical variable named warn,ind,
#   which is the indicator for warning messges, and a data frame named warn.df,
#   which contains warning messages.
# Other Functions Required:
#   localmean.weight - calculate the weighting matrix for the local mean
#     variance/covariance estimator
#   localmean.cov - calculate the variance/covariance matrix using the local
#     mean estimator
################################################################################

# Assign the function name
fname <- "changevar.mean"

#
# Calculate covariance or correlation using the repeat visit sites
#

# Begin the section for a two-stage sample
if(cluster.ind) {

# Calculate additional required values
   cluster <- factor(cluster)
   cluster.levels <- levels(cluster)
   ncluster <- length(cluster.levels)
   z1.lst <- split(z1, cluster)
   z2.lst <- split(z2, cluster)
   wgt2.lst <- split(wgt, cluster)
   wgt1.u <- as.vector(tapply(wgt1, cluster, unique))
   tw2 <- (sum(wgt1*wgt))^2
   if(vartype == "Local") {
      x2.lst <- split(x, cluster)
      y2.lst <- split(y, cluster)
      x1.u <- as.vector(tapply(x1, cluster, unique))
      y1.u <- as.vector(tapply(y1, cluster, unique))
   }
   if(pcfactor.ind) {
      support.lst <- split(support, cluster)
   } else {
      support.lst <- NULL
   }
   var.ind <- sapply(split(cluster, cluster), length) > 1

# Determine whether the mean could be calculated for both surveys
   if(is.na(mean1) | is.na(mean2)) {
      warn.ind <- TRUE
      act <- "Covariance among the repeat visit sites was not included in calculation of \nthe standard error estimate.\n"
      if(stratum.ind) {
         warn <- paste("The repeat visit sites mean in stratum \"", stratum.level, "\" \ncould not be calculated in one of the surveys.\n", sep="")
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2],
            indicator=warn.vec[3], stratum=I(stratum.level), warning=I(warn),
            action=I(act)))
      } else {
         warn <- paste("The repeat visit sites mean could not be calculated in one of the surveys.\n", sep="")
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2],
            indicator=warn.vec[3], stratum=NA, warning=I(warn),
            action=I(act)))
      }
      rslt <- NA
   } else {

# Calculate estimates of the total of the stage two sampling unit residuals 
# and the variance/covariance of those totals for each stage one sampling unit
      total2est <- matrix(0, ncluster, 2)
      var2est <- matrix(0, ncluster, 4)
      phat <- c(mean1, mean2)
      for(i in 1:ncluster) {

# Calculate the weighted residuals matrix
         n <- length(z1.lst[[i]])
         z1 <- z1.lst[[i]]
         z2 <- z2.lst[[i]]
         rm <- (cbind(z1, z2) - matrix(rep(phat, n), nrow=n, byrow=TRUE)) *
            matrix(rep(wgt2.lst[[i]], 2), nrow=n)

# Calculate total estimates for the stage one sampling unit
         total2est[i,] <- apply(rm, 2, sum)

# Adjust the variance/covariance estimator for small sample size
         SRSind <- FALSE
         if(vartype == "Local" && n < 4) {
            warn.ind <- TRUE
            act <- "The simple random sampling variance estimator was used.\n"
            if(stratum.ind) {
               warn <- paste("There are less than four response values for stage one sampling unit ", cluster.levels[i], "\nin stratum ", stratum.level, ", the simple random sampling variance estimator \nwas used to calculate covariance of the estimate.\n", sep="")
               warn.df <- rbind(warn.df, data.frame(func=I(fname),
                  subpoptype=warn.vec[1], subpop=warn.vec[2],
                  indicator=warn.vec[3], stratum=I(stratum.level),
                  warning=I(warn), action=I(act)))
            } else {
               warn <- paste("There are less than four response values for stage one sampling unit ", cluster.levels[i], ", \nthe simple random sampling variance estimator was used to calculate covariance of \nthe estimate.\n", sep="")
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

# Calculate variance/covariance estimates for the stage one sampling unit
         if(var.ind[i]) {
            if(vartype == "Local") {
               weight.lst <- localmean.weight(x2.lst[[i]], y2.lst[[i]],
                  1/wgt2.lst[[i]])
               var2est[i,] <- as.vector(pcfactor*localmean.cov(rm, weight.lst))
            } else {
               var2est[i,] <- as.vector(pcfactor*n*var(rm))
               if(SRSind)
                  vartype <- "Local"
            }
         }
      }

# Assign the mean variance to stage one sampling units with a single stage two
# sampling unit
      for(j in 1:2) {
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
            warn <- paste("There are less than four stage one sampling units in stratum ", stratum.level, ", \nthe simple random sampling variance estimator was used to calculate covariance of \nthe estimate.\n", sep="")
            warn.df <- rbind(warn.df, data.frame(func=I(fname),
               subpoptype=warn.vec[1], subpop=warn.vec[2],
               indicator=warn.vec[3], stratum=I(stratum.level), warning=I(warn),
               action=I(act)))
         } else {
            warn <- paste("There are less than four stage one sampling units, the simple random sampling \nvariance estimator was used to calculate covariance of the estimate.\n", sep="")
            warn.df <- rbind(warn.df, data.frame(func=I(fname),
               subpoptype=warn.vec[1], subpop=warn.vec[2],
               indicator=warn.vec[3], stratum=NA, warning=I(warn),
               action=I(act)))
         }
         vartype <- "SRS"
      }

# Calculate the population correction factor for the stage one sample
      pcfactor <- ifelse(pcfactor.ind, (N.cluster - ncluster)/N.cluster, 1)

# Calculate the covariance or correlation estimates
      if(vartype == "Local") {
         weight.lst <- localmean.weight(x1.u, y1.u, 1/wgt1.u)
         varest <- (pcfactor*localmean.cov(total2est * matrix(rep(wgt1.u, 2),
            nrow = ncluster), weight.lst) + matrix(apply(var2est *
            matrix(rep(wgt1.u, 4), nrow=ncluster), 2, sum), nrow=2)) / tw2
      } else {
         varest <- (pcfactor*ncluster*var(total2est * matrix(rep(wgt1.u, 2),
            nrow=ncluster)) + matrix(apply(var2est * matrix(rep(wgt1.u, 4),
            nrow=ncluster), 2, sum), nrow=2))/ tw2
      }
      if(revisitwgt) {
         rslt <- varest[1,2]
      } else {
         if(varest[1,1] == 0 | varest[2,2] == 0) {
            warn.ind <- TRUE
            act <- "Covariance among the repeat visit sites was not included in calculation of \nthe standard error estimate.\n"
            if(stratum.ind) {
               warn <- paste("The variance estimate for the repeat visit sites mean in stratum \"", stratum.level, "\" \nwas equal to zero for at least one of the surveys.\n", sep="")
               warn.df <- rbind(warn.df, data.frame(func=I(fname),
                  subpoptype=warn.vec[1], subpop=warn.vec[2],
                  indicator=warn.vec[3], stratum=I(stratum.level),
                  warning=I(warn), action=I(act)))
            } else {
               warn <- paste("The variance estimate for the repeat visit sites mean was equal to zero \nfor at least \none of the surveys.\n", sep="")
               warn.df <- rbind(warn.df, data.frame(func=I(fname),
                  subpoptype=warn.vec[1], subpop=warn.vec[2],
                  indicator=warn.vec[3], stratum=NA, warning=I(warn),
                  action=I(act)))
            }
            rslt <- NA
         }
         rslt <- varest[1,2]/sqrt(varest[1,1]*varest[2,2])
      }
   }

# End the section for a two-stage sample
} else {

# Begin the section for a single-stage sample

# Calculate additional required values
   n <- length(z1)
   tw2 <- (sum(wgt))^2

# Determine whether the mean could be calculated for both surveys
   if(is.na(mean1) | is.na(mean2)) {
      warn.ind <- TRUE
      act <- "Covariance among the repeat visit sites was not included in calculation of \nthe standard error estimate.\n"
      if(stratum.ind) {
         warn <- paste("The repeat visit sites mean in stratum \"", stratum.level, "\" \ncould not be calculated in one of the surveys.\n", sep="")
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2],
            indicator=warn.vec[3], stratum=I(stratum.level), warning=I(warn),
            action=I(act)))
      } else {
         warn <- paste("The repeat visit sites mean could not be calculated in one of the surveys.\n", sep="")
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2],
            indicator=warn.vec[3], stratum=NA, warning=I(warn),
            action=I(act)))
      }
      rslt <- NA
   } else {

# Calculate the weighted residuals matrix
      phat <- c(mean1, mean2)
      rm <- (cbind(z1, z2) - matrix(rep(phat, n), nrow=n, byrow=TRUE)) *
         matrix(rep(wgt, 2), nrow=n)

# Adjust the variance estimator for small sample size
      if(vartype == "Local" && n < 4) {
         warn.ind <- TRUE
         act <- "The simple random sampling variance estimator was used.\n"
         if(stratum.ind) {
            warn <- paste("There are less than four response values in stratum ", stratum.level, ", \nthe simple random sampling variance estimator was used to calculate covariance of \nthe estimate.\n", sep="")
            warn.df <- rbind(warn.df, data.frame(func=I(fname),
               subpoptype=warn.vec[1], subpop=warn.vec[2],
               indicator=warn.vec[3], stratum=I(stratum.level), warning=I(warn),
               action=I(act)))
         } else {
            warn <- "\nThere are less than four response values, the simple random sampling variance \nestimator was used to calculate covariance of the estimate.\n"
            warn.df <- rbind(warn.df, data.frame(func=I(fname),
               subpoptype=warn.vec[1], subpop=warn.vec[2],
               indicator=warn.vec[3], stratum=NA, warning=I(warn),
               action=I(act)))
         }
         vartype <- "SRS"
      }

# Calculate the population correction factor
      pcfactor <- ifelse(pcfactor.ind, (pcfsize - sum(support))/pcfsize, 1)

# Calculate covariance or correlation estimates
      if(vartype == "Local") {
         weight.lst <- localmean.weight(x, y, 1/wgt)
         varest <- pcfactor*localmean.cov(rm, weight.lst) / tw2
      } else {
         varest <- pcfactor*n*var(rm) / tw2
      }
      if(revisitwgt) {
         rslt <- varest[1,2]
      } else {
         if(varest[1,1] == 0 | varest[2,2] == 0) {
            warn.ind <- TRUE
            act <- "Covariance among the repeat visit sites was not included in calculation of \nthe standard error estimate.\n"
            if(stratum.ind) {
               warn <- paste("The variance estimate for the repeat visit sites mean in stratum \"", stratum.level, "\" \nwas equal to zero for at least one of the surveys.\n", sep="")
               warn.df <- rbind(warn.df, data.frame(func=I(fname),
                  subpoptype=warn.vec[1], subpop=warn.vec[2],
                  indicator=warn.vec[3], stratum=I(stratum.level),
                  warning=I(warn), action=I(act)))
            } else {
               warn <- paste("The variance estimate for the repeat visit sites mean was equal to zero \nfor at least \none of the surveys.\n", sep="")
               warn.df <- rbind(warn.df, data.frame(func=I(fname),
                  subpoptype=warn.vec[1], subpop=warn.vec[2],
                  indicator=warn.vec[3], stratum=NA, warning=I(warn),
                  action=I(act)))
            }
            rslt <- NA
         }
         rslt <- varest[1,2]/sqrt(varest[1,1]*varest[2,2])
      }
   }

# End the section for a single-stage sample
}

# Return the covariance or correlation estimate, the warning message indicator,
# and the warn.df data frame
   list(rslt=rslt, warn.ind=warn.ind, warn.df=warn.df)
}
