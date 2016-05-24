change.est <- function(resp.ind, z_1, wgt_1, x_1=NULL, y_1=NULL, repeat_1, z_2,
   wgt_2, x_2=NULL, y_2=NULL, repeat_2, revisitwgt=FALSE, stratum_1=NULL,
   stratum_2=NULL, cluster_1=NULL, cluster_2=NULL, wgt1_1=NULL, x1_1=NULL,
   y1_1=NULL, wgt1_2=NULL, x1_2=NULL, y1_2=NULL, popsize_1=NULL, popsize_2=NULL,
   popcorrect_1=FALSE, pcfsize_1=NULL, N.cluster_1=NULL, stage1size_1=NULL,
   support_1=NULL, popcorrect_2=FALSE, pcfsize_2=NULL, N.cluster_2=NULL,
   stage1size_2=NULL, support_2=NULL, sizeweight_1=FALSE, swgt_1=NULL,
   swgt1_1=NULL, sizeweight_2=FALSE, swgt_2=NULL, swgt1_2=NULL,
   vartype_1="Local", vartype_2="Local", conf=95, check.ind=TRUE, warn.ind=NULL,
   warn.df=NULL, warn.vec=NULL) {

################################################################################
# Function: change.est
# Purpose: Estimate change between two probability surveys
# Programmer: Tom Kincaid
# Date: January 27, 2012
# Last Revised: December 7, 2012
# Description:
#   This function estimates change between two probability surveys.  The
#   function can accommodate both categorical and continuous response variables.
#   For a categorical variable, a change estimate is calculated for each
#   category.  In addition the standard error of the change estimates and
#   confidence bounds are calculated.  Variance estimates are calculated using 
#   either the local mean variance estimator or the simple random sampling (SRS) 
#   variance estimator.  The choice of variance estimator is subject to user 
#   control.  The local mean variance estimator requires the x-coordinate and
#   y-coordinate of each site.  The SRS variance estimator uses the independent 
#   random sample approximation to calculate joint inclusion probabilities.  
#   Confidence bounds are calculated using a Normal distribution multiplier.
#   The function can accommodate a stratified sample.  For a stratified sample, 
#   separate estimates and standard errors are calculated for each stratum,
#   which are used to produce estimates and standard errors for all strata
#   combined.  Strata that contain a single value are removed.  For a stratified
#   sample, when either the size of the resource or the sum of the size-weights
#   of the resource is provided for each stratum, those values are used as
#   stratum weights for calculating the estimates and standard errors for all
#   strata combined.  For a stratified sample when neither the size of the
#   resource nor the sum of the size-weights of the resource is provided for
#   each stratum, estimated values are used as stratum weights for calculating
#   the estimates and standard errors for all strata combined.  The function can
#   accommodate single-stage and two-stage samples for both stratified and
#   unstratified sampling designs.  It is assumed that both surveys employ the
#   same type of survey design.  Finite population and continuous population
#   correction factors can be utilized in variance estimation.  The function
#   checks for compatibility of input values and removes missing values.
# Arguments:
#   resp.ind = a character value that indicates the type of response variable,
#     where "cat" indicates a categorical variable and "cont" indicates a
#     continuous variable.
#   z_1 = response value for each survey one site.
#   wgt_1 = final adjusted weight (inverse of the sample inclusion probability)
#     for each survey one site, which is either the weight for a single- stage
#     sample or the stage two weight for a two-stage sample.
#   x_1 = x-coordinate for location for each survey one site, which is either
#     the x-coordinate for a single-stage sample or the stage two x-coordinate
#     for a two-stage sample.  The default is NULL.
#   y_1 = y-coordinate for location for each survey one site, which is either
#     the y-coordinate for a single-stage sample or the stage two y-coordinate
#     for a two-stage sample.  The default is NULL.
#   repeat_1 = a logical variable that identifies repeat visit sites for survey
#     one.
#   z_2 = response value for each survey two site.
#   wgt_2 = final adjusted weight  for each survey two site.
#   x_2 = x-coordinate for location for each survey two site.  The default is
#     NULL.
#   y_2 = y-coordinate for location for each survey two site.  The default is
#     NULL.
#   repeat_2 = a logical variable that identifies repeat visit sites for survey
#     two.
#   revisitwgt = a logical value that indicates whether each repeat visit site
#     has the same survey design weight in the two surveys, where TRUE = the
#     weight for each repeat visit site is the same and FALSE = the weight for
#     each repeat visit site is not the same.  When this argument is FALSE, all
#     of the repeat visit sites are assigned equal weights when calculating the
#     covariance component of the change estimate standard error.  The default
#     is FALSE.
#   stratum_1 = the stratum for each survey one site.  The default is NULL.
#   stratum_2 = the stratum for each survey two site.  The default is NULL.
#   cluster_1 = the stage one sampling unit (primary sampling unit or cluster)
#     code for each survey one site.  The default is NULL.
#   cluster_2 = the stage one sampling unit (primary sampling unit or cluster)
#     code for each survey two site.  The default is NULL.
#   wgt1_1 = the final adjusted stage one weight for each survey one site.  The
#     default is NULL.
#   x1_1 = the stage one x-coordinate for location for each survey one site.
#     The default is NULL.
#   y1_1 = the stage one y-coordinate for location for each survey one site.
#     The default is NULL.
#   wgt1_2 = the final adjusted stage one weight for each survey two site.  The
#     default is NULL.
#   x1_2 = the stage one x-coordinate for location for each survey two site.
#     The default is NULL.
#   y1_2 = the stage one y-coordinate for location for each survey two site.
#     The default is NULL.
#   popsize_1 = the known size of the survey one resource - the total number of
#     sampling units of a finite resource or the measure of an extensive
#     resource, which is required for calculation of finite and continuous
#     population correction factors for a single-stage sample.  For a stratified
#     sample, this variable also is used to calculate strata weights.  For a
#     stratified sample, this variable must be a vector containing a value for
#     each stratum and must have the names attribute set to identify the stratum
#     codes.  The default is NULL.
#   popsize_2 = the known size of the survey two resource.  The default is NULL.
#   popcorrect_1 = a logical value that indicates whether finite or continuous
#     population correction factors should be employed during variance
#     estimation for survey one, where TRUE = use the correction factors and
#     FALSE = do not use the correction factors.  The default is FALSE.
#   pcfsize_1 = size of the survey one resource, which is required for
#     calculation of finite and continuous population correction factors for a
#     single-stage sample.  For a stratified sample this argument must be a
#     vector containing a value for each stratum and must have the names
#     attribute set to identify the stratum codes.  The default is NULL.
#   N.cluster_1 = the number of stage one sampling units in the survey one
#     resource, which is required for calculation of finite and continuous
#     population correction factors for a two-stage sample.  For a stratified
#     sample this variable must be a vector containing a value for each stratum
#     and must have the names attribute set to identify the stratum codes.  The
#     default is NULL.
#   stage1size_1 = size of the stage one sampling units of a two-stage sample
#     for survey one, which is required for calculation of finite and continuous
#     population correction factors for a two-stage sample and must have the
#     names attribute set to identify the stage one sampling unit codes.  For a
#     stratified sample, the names attribute must be set to identify both
#     stratum codes and stage one sampling unit codes using a convention where
#     the two codes are separated by the & symbol, e.g., "Stratum 1&Cluster 1".
#     The default is NULL.
#   support_1 = the support value for each survey one site - the value one (1)
#     for a site from a finite resource or the measure of the sampling unit
#     associated with a site from an extensive resource, which is required for
#     calculation of finite and continuous population correction factors.  The
#     default is NULL.
#   popcorrect_2 = a logical value that indicates whether finite or continuous
#     population correction factors should be employed during variance
#     estimation for survey two.  The default is FALSE.
#   pcfsize_2 = size of the survey two resource.  The default is NULL.
#   N.cluster_2 = the number of stage one sampling units in the survey two
#     resource.  The default is NULL.
#   stage1size_2 = size of the stage one sampling units of a two-stage survey
#     for survey two.  The default is NULL.
#   support_2 = the support value for each survey two site.  The default is
#     NULL.
#   sizeweight_1 = a logical value that indicates whether size-weights should
#     be used in the analysis for survey one, where TRUE = use the size-weights
#     and FALSE = do not use the size-weights.  The default is FALSE.
#   swgt_1 = the size-weight for each survey one site, which is the stage two
#     size-weight for a two-stage sample.  The default is NULL.
#   swgt1_1 = the stage one size-weight for each survey one site.  The default
#     is NULL.
#   sizeweight_2 = a logical value that indicates whether size-weights should
#     be used in the analysis for survey two.  The default is FALSE.
#   swgt_2 = the size-weight for each survey two site.  The default is NULL.
#   swgt1_2 = the stage one size-weight for each survey two site.  The default
#     is NULL.
#   vartype_1 = the choice of variance estimator for survey one, where "Local" =
#     local mean estimator and "SRS" = SRS estimator.  The default is "Local".
#   vartype_2 = the choice of variance estimator for survey two.  The default is
#     "Local".
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
#   If the function was called by the change.analysis function, then output is
#   an object in list format composed of the Results data frame, which contains
#   estimates and confidence bounds, and the warn.df data frame, which contains
#   warning messages.  If the function was called directly, then output is the
#   Results data frame.
# Other Functions Required:
#   input.check - check input values for errors, consistency, and compatibility
#     with analytical functions
#   wnas - remove missing values
#   vecprint - takes an input vector and outputs a character string with line
#     breaks inserted
#   category.est - estimate proportion (expressed as percent) and size of a
#     resource in each of a set of categories
#   changevar.prop - calculate covariance or correlation estimates of the
#     estimated change in class proportion estimates between two probability
#     surveys
#   total.est - estimate the population total, mean, variance, and standard
#     deviation of a response variable
#   changevar.mean - calculate the covariance or correlation estimate of the
#     estimated change in means between two probability surveys
# Examples:
#   z_1 <- sample(c("Good","Fair","Poor"), 100, replace=TRUE)
#   z_2 <- sample(c("Good","Fair","Poor"), 100, replace=TRUE)
#   wgt_1 <- runif(100, 10, 100)
#   wgt_2 <- runif(100, 10, 100)
#   repeat_1 <- rep(c(TRUE,FALSE), c(20,80))
#   repeat_2 <- rep(c(TRUE,FALSE), c(20,80))
#   change.est(resp.ind="cat", z_1=z_1, wgt_1=wgt_1, repeat_1=repeat_1,
#      z_2=z_2, wgt_2=wgt_2, repeat_2=repeat_2, vartype_1="SRS",
#      vartype_2="SRS")
#
#   z_1 <- rnorm(100, 10,10)
#   z_2 <- rnorm(100, 12, 10)
#   change.est(resp.ind="cont", z_1=z_1, wgt_1=wgt_1, repeat_1=repeat_1,
#      z_2=z_2, wgt_2=wgt_2, repeat_2=repeat_2, vartype_1="SRS",
#      vartype_2="SRS")
################################################################################

# As necessary, create a data frame for warning messages
if(is.null(warn.ind)) {
   warn.ind <- FALSE
   warn.df <- NULL
   warn.vec <- c(0, NA, NA)
}
fname <- "change.est"

# Check the logical variables for repeat visit sites
if(check.ind) {
   if(!is.logical(repeat_1))
      stop("\nThe repeat_1 argument must be a logical variable.")
   if(!is.logical(repeat_2))
      stop("\nThe repeat_2 argument must be a logical variable.")
   n1 <- sum(repeat_1)
   n2 <- sum(repeat_2)
   if(n1 != n2)
      stop(paste("\nThe number of repeat visit sites for survey one, ", n1, ", does not equal the number \nof repeat visit sites for survey two, ", n2, ".", sep=""))
}

#
# Begin the section for a categorical variable
#

if(resp.ind == "cat") {

# Calculate estimates for all sites from survey one
   temp <- category.est(catvar=z_1, wgt=wgt_1, x=x_1, y=y_1, stratum=stratum_1,
      cluster=cluster_1, wgt1=wgt1_1, x1=x1_1, y1=y1_1, popsize=popsize_1,
      popcorrect=popcorrect_1, pcfsize=pcfsize_1, N.cluster=N.cluster_1,
      stage1size=stage1size_1, support=support_1, sizeweight=sizeweight_1,
      swgt=swgt_1, swgt1=swgt1_1, vartype=vartype_1, conf=conf,
      check.ind=check.ind, warn.ind=warn.ind, warn.df=warn.df,
      warn.vec=warn.vec)
   temp.cat_1 <- temp$Results[temp$Results$Category != "Total",]
   tw_1 <- temp$Results$Estimate.U[temp$Results$Category == "Total"]
   warn.ind <- temp$warn.ind
   warn.df <- temp$warn.df

# Calculate estimates for all sites from survey two
   temp <- category.est(catvar=z_2, wgt=wgt_2, x=x_2, y=y_2, stratum=stratum_2,
      cluster=cluster_2, wgt1=wgt1_2, x1=x1_2, y1=y1_2, popsize=popsize_2,
      popcorrect=popcorrect_2, pcfsize=pcfsize_2, N.cluster=N.cluster_2,
      stage1size=stage1size_2, support=support_2, sizeweight=sizeweight_2,
      swgt=swgt_2, swgt1=swgt1_2, vartype=vartype_2, conf=conf,
      check.ind=check.ind, warn.ind=warn.ind, warn.df=warn.df,
      warn.vec=warn.vec)
   temp.cat_2 <- temp$Results[temp$Results$Category != "Total",]
   tw_2 <- temp$Results$Estimate.U[temp$Results$Category == "Total"]
   warn.ind <- temp$warn.ind
   warn.df <- temp$warn.df

# Merge results for the two surveys
   Results <- merge(temp.cat_1, temp.cat_2, by="Category", suffix=c("_1", "_2"),
      all=TRUE, sort=FALSE)

# Calculate the change estimates
   Results$DiffEst.P <- (Results$Estimate.P_2 - Results$Estimate.P_1)/100
   Results$DiffEst.U <- Results$Estimate.U_2 - Results$Estimate.U_1

# Express the standard error estimates for the two surveys on the proportion
# scale
   Results$StdError.P_1 <- Results$StdError.P_1/100
   Results$StdError.P_2 <- Results$StdError.P_2/100

# Calculate confidence bound multiplier
   mult <- qnorm(0.5 + (conf/100)/2)

#
# Calculate standard error of the change estimates
#

# Section for surveys with no repeat visit sites
   if(sum(repeat_1) == 0) {
      Results$StdError.P <- sqrt(Results$StdError.P_1^2 + Results$StdError.P_2^2)
      Results$StdError.U <- sqrt(Results$StdError.U_1^2 + Results$StdError.U_2^2)
      eval(parse(text=paste("Results$LCB", conf, "Pct.P <- 100*pmax(Results$DiffEst.P - mult*Results$StdError.P, -1)", sep="")))
      eval(parse(text=paste("Results$UCB", conf, "Pct.P <- 100*pmin(Results$DiffEst.P + mult*Results$StdError.P, 1)", sep="")))
      Results$DiffEst.P <- 100*Results$DiffEst.P
      Results$StdError.P <- 100*Results$StdError.P
      Results$StdError.P_1 <- 100*Results$StdError.P_1
      Results$StdError.P_2 <- 100*Results$StdError.P_2
      if(!is.null(popsize_1)) {
         if(!is.null(popsize_2)) {
            eval(parse(text=paste("Results$LCB", conf, "Pct.U <- pmax(Results$DiffEst.U - mult*Results$StdError.U, -popsize_1)", sep="")))
            eval(parse(text=paste("Results$UCB", conf, "Pct.U <- pmin(Results$DiffEst.U + mult*Results$StdError.U, popsize_2)", sep="")))
         } else {
            eval(parse(text=paste("Results$LCB", conf, "Pct.U <- pmax(Results$DiffEst.U - mult*Results$StdError.U, -popsize_1)", sep="")))
            eval(parse(text=paste("Results$UCB", conf, "Pct.U <- Results$DiffEst.U + mult*Results$StdError.U", sep="")))
         }
      } else {
         if(!is.null(popsize_2)) {
            eval(parse(text=paste("Results$LCB", conf, "Pct.U <- Results$DiffEst.U - mult*Results$StdError.U", sep="")))
            eval(parse(text=paste("Results$UCB", conf, "Pct.U <- pmin(Results$DiffEst.U + mult*Results$StdError.U, popsize_2)", sep="")))
         } else {
            eval(parse(text=paste("Results$LCB", conf, "Pct.U <- Results$DiffEst.U - mult*Results$StdError.U", sep="")))
            eval(parse(text=paste("Results$UCB", conf, "Pct.U <- Results$DiffEst.U + mult*Results$StdError.U", sep="")))
         }
      }

   } else {

# Section for surveys with repeat visit sites

# Assign values for the categorical variables
   catvar <- 1:sum(repeat_1)
   catvar_1 <- z_1[repeat_1]
   catvar_2 <- z_2[repeat_2]
   catvar[is.na(catvar_1) | is.na(catvar_2)] <- NA

# Assign values for survey design variables using the survey one sites
   if(revisitwgt) {
      wgt <- wgt_1[repeat_1]
   } else {
      wgt <- rep(1, length(catvar_1))
   }
   x <- x_1[repeat_1]
   y <- y_1[repeat_1] 
   stratum <- stratum_1[repeat_1]
   cluster <- cluster_1[repeat_1]
   wgt1 <- wgt1_1[repeat_1]
   x1 <- x1_1[repeat_1]
   y1 <- y1_1[repeat_1]
   popsize <- popsize_1
   pcfsize <- pcfsize_1
   N.cluster <- N.cluster_1
   stage1size <- stage1size_1
   support <- support_1[repeat_1]

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
   pcfactor.ind <- popcorrect_1

# Assign the value of sizeweight to the indicator variable for use of size
# weights
   swgt.ind <- sizeweight_1

# Assign the value of vartype to the character variable for type of variance
# estimator
   vartype <- vartype_1

# Remove missing values
   if(vartype == "Local") {
      if(swgt.ind) {
         if(stratum.ind) {
            if(cluster.ind)
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y,
                  stratum=stratum, cluster=cluster, wgt1=wgt1, x1=x1, y1=y1,
                  swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y,
                  stratum=stratum, swgt=swgt))
         } else {
            if(cluster.ind)
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y,
                  cluster=cluster, wgt1=wgt1, x1=x1, y1=y1, swgt=swgt,
                  swgt1=swgt1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y, swgt=swgt))
         }
         catvar <- temp$catvar
         catvar_1 <- catvar_1[catvar]
         catvar_2 <- catvar_2[catvar]
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
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y,
                  stratum=stratum, cluster=cluster, wgt1=wgt1, x1=x1, y1=y1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y,
                  stratum=stratum))
         } else {
            if(cluster.ind)
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y,
                  cluster=cluster, wgt1=wgt1, x1=x1, y1=y1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt, x=x, y=y))
         }
         catvar <- temp$catvar
         catvar_1 <- catvar_1[catvar]
         catvar_2 <- catvar_2[catvar]
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
               temp <- wnas(list(catvar=catvar, wgt=wgt, stratum=stratum,
                  cluster=cluster, wgt1=wgt1, swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt, stratum=stratum,
                  swgt=swgt))
         } else {
            if(cluster.ind)
               temp <- wnas(list(catvar=catvar, wgt=wgt, cluster=cluster,
                  wgt1=wgt1, swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt, swgt=swgt))
         }
         catvar <- temp$catvar
         catvar_1 <- catvar_1[catvar]
         catvar_2 <- catvar_2[catvar]
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
               temp <- wnas(list(catvar=catvar, wgt=wgt, stratum=stratum,
                  cluster=cluster, wgt1=wgt1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt, stratum=stratum))
         } else {
            if(cluster.ind)
               temp <- wnas(list(catvar=catvar, wgt=wgt, cluster=cluster,
                  wgt1=wgt1))
            else
               temp <- wnas(list(catvar=catvar, wgt=wgt))
         }
         catvar <- temp$catvar
         catvar_1 <- catvar_1[catvar]
         catvar_2 <- catvar_2[catvar]
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
      if(nstrata > 0 & nstrata < nstrata.old) {
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
         warn.df <- rbind(warn.df, data.frame(func=I(fname), subpoptype=warn.vec[1],
            subpop=warn.vec[2], indicator=warn.vec[3], stratum=NA, warning=I(warn),
            action=I(act)))
      }
   }

# Assign levels of the categorical variables
   catvar.levels <- Results$Category
   nlevels <- length(catvar.levels)

# Calculate additional required values
   if(length(catvar_1) > 1) {
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
   }

# Branch to handle stratified and unstratified data

   if(stratum.ind) {

# Begin the section for stratified data

# Create the vector of covariance or correlation estimates for all strata
# combined
   rslt.p <- rep(NA, nlevels)
   rslt.u <- rep(NA, nlevels)

# Check whether the vectors of categorical variable values for revisit sites
# are empty or contain a single value
   if(length(catvar_1) <= 1) {
      warn.ind <- TRUE
      act <- "Covariance among the revisited sites was not included in calculation of \nthe standard error estimate.\n"
      warn <- paste("The number of nonmissing repeat visit sites was less than two in one of the \nsurveys.\n", sep="")
      warn.df <- rbind(warn.df, data.frame(func=I(fname),
         subpoptype=warn.vec[1], subpop=warn.vec[2], indicator=warn.vec[3],
         stratum=NA, warning=I(warn), action=I(act)))

# Begin section for nonempty vectors of categorical variables values for revisit
# sites
   } else {

# Begin the loop for individual strata

   for(i in 1:nstrata) {

# Check whether the vectors of categorical variable values for revisit sites
# are empty or contain a single value for a stratum
      stratum.i <- stratum == stratum.levels[i]
      if(length(catvar_1[stratum.i]) <= 1) {
         warn.ind <- TRUE
         act <- "Due to insufficient number of sites, the stratum was not included in \ncalculation of covariance among the revisited sites.\n"
         warn <- paste("The number of nonmissing repeat visit sites  in one of the surveys was less \nthan two for stratum \"", stratum.levels[i], "\".\n", sep="")
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2],
            indicator=warn.vec[3], stratum=I(stratum.levels[i]),
            warning=I(warn), action=I(act)))

# Begin section for nonempty vectors of categorical variables values for revisit
# sites for a stratum
      } else {
# Calculate proportion estimates
         z1 <- factor(catvar_1[stratum.i], levels=catvar.levels)
         z2 <- factor(catvar_2[stratum.i], levels=catvar.levels)
         m <- length(catvar.levels)
         if(swgt.ind) {
            if(cluster.ind) {
               w2 <- wgt[stratum.i]*swgt[stratum.i]
               w1 <- wgt1[stratum.i]*swgt1[stratum.i]
               prop1 <- tapply(w2*w1, z1, sum) / popsize.hat[i]
               prop2 <- tapply(w2*w1, z2, sum) / popsize.hat[i]
            } else {
               w <- wgt[stratum.i]*swgt[stratum.i]
               prop1 <- tapply(w, z1, sum) / popsize.hat[i]
               prop2 <- tapply(w, z2, sum) / popsize.hat[i]
            }
         } else {
            if(cluster.ind) {
               w2 <- wgt[stratum.i]
               w1 <- wgt1[stratum.i]
               prop1 <- tapply(w2*w1, z1, sum) / popsize.hat[i]
               prop2 <- tapply(w2*w1, z2, sum) / popsize.hat[i]
            } else {
               w <- wgt[stratum.i]
               prop1 <- tapply(w, z1, sum) / popsize.hat[i]
               prop2 <- tapply(w, z2, sum) / popsize.hat[i]
            }
         }

# Calculate covariance or correlation estimates
         if(cluster.ind) {
            temp <- changevar.prop(catvar.levels, z1, z2, w2, x[stratum.i],
               y[stratum.i], revisitwgt, prop1, prop2, stratum.ind,
               stratum.levels[i], cluster.ind, cluster[stratum.i], w1,
               x1[stratum.i], y1[stratum.i], pcfactor.ind, NULL, N.cluster[i],
               stage1size[[i]], support[stratum.i], vartype, warn.ind, warn.df,
               warn.vec)
         } else {
            temp <- changevar.prop(catvar.levels, z1, z2, w, x[stratum.i],
               y[stratum.i], revisitwgt, prop1, prop2, stratum.ind,
               stratum.levels[i], cluster.ind, pcfactor.ind=pcfactor.ind,
               pcfsize=pcfsize[i], support=support[stratum.i], vartype=vartype,
               warn.ind=warn.ind, warn.df=warn.df, warn.vec=warn.vec)
         }
         correst <- temp$rslt
         warn.ind <- temp$warn.ind
         warn.df <- temp$warn.df

# Add estimates to the vector for all strata combined
         if(!is.null(popsize)) {
            rslt.p[!is.na(correst)] <- rslt.p[!is.na(correst)] +
               (popsize[i]/sum.popsize)*correst[!is.na(correst)]
         } else {
            rslt.p[!is.na(correst)] <- rslt.p[!is.na(correst)] +
               (popsize.hat[i]/sum.popsize.hat)*correst[!is.na(correst)]
         }

# Estimate the size of each category
         if(!is.null(popsize)) {
            size1 <- popsize[i]*prop1
            size2 <- popsize[i]*prop2
         } else {
            size1 <- popsize.hat[i]*prop1
            size2 <- popsize.hat[i]*prop2
         }

# Calculate covariance or correlation estimates
         if(cluster.ind) {
            temp <- changevar.size(catvar.levels, z1, z2, w2, x[stratum.i],
               y[stratum.i], revisitwgt, size1, size2, stratum.ind,
               stratum.levels[i], cluster.ind, cluster[stratum.i], w1,
               x1[stratum.i], y1[stratum.i], pcfactor.ind, NULL, N.cluster[i],
               stage1size[[i]], support[stratum.i], vartype, warn.ind, warn.df,
               warn.vec)
         } else {
            temp <- changevar.size(catvar.levels, z1, z2, w, x[stratum.i],
               y[stratum.i], revisitwgt, size1, size2, stratum.ind,
               stratum.levels[i], cluster.ind, pcfactor.ind=pcfactor.ind,
               pcfsize=pcfsize[i], support=support[stratum.i], vartype=vartype,
               warn.ind=warn.ind, warn.df=warn.df, warn.vec=warn.vec)
         }
         correst <- temp$rslt
         warn.ind <- temp$warn.ind
         warn.df <- temp$warn.df

# Add estimates to the vector for all strata combined
         if(!is.null(popsize)) {
            rslt.u[!is.na(correst)] <- rslt.u[!is.na(correst)] +
               correst[!is.na(correst)]
         } else {
            rslt.u[!is.na(correst)] <- rslt.u[!is.na(correst)] +
               correst[!is.na(correst)]
         }

# End the section for nonempty vectors of categorical variables values for
# revisit sites for a stratum
      }

# End the loop for individual strata
   }

# End the section for nonempty vectors of categorical variables values for
# revisit sites
   }

# End the section for stratified data
   } else {

# Begin the section for unstratified data

# Check whether the vectors of categorical variable values for revisit sites
# are empty or contain a single value
   if(length(catvar_1) <= 1) {
      rslt.p <- rep(NA, nlevels)
      rslt.u <- rep(NA, nlevels)
      warn.ind <- TRUE
      act <- "Covariance among the revisited sites was not included in calculation of \nthe standard error estimate.\n"
      warn <- paste("The number of nonmissing repeat visit sites was less than two in one of the \nsurveys.\n", sep="")
      warn.df <- rbind(warn.df, data.frame(func=I(fname),
         subpoptype=warn.vec[1], subpop=warn.vec[2], indicator=warn.vec[3],
         stratum=NA, warning=I(warn), action=I(act)))

# Begin section for nonempty vectors of categorical variables values for revisit
# sites
   } else {

# Calculate proportion estimates
      z1 <- factor(catvar_1, levels=catvar.levels)
      z2 <- factor(catvar_2, levels=catvar.levels)
      if(swgt.ind) {
         if(cluster.ind) {
            w2 <- wgt*swgt
            w1 <- wgt1*swgt1
            prop1 <- tapply(w2*w1, z1, sum) / popsize.hat
            prop2 <- tapply(w2*w1, z2, sum) / popsize.hat
         } else {
            w <- wgt*swgt
            prop1 <- tapply(w, z1, sum) / popsize.hat
            prop2 <- tapply(w, z2, sum) / popsize.hat
         }
      } else {
         if(cluster.ind) {
            w2 <- wgt
            w1 <- wgt1
            prop1 <- tapply(w2*w1, z1, sum) / popsize.hat
            prop2 <- tapply(w2*w1, z2, sum) / popsize.hat
         } else {
            w <- wgt
            prop1 <- tapply(w, z1, sum) / popsize.hat
            prop2 <- tapply(w, z2, sum) / popsize.hat
         }
      }

# Calculate covariance or correlation estimates
      if(cluster.ind) {
         temp <- changevar.prop(catvar.levels, z1, z2, w2, x, y, revisitwgt,
            prop1, prop2, stratum.ind, NULL, cluster.ind, cluster, w1, x1, y1,
            pcfactor.ind, NULL, N.cluster, stage1size, support, vartype,
            warn.ind, warn.df, warn.vec)
      } else {
         temp <- changevar.prop(catvar.levels, z1, z2, w, x, y, revisitwgt,
            prop1, prop2, stratum.ind, NULL, cluster.ind,
            pcfactor.ind=pcfactor.ind, pcfsize=pcfsize, support=support,
            vartype=vartype, warn.ind=warn.ind, warn.df=warn.df,
            warn.vec=warn.vec)
      }
      rslt.p <- temp$rslt
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df

# Estimate the size of each category
      if(!is.null(popsize)) {
         size1 <- popsize*prop1
         size2 <- popsize*prop2
      } else {
         size1 <- popsize.hat*prop1
         size2 <- popsize.hat*prop2
      }

# Calculate covariance or correlation estimates
      if(cluster.ind) {
         temp <- changevar.size(catvar.levels, z1, z2, w2, x, y, revisitwgt,
            size1, size2, stratum.ind, NULL, cluster.ind, cluster, w1, x1, y1,
            pcfactor.ind, NULL, N.cluster, stage1size, support, vartype,
            warn.ind, warn.df, warn.vec)
      } else {
         temp <- changevar.size(catvar.levels, z1, z2, w, x, y, revisitwgt,
            size1, size2, stratum.ind, NULL, cluster.ind,
            pcfactor.ind=pcfactor.ind, pcfsize=pcfsize, support=support,
            vartype=vartype, warn.ind=warn.ind, warn.df=warn.df,
            warn.vec=warn.vec)
      }
      rslt.u <- temp$rslt
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df

# End the section for nonempty vectors of categorical variables values for
# revisit sites
   }

# End the section for unstratified data
   }

# Calculate standard errors
   Results$StdError.P <- rep(NA, nlevels)
   Results$StdError.U <- rep(NA, nlevels)
   ind <- is.na(rslt.p)
   Results$StdError.P[ind] <- sqrt(Results$StdError.P_1[ind]^2 +
      Results$StdError.P_2[ind]^2)
   Results$StdError.U[ind] <- sqrt(Results$StdError.U_1[ind]^2 +
      Results$StdError.U_2[ind]^2)
   if(any(!ind)) {
      tw_1r <- sum(wgt_1[repeat_1])
      tw_2r <- sum(wgt_2[repeat_2])
      if(revisitwgt) {
         temp <- Results$StdError.P_1^2 + Results$StdError.P_2^2 -
            ((2*tw_1r*tw_2r)/(tw_1*tw_2))*rslt.p
         ind <- !is.na(rslt.p) & temp <= 0
         Results$StdError.P[ind] <- sqrt(Results$StdError.P_1[ind]^2 +
            Results$StdError.P_2[ind]^2)
         ind <- !is.na(rslt.p) & temp > 0
         Results$StdError.P[ind] <- sqrt(temp[ind])
         temp <- Results$StdError.U_1^2 + Results$StdError.U_2^2 -2*rslt.u
         ind <- !is.na(rslt.u) & temp <= 0
         Results$StdError.U[ind] <- sqrt(Results$StdError.U_1[ind]^2 +
            Results$StdError.U_2[ind]^2)
         Results$StdError.U[ind] <- sqrt(Results$StdError.U_1[ind]^2 +
            Results$StdError.U_2[ind]^2)
         ind <- !is.na(rslt.u) & temp > 0
         Results$StdError.U[ind] <- sqrt(temp[ind])
      } else {
         temp <- category.est(catvar=z_1[repeat_1], wgt=wgt_1[repeat_1],
            x=x_1[repeat_1], y=y_1[repeat_1], stratum=stratum_1[repeat_1],
            cluster=cluster_1[repeat_1], wgt1=wgt1_1[repeat_1],
            x1=x1_1[repeat_1], y1=y1_1[repeat_1], popsize=popsize_1,
            popcorrect=popcorrect_1, pcfsize=pcfsize_1, N.cluster=N.cluster_1,
            stage1size=stage1size_1, support=support_1[repeat_1],
            sizeweight=sizeweight_1, swgt=swgt_1[repeat_1],
            swgt1=swgt1_1[repeat_1], vartype=vartype_1, conf=conf,
            check.ind=check.ind, warn.ind=warn.ind, warn.df=warn.df,
            warn.vec=warn.vec)
         se_1.p <- rep(NA, nlevels)
         se_1.u <- rep(NA, nlevels)
         temp$Results <- temp$Results[temp$Results$Category != "Total",]
         ind <- match(temp$Results$Category, catvar.levels, nomatch=0)
         se_1.p[ind] <- temp$Results$StdError.P/100
         se_1.u[ind] <- temp$Results$StdError.U
         temp <- category.est(catvar=z_2[repeat_2], wgt=wgt_2[repeat_2],
            x=x_2[repeat_2], y=y_2[repeat_2], stratum=stratum_2[repeat_2],
            cluster=cluster_2[repeat_2], wgt1=wgt1_2[repeat_2],
            x1=x1_2[repeat_2], y1=y1_2[repeat_2], popsize=popsize_2,
            popcorrect=popcorrect_2, pcfsize=pcfsize_2, N.cluster=N.cluster_2,
            stage1size=stage1size_2, support=support_2[repeat_2],
            sizeweight=sizeweight_2, swgt=swgt_2[repeat_2],
            swgt1=swgt1_2[repeat_2], vartype=vartype_2, conf=conf,
            check.ind=check.ind, warn.ind=warn.ind, warn.df=warn.df,
            warn.vec=warn.vec)
         se_2.p <- rep(NA, nlevels)
         se_2.u <- rep(NA, nlevels)
         temp$Results <- temp$Results[temp$Results$Category != "Total",]
         ind <- match(temp$Results$Category, catvar.levels, nomatch=0)
         se_2.p[ind] <- temp$Results$StdError.P/100
         se_2.u[ind] <- temp$Results$StdError.U
         covest <- rslt.p * se_1.p * se_2.p
         temp <- Results$StdError.P_1^2 + Results$StdError.P_2^2 -
            ((2*tw_1r*tw_2r)/(tw_1*tw_2))*covest
         ind <- !is.na(rslt.p) & temp <= 0
         Results$StdError.P[ind] <- sqrt(Results$StdError.P_1[ind]^2 +
            Results$StdError.P_2[ind]^2)
         ind <- !is.na(rslt.p) & temp > 0
         Results$StdError.P[ind] <- sqrt(temp[ind])
         covest <- rslt.u * se_1.u * se_2.u
         temp <- Results$StdError.U_1^2 + Results$StdError.U_2^2 -2*covest
         ind <- !is.na(rslt.u) & temp <= 0
         Results$StdError.U[ind] <- sqrt(Results$StdError.U_1[ind]^2 +
            Results$StdError.U_2[ind]^2)
         ind <- !is.na(rslt.u) & temp > 0
         Results$StdError.U[ind] <- sqrt(temp[ind])
      }
   }

# Calculate confidence bounds
   eval(parse(text=paste("Results$LCB", conf, "Pct.P <- 100*pmax(Results$DiffEst.P - mult*Results$StdError.P, -1)", sep="")))
   eval(parse(text=paste("Results$UCB", conf, "Pct.P <- 100*pmin(Results$DiffEst.P + mult*Results$StdError.P, 1)", sep="")))
   Results$DiffEst.P <- 100*Results$DiffEst.P
   Results$StdError.P <- 100*Results$StdError.P
   Results$StdError.P_1 <- 100*Results$StdError.P_1
   Results$StdError.P_2 <- 100*Results$StdError.P_2
   if(!is.null(popsize_1)) {
      if(!is.null(popsize_2)) {
         eval(parse(text=paste("Results$LCB", conf, "Pct.U <- pmax(Results$DiffEst.U - mult*Results$StdError.U, -popsize_1)", sep="")))
         eval(parse(text=paste("Results$UCB", conf, "Pct.U <- pmin(Results$DiffEst.U + mult*Results$StdError.U, popsize_2)", sep="")))
      } else {
         eval(parse(text=paste("Results$LCB", conf, "Pct.U <- pmax(Results$DiffEst.U - mult*Results$StdError.U, -popsize_1)", sep="")))
         eval(parse(text=paste("Results$UCB", conf, "Pct.U <- Results$DiffEst.U + mult*Results$StdError.U", sep="")))
      }
   } else {
      if(!is.null(popsize_2)) {
         eval(parse(text=paste("Results$LCB", conf, "Pct.U <- Results$DiffEst.U - mult*Results$StdError.U", sep="")))
         eval(parse(text=paste("Results$UCB", conf, "Pct.U <- pmin(Results$DiffEst.U + mult*Results$StdError.U, popsize_2)", sep="")))
      } else {
         eval(parse(text=paste("Results$LCB", conf, "Pct.U <- Results$DiffEst.U - mult*Results$StdError.U", sep="")))
         eval(parse(text=paste("Results$UCB", conf, "Pct.U <- Results$DiffEst.U + mult*Results$StdError.U", sep="")))
      }
   }

# End the section for surveys with repeat visit sites
   }

# Rearrange columns of the Results data frame
   Results <- Results[,c(1, 20, 22, 24:25, 21, 23, 26:27, 2:19)]

#
# End the section for a categorical variable
#

   } else if(resp.ind == "cont") {

#
# Begin the section for a continuous variable
#

# Calculate estimates for all sites from survey one
   temp <- total.est(z=z_1, wgt=wgt_1, x=x_1, y=y_1, stratum=stratum_1,
      cluster=cluster_1, wgt1=wgt1_1, x1=x1_1, y1=y1_1, popsize=popsize_1,
      popcorrect=popcorrect_1, pcfsize=pcfsize_1, N.cluster=N.cluster_1,
      stage1size=stage1size_1, support=support_1, sizeweight=sizeweight_1,
      swgt=swgt_1, swgt1=swgt1_1, vartype=vartype_1, conf=conf,
      check.ind=check.ind, warn.ind=warn.ind, warn.df=warn.df,
      warn.vec=warn.vec)
   temp.cont_1 <- temp$Results[temp$Results$Statistic == "Mean",]
   warn.ind <- temp$warn.ind
   warn.df <- temp$warn.df
   tw_1 <- sum(wgt_1)

# Calculate estimates for all sites from survey two
   temp <- total.est(z=z_2, wgt=wgt_2, x=x_2, y=y_2, stratum=stratum_2,
      cluster=cluster_2, wgt1=wgt1_2, x1=x1_2, y1=y1_2, popsize=popsize_2,
      popcorrect=popcorrect_2, pcfsize=pcfsize_2, N.cluster=N.cluster_2,
      stage1size=stage1size_2, support=support_2, sizeweight=sizeweight_2,
      swgt=swgt_2, swgt1=swgt1_2, vartype=vartype_2, conf=conf,
      check.ind=check.ind, warn.ind=warn.ind, warn.df=warn.df,
      warn.vec=warn.vec)
   temp.cont_2 <- temp$Results[temp$Results$Statistic == "Mean",1:6]
   warn.ind <- temp$warn.ind
   warn.df <- temp$warn.df
   tw_2 <- sum(wgt_2)

# Merge results for the two surveys
   Results <- merge(temp.cont_1, temp.cont_2, by="Statistic", suffix=c("_1", "_2"))

# Calculate the change estimates
   Results$DiffEst <- Results$Estimate_2 - Results$Estimate_1

# Calculate confidence bound multiplier
   mult <- qnorm(0.5 + (conf/100)/2)

#
# Calculate standard error of the change estimates
#

# Section for surveys with no repeat visit sites
   if(sum(repeat_1) == 0) {
      Results$StdError <- sqrt(Results$StdError_1^2 + Results$StdError_2^2)
      eval(parse(text=paste("Results$LCB", conf, "Pct <- Results$DiffEst - mult*Results$StdError", sep="")))
      eval(parse(text=paste("Results$UCB", conf, "Pct <- Results$DiffEst + mult*Results$StdError", sep="")))
   } else {

# Section for surveys with repeat visit sites

# Assign values for the continuous variables
   contvar <- 1:sum(repeat_1)
   contvar_1 <- z_1[repeat_1]
   contvar_2 <- z_2[repeat_2]
   contvar[is.na(contvar_1) | is.na(contvar_2)] <- NA

# Assign values for survey design variables using the survey one sites
   if(revisitwgt) {
      wgt <- wgt_1[repeat_1]
   } else {
      wgt <- rep(1, length(contvar_1))
   }
   x <- x_1[repeat_1]
   y <- y_1[repeat_1] 
   stratum <- stratum_1[repeat_1]
   cluster <- cluster_1[repeat_1]
   wgt1 <- wgt1_1[repeat_1]
   x1 <- x1_1[repeat_1]
   y1 <- y1_1[repeat_1]
   popsize <- popsize_1
   pcfsize <- pcfsize_1
   N.cluster <- N.cluster_1
   stage1size <- stage1size_1
   support <- support_1[repeat_1]

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
   pcfactor.ind <- popcorrect_1

# Assign the value of sizeweight to the indicator variable for use of size
# weights
   swgt.ind <- sizeweight_1

# Assign the value of vartype to the character variable for type of variance
# estimator
   vartype <- vartype_1

# Remove missing values
   if(vartype == "Local") {
      if(swgt.ind) {
         if(stratum.ind) {
            if(cluster.ind)
               temp <- wnas(list(contvar=contvar, wgt=wgt, x=x, y=y,
                  stratum=stratum, cluster=cluster, wgt1=wgt1, x1=x1, y1=y1,
                  swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(contvar=contvar, wgt=wgt, x=x, y=y,
                  stratum=stratum, swgt=swgt))
         } else {
            if(cluster.ind)
               temp <- wnas(list(contvar=contvar, wgt=wgt, x=x, y=y,
                  cluster=cluster, wgt1=wgt1, x1=x1, y1=y1, swgt=swgt,
                  swgt1=swgt1))
            else
               temp <- wnas(list(contvar=contvar, wgt=wgt, x=x, y=y, swgt=swgt))
         }
         contvar <- temp$contvar
         contvar_1 <- contvar_1[contvar]
         contvar_2 <- contvar_2[contvar]
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
               temp <- wnas(list(contvar=contvar, wgt=wgt, x=x, y=y,
                  stratum=stratum, cluster=cluster, wgt1=wgt1, x1=x1, y1=y1))
            else
               temp <- wnas(list(contvar=contvar, wgt=wgt, x=x, y=y,
                  stratum=stratum))
         } else {
            if(cluster.ind)
               temp <- wnas(list(contvar=contvar, wgt=wgt, x=x, y=y,
                  cluster=cluster, wgt1=wgt1, x1=x1, y1=y1))
            else
               temp <- wnas(list(contvar=contvar, wgt=wgt, x=x, y=y))
         }
         contvar <- temp$contvar
         contvar_1 <- contvar_1[contvar]
         contvar_2 <- contvar_2[contvar]
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
               temp <- wnas(list(contvar=contvar, wgt=wgt, stratum=stratum,
                  cluster=cluster, wgt1=wgt1, swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(contvar=contvar, wgt=wgt, stratum=stratum,
                  swgt=swgt))
         } else {
            if(cluster.ind)
               temp <- wnas(list(contvar=contvar, wgt=wgt, cluster=cluster,
                  wgt1=wgt1, swgt=swgt, swgt1=swgt1))
            else
               temp <- wnas(list(contvar=contvar, wgt=wgt, swgt=swgt))
         }
         contvar <- temp$contvar
         contvar_1 <- contvar_1[contvar]
         contvar_2 <- contvar_2[contvar]
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
               temp <- wnas(list(contvar=contvar, wgt=wgt, stratum=stratum,
                  cluster=cluster, wgt1=wgt1))
            else
               temp <- wnas(list(contvar=contvar, wgt=wgt, stratum=stratum))
         } else {
            if(cluster.ind)
               temp <- wnas(list(contvar=contvar, wgt=wgt, cluster=cluster,
                  wgt1=wgt1))
            else
               temp <- wnas(list(contvar=contvar, wgt=wgt))
         }
         contvar <- temp$contvar
         contvar_1 <- contvar_1[contvar]
         contvar_2 <- contvar_2[contvar]
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
      if(nstrata > 0 & nstrata < nstrata.old) {
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
            contvar <- contvar[!stratum.i]
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
         warn.df <- rbind(warn.df, data.frame(func=I(fname), subpoptype=warn.vec[1],
            subpop=warn.vec[2], indicator=warn.vec[3], stratum=NA, warning=I(warn),
            action=I(act)))
      }
   }

# Calculate additional required values
   if(length(contvar_1) > 1) {
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
   }

# Branch to handle stratified and unstratified data

   if(stratum.ind) {

# Begin the section for stratified data

# Create the object for covariance or correlation estimates for all strata
# combined
   rslt <- NA

# Check whether the vectors of continuous variable values for revisit sites
# are empty or contain a single value
   if(length(contvar_1) <= 1) {
      warn.ind <- TRUE
      act <- "Covariance among the revisited sites was not included in calculation of \nthe standard error estimate.\n"
      warn <- paste("The number of nonmissing repeat visit sites was less than two in one of the \nsurveys.\n", sep="")
      warn.df <- rbind(warn.df, data.frame(func=I(fname),
         subpoptype=warn.vec[1], subpop=warn.vec[2], indicator=warn.vec[3],
         stratum=NA, warning=I(warn), action=I(act)))

# Begin section for nonempty vectors of continuous variable values for revisit
# sites
   } else {

# Begin the loop for individual strata

   for(i in 1:nstrata) {

# Check whether the vectors of continuous variable values for revisit sites
# are empty or contain a single value for a stratum
      stratum.i <- stratum == stratum.levels[i]
      if(length(contvar_1[stratum.i]) <= 1) {
         warn.ind <- TRUE
         act <- "Due to insufficient number of sites, the stratum was not included in \ncalculation of covariance among the revisited sites.\n"
         warn <- paste("The number of nonmissing repeat visit sites  in one of the surveys was less \nthan two for stratum \"", stratum.levels[i], "\".\n", sep="")
         warn.df <- rbind(warn.df, data.frame(func=I(fname),
            subpoptype=warn.vec[1], subpop=warn.vec[2],
            indicator=warn.vec[3], stratum=I(stratum.levels[i]),
            warning=I(warn), action=I(act)))

# Begin section for nonempty vectors of continuous variables values for revisit
# sites for a stratum
      } else {
# Calculate mean estimates
         stratum.i <- stratum == stratum.levels[i]
         z1 <- contvar_1[stratum.i]
         z2 <- contvar_2[stratum.i]
         if(swgt.ind) {
            if(cluster.ind) {
               w2 <- wgt[stratum.i]*swgt[stratum.i]
               w1 <- wgt1[stratum.i]*swgt1[stratum.i]
               mean1 <- sum(w2*w1, z1) / popsize.hat[i]
               mean2 <- sum(w2*w1, z2) / popsize.hat[i]
            } else {
               w <- wgt[stratum.i]*swgt[stratum.i]
               mean1 <- sum(w, z1) / popsize.hat[i]
               mean2 <- sum(w, z2) / popsize.hat[i]
            }
         } else {
            if(cluster.ind) {
               w2 <- wgt[stratum.i]
               w1 <- wgt1[stratum.i]
               mean1 <- sum(w2*w1, z1) / popsize.hat[i]
               mean2 <- sum(w2*w1, z2) / popsize.hat[i]
            } else {
               w <- wgt[stratum.i]
               mean1 <- sum(w, z1) / popsize.hat[i]
               mean2 <- sum(w, z2) / popsize.hat[i]
            }
         }

# Calculate covariance or correlation estimates
         if(cluster.ind) {
            temp <- changevar.mean(z1, z2, w2, x[stratum.i], y[stratum.i],
               revisitwgt, mean1, mean2, stratum.ind, stratum.levels[i],
               cluster.ind, cluster[stratum.i], w1, x1[stratum.i],
               y1[stratum.i], pcfactor.ind, NULL, N.cluster[i], stage1size[[i]],
               support[stratum.i], vartype, warn.ind, warn.df, warn.vec)
         } else {
            temp <- changevar.mean(z1, z2, w, x[stratum.i], y[stratum.i],
               revisitwgt, mean1, mean2, stratum.ind, stratum.levels[i],
               cluster.ind, pcfactor.ind=pcfactor.ind, pcfsize=pcfsize[i],
               support=support[stratum.i], vartype=vartype, warn.ind=warn.ind,
               warn.df=warn.df, warn.vec=warn.vec)
         }
         correst <- temp$rslt
         warn.ind <- temp$warn.ind
         warn.df <- temp$warn.df

# Add estimates to the object for all strata combined
         if(!is.null(popsize)) {
            rslt[!is.na(correst)] <- rslt[!is.na(correst)] +
               (popsize[i]/sum.popsize)*correst[!is.na(correst)]
         } else {
            rslt[!is.na(correst)] <- rslt[!is.na(correst)] +
               (popsize.hat[i]/sum.popsize.hat)*correst[!is.na(correst)]
         }

# End the section for nonempty vectors of categorical variables values for
# revisit sites for a stratum
      }

# End the loop for individual strata

   }

# End the section for nonempty vectors of categorical variables values for
# revisit sites
   }

# End the section for stratified data

   } else {

# Begin the section for unstratified data

# Check whether the vectors of categorical variable values for revisit sites
# are empty or contain a single value
   if(length(contvar_1) <= 1) {
      rslt <- NA
      warn.ind <- TRUE
      act <- "Covariance among the revisited sites was not included in calculation of \nthe standard error estimate.\n"
      warn <- paste("The number of nonmissing repeat visit sites was less than two in one of the \nsurveys.\n", sep="")
      warn.df <- rbind(warn.df, data.frame(func=I(fname),
         subpoptype=warn.vec[1], subpop=warn.vec[2], indicator=warn.vec[3],
         stratum=NA, warning=I(warn), action=I(act)))

# Begin section for nonempty vectors of categorical variables values for revisit
# sites
   } else {

# Calculate mean estimates
      z1 <- contvar_1
      z2 <- contvar_2
      if(swgt.ind) {
         if(cluster.ind) {
            w2 <- wgt*swgt
            w1 <- wgt1*swgt1
            mean1 <- sum(z1*w2*w1) / popsize.hat
            mean2 <- sum(z2*w2*w1) / popsize.hat
         } else {
            w <- wgt*swgt
            mean1 <- sum(z1*w) / popsize.hat
            mean2 <- sum(z2*w) / popsize.hat
         }
      } else {
         if(cluster.ind) {
            w2 <- wgt
            w1 <- wgt1
            mean1 <- sum(z1*w2*w1) / popsize.hat
            mean2 <- sum(z2*w2*w1) / popsize.hat
         } else {
            w <- wgt
            mean1 <- sum(z1*w) / popsize.hat
            mean2 <- sum(z2*w) / popsize.hat
         }
      }

# Calculate covariance or correlation estimates
      if(cluster.ind) {
         temp <- changevar.mean(z1, z2, w2, x, y, revisitwgt, mean1, mean2,
            stratum.ind, NULL, cluster.ind, cluster, w1, x1, y1, pcfactor.ind,
            NULL, N.cluster, stage1size, support, vartype, warn.ind, warn.df,
            warn.vec)
      } else {
         temp <- changevar.mean(z1, z2, w, x, y, revisitwgt, mean1, mean2,
            stratum.ind, NULL, cluster.ind, pcfactor.ind=pcfactor.ind,
            pcfsize=pcfsize, support=support, vartype=vartype,
            warn.ind=warn.ind, warn.df=warn.df, warn.vec=warn.vec)
      }
      rslt <- temp$rslt
      warn.ind <- temp$warn.ind
      warn.df <- temp$warn.df

# End the section for nonempty vectors of categorical variables values for
# revisit sites
   }

# End the section for unstratified data
   }

# Calculate standard errors
   if(is.na(rslt)) {
      Results$StdError <- sqrt(Results$StdError_1^2 + Results$StdError_2^2)
   } else {
      tw_1r <- sum(wgt_1[repeat_1])
      tw_2r <- sum(wgt_2[repeat_2])
      if(revisitwgt) {
         temp <- Results$StdError_1^2 + Results$StdError_2^2 -
            ((2*tw_1r*tw_2r)/(tw_1*tw_2))*rslt
         if(temp <= 0) {
            Results$StdError <- sqrt(Results$StdError_1^2 +
               Results$StdError_2^2)
         } else {
            Results$StdError <- sqrt(temp)
         }
      } else {
         temp <- total.est(z=z_1[repeat_1], wgt=wgt_1[repeat_1],
            x=x_1[repeat_1], y=y_1[repeat_1], stratum=stratum_1[repeat_1],
            cluster=cluster_1[repeat_1], wgt1=wgt1_1[repeat_1],
            x1=x1_1[repeat_1], y1=y1_1[repeat_1], popsize=popsize_1,
            popcorrect=popcorrect_1, pcfsize=pcfsize_1, N.cluster=N.cluster_1,
            stage1size=stage1size_1, support=support_1[repeat_1],
            sizeweight=sizeweight_1, swgt=swgt_1[repeat_1],
            swgt1=swgt1_1[repeat_1], vartype=vartype_1, conf=conf,
            check.ind=check.ind, warn.ind=warn.ind, warn.df=warn.df,
            warn.vec=warn.vec)
         se_1 <- temp$Results$StdError[temp$Results$Statistic == "Mean"]
         temp <- total.est(z=z_2[repeat_2], wgt=wgt_2[repeat_2],
            x=x_2[repeat_2], y=y_2[repeat_2], stratum=stratum_2[repeat_2],
            cluster=cluster_2[repeat_2], wgt1=wgt1_2[repeat_2],
            x1=x1_2[repeat_2], y1=y1_2[repeat_2], popsize=popsize_2,
            popcorrect=popcorrect_2, pcfsize=pcfsize_2, N.cluster=N.cluster_2,
            stage1size=stage1size_2, support=support_2[repeat_2],
            sizeweight=sizeweight_2, swgt=swgt_2[repeat_2],
            swgt1=swgt1_2[repeat_2], vartype=vartype_2, conf=conf,
            check.ind=check.ind, warn.ind=warn.ind, warn.df=warn.df,
            warn.vec=warn.vec)
         se_2 <- temp$Results$StdError[temp$Results$Statistic == "Mean"]
         covest <- rslt * se_1 * se_2
         temp <- Results$StdError_1^2 + Results$StdError_2^2 -
            ((2*tw_1r*tw_2r)/(tw_1*tw_2))*covest
         if(temp <= 0) {
            Results$StdError <- sqrt(Results$StdError_1^2 +
               Results$StdError_2^2)
         } else {
            Results$StdError <- sqrt(temp)
         }
      }
   }

# Calculate confidence bounds
   eval(parse(text=paste("Results$LCB", conf, "Pct <- Results$DiffEst - mult*Results$StdError", sep="")))
   eval(parse(text=paste("Results$UCB", conf, "Pct <- Results$DiffEst + mult*Results$StdError", sep="")))

# End the section for surveys with repeat visit sites
   }

# Rearrange columns of the Results data frame
   Results <- Results[,c(1, 12:15, 2:11)]

#
# End the section for a continuous variable
#

   } else {

# Print an error message for an unrecognized type of response variable
      stop(paste("\nThe value provided for argument resp.ind, ", resp.ind, ", is not a valid value", sep=""))
   }

# Depending on whether the function was called directly or was called by
# change.analysis, return appropriate results
   if(warn.vec[1] == 0) {

# As necessary, output a message indicating that warning messages were generated
# during execution of the program
      if(warn.ind) {
         warn.df$subpoptype[warn.df$subpoptype == 0] <- NA
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
