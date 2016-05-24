#' Density and abundance estimates and variances 
#'
#' Computes density and abundance estimates and variances based on
#' Horvitz-Thompson-like estimator
#'
#' Density and abundance within the sampled region is computed based on a
#' Horvitz-Thomspon-like estimator for groups and individuals (if a clustered
#' population) and this is extrapolated to the entire survey region based on
#' any defined regional stratification.  The variance is based on replicate
#' samples within any regional stratification.  For clustered populations, E(s)
#' and its standard error are also output.
#'
#' Abundance is estimated with a Horvitz-Thompson-like estimator (Huggins
#' 1989,1991; Borchers et al 1998; Borchers and Burnham 2004).  The abundance
#' in the sampled region is simply 1/p_1 + 1/p_2 + ... + 1/p_n where p_i is the
#' estimated detection probability for the ith detection of n total
#' observations. It is not strictly a Horvitz-Thompson estimator because the
#' p_i are estimated and not known. For animals observed in tight clusters,
#' that estimator gives the abundance of groups (\code{group=TRUE} in
#' \code{options}) and the abundance of individuals is estimated as s_1/p_1 +
#' s_2/p_2 + ... + s_n/p_n, where s_i is the size (e.g., number of animals in
#' the group) of each observation(\code{group=FALSE} in \code{options}).
#'
#' Extrapolation and estimation of abundance to the entire survey region is
#' based on either a random sampling design or a stratified random sampling
#' design.  Replicate samples(lines)(\code{sample.table} are specified within
#' regional strata \code{region.table}, if any.  If there is no stratification,
#' \code{region.table} should contain only a single record with the \code{Area}
#' for the entire survey region.  The \code{sample.table} is linked to the
#' \code{region.table} with the \code{Region.Label}.  The \code{obs.table} is
#' linked to the \code{sample.table} with the \code{Sample.Label} and
#' \code{Region.Label}.  Abundance can be restricted to a subset (e.g., for a
#' particular species) of the population by limiting the list the observations
#' in \code{obs.table} to those in the desired subset. Alternatively, if
#' \code{Sample.Label} and \code{Region.Label} are in the dataframe used to fit
#' the model, then a \code{subset} argument can be given in place of the
#' \code{obs.table}. To use the \code{subset} argument but include all of the
#' observations, use \code{subset=1==1} to avoid creating an \code{obs.table}.
#'
#' In extrapolating to the entire survey region it is important that the unit
#' measurements be consistent or converted for consistency.  A conversion
#' factor can be specified with the \code{convert.units} variable in the
#' \code{options} list.  The values of \code{Area} in \code{region.table}, must
#' be made consistent with the units for \code{Effort} in \code{sample.table}
#' and the units of \code{distance} in the dataframe that was analyzed.  It is
#' easiest to do if the units of \code{Area} is the square of the units of
#' \code{Effort} and then it is only necessary to convert the units of
#' \code{distance} to the units of \code{Effort}. For example, if \code{Effort}
#' was entered in kilometers and \code{Area} in square kilometers and
#' \code{distance} in meters then using
#' \code{options=list(convert.units=0.001)} would convert meters to kilometers,
#' density would be expressed in square kilometers which would then be
#' consistent with units for \code{Area}.  However, they can all be in
#' different units as long as the appropriate composite value for
#' \code{convert.units} is chosen.  Abundance for a survey region can be
#' expressed as: \code{A*N/a} where \code{A} is \code{Area} for the survey
#' region, \code{N} is the abundance in the covered (sampled) region, and
#' \code{a} is the area of the sampled region and is in units of \code{Effort *
#' distance}.  The sampled region \code{a} is multiplied by
#' \code{convert.units}, so it should be chosen such that the result is in the
#' same units of \code{Area}.  For example, if \code{Effort} was entered in
#' kilometers, \code{Area} in hectares (100m x 100m) and \code{distance} in
#' meters, then using \code{options=list(convert.units=10)} will convert
#' \code{a} to units of hectares (100 to convert meters to 100 meters for
#' distance and .1 to convert km to 100m units).
#'
#' If the argument \code{se} is set to \code{TRUE}, a standard error for
#' density and abundance is computed and the coefficient of variation and
#' log-normal confidence intervals are constructed using a Satterthwaite
#' approximation for degrees of freedom (Buckland et al. 2001 pg 90).  The
#' function \code{\link{dht.se}} computes the variance and interval estimates.
#' The variance has two components: 1) variation due to uncertanity from
#' estimation of the detection function and 2) variation in abundance due to
#' random sample selection.  The first component is computed using a delta
#' method estimate of variance (\code{\link{DeltaMethod}} (Huggins 1989, 1991,
#' Borchers et al. 1998) in which the first derivatives of the abundance
#' estimator with respect to the parameters in the detection function are
#' computed numerically.  The second component can be computed in one of three
#' ways as set by the option \code{varflag} with values 0,1,2. A value of 0 is
#' to use a binomial variance for the number of observations and it is only
#' useful if the sampled region is the survey region and the objects are not
#' clustered which will not occur very often.  A value of 1 uses the standard
#' variance for the encounter rate (Buckland et al. 2001 pg 78-79, although the
#' actual encounter rate formula used by default is now estimator R2 from
#' Fewster et al. (2009) - see \link{varn} for details).  If the population is
#' clustered the mean group size and standard error is also included.  This
#' variance estimator is not appropriate if \code{size} or a derivative of
#' \code{size} is used in the any of the detection function models. In general
#' if any covariates are used in the models, the default option 2 is
#' preferable.  It uses the variance estimator suggested by Innes et al (2002)
#' which used the formula for the variance ecounter rate but replaces the
#' number of observations per sample with the estimated abundance per sample.
#' This latter variance is also given in Marques and Buckland (2004).
#'
#' The argument \code{options} is a list of variable=value pairs that set
#' options for the analysis. All but one of these has been described so far.
#' The remaining variable \code{pdelta} should not need to be changed but was
#' included for completeness.  It controls the precision of the first
#' derivative calculation for the delta method variance.
#'
#' @param model ddf model object
#' @param region.table \code{data.frame} of region records. Two columns:
#'  \code{Region.Label} and \code{Area}.
#' @param sample.table \code{data.frame} of sample records. Three columns:
#'   \code{Region.Label}, \code{Sample.Label}, \code{Effort}.
#' @param obs.table \code{data.frame} of observation records with fields:
#'  \code{object}, \code{Region.Label}, and \code{Sample.Label} which give
#'  links to \code{sample.table}, \code{region.table} and the data records used
#'  in \code{model}. Not necessary if the \code{data.frame} used to create
#'  the model contains \code{Region.Label}, \code{Sample.Label} columns.
#' @param subset subset statement to create \code{obs.table}
#' @param se if \code{TRUE} computes std errors, cv and confidence interval
#'  based on log-normal
#' @param bootstrap if \code{TRUE} uses bootstrap approach (currently not
#'  implemented)
#' @param options a list of options that can be set, see "\code{dht} options", beow.
#' @export
#' @return list object of class \code{dht} with elements:
#' \item{clusters}{result list for object clusters}
#' \item{individuals}{result list for individuals}
#' \item{Expected.S}{\code{data.frame} of estimates of expected cluster size
#'  with fields \code{Region}, \code{Expected.S} and \code{se.Expected.S}
#'  If each cluster \code{size=1}, then the result only includes individuals
#'  and not clusters and \code{Expected.S}.}
#'
#' The list structure of clusters and individuals are the same:
#' \item{bysample}{\code{data.frame} giving results for each sample; Nchat is the
#'  estimated abundance within the sample and Nhat is scaled by surveyed area/
#'  covered area within that region}
#' \item{summary}{\code{data.frame} of summary statistics for each region and
#'  total}
#' \item{N}{\code{data.frame} of estimates of abundance for each region and
#'  total}
#' \item{D}{\code{data.frame} of estimates of density for each region and total}
#' \item{average.p}{average detection probability estimate}
#' \item{cormat}{correlation matrix of regional abundance/density estimates and
#'  total (if more than one region)}
#' \item{vc}{list of 3: total v-c matrix and detection and er (encounter rate)
#'  components of variance; for detection the v-c matrix and partial vector
#'  are returned}
#' \item{Nhat.by.sample}{another summary of \code{Nhat} by sample used by
#'  \code{dht.se}}
#'
#'
#' @section \code{dht} options:
#'  Several options are available to control calculations and output:
#'
#' \describe{
#'  \item{\code{ci.width}}{Confidence iterval width, expressed as a decimal between 0 and 1 (default 0.95, giving a 95\% CI)}
#'  \item{\code{pdelta}}{ delta value for computing numerical first derivatives (Default: 0.001)}
#'  \item{\code{varflag}}{ 0,1,2 (see Details) (Default: 2)}
#'  \item{\code{convert.units}}{ multiplier for width to convert to units of length (Default: 1)}
#'  \item{\code{ervar}}{ encounter rate variance type - see type argument to \code{\link{varn}} (Default: "R2")}
#'}
#'
#' @author Jeff Laake
#' @seealso print.dht dht.se
#' @references
#'
#' Borchers, D.L., S.T. Buckland, P.W. Goedhart, E.D. Clarke, and S.L. Hedley.
#'   1998.  Horvitz-Thompson estimators for double-platform line transect
#'   surveys.  Biometrics 54: 1221-1237.
#'
#' Borchers, D.L. and K.P. Burnham. General formulation for distance sampling
#'   pp 10-11 In: Advanced Distance Sampling, eds. S.T. Buckland, D.R.Anderson,
#'   K.P. Burnham, J.L. Laake, D.L. Borchers, and L. Thomas. Oxford University
#'   Press.
#'
#' Buckland, S.T., D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and
#'   L. Thomas. 2001.  Introduction to Distance Sampling: Estimating Abundance
#'   of Biological Populations. Oxford University Press.
#'
#' Fewster, R.M., S.T. Buckland, K.P. Burnham, D.L. Borchers, P.E. Jupp, J.L.
#'   Laake and L. Thomas. 2009.  Estimating the encounter rate variance in
#'   distance sampling. Biometrics 65: 225-236.
#'
#' Huggins, R.M. 1989.  On the statistical analysis of capture experiments.
#'   Biometrika 76:133-140.
#'
#' Huggins, R.M. 1991. Some practical aspects of a conditional likelihood
#'   approach to capture experiments.  Biometrics 47: 725-732.
#'
#' Innes, S.  M.P. Heide-Jorgensen, J.L. Laake, K.L. Laidre, H.J. Cleator, P.
#'   Richard, and R.E.A. Stewart.  2002. Surveys of belugas and narwhals in the
#'   Canadian High Arctic in 1996.  NAMMCO Scientific Publications 4: 169-190.
#'
#' Marques, F.F.C. and S.T. Buckland. 2004. Covariate models for the detection
#'   function. In: Advanced Distance Sampling, eds. S.T. Buckland,
#'   D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L. Thomas.
#'   Oxford University Press.
#' @keywords utility
dht <- function(model,region.table,sample.table, obs.table=NULL, subset=NULL,
                se=TRUE, bootstrap=FALSE, options=list()){
  # Functions Used:  assign.default.values, create.varstructure,
  #                  covered.region.dht, survey.region.dht, dht.se, varn,
  #                  covn(in varn.R), solvecov (in coef.ds.R).

  tables.dht <- function(group){
    # Internal function to create summary tables for clusters (group=TRUE) or
    # individuals (group=FALSE).
    options$group <- group

    # Compute covered region abundances by sample depending on value of group
    Nhat.by.sample <- covered.region.dht(obs, samples, group)

    # Mod 18-Aug-05 jll; added computation of avergage detection probability
    # which is simply n/Nhat in the covered region
    average.p <- nrow(obs)/sum(Nhat.by.sample$Nhat)

    # Scale up abundances to survey region
    # jll 19-Jan-05 - sort Nhat.by.sample by Region.Label and Sample.Label
    width <- model$meta.data$width * options$convert.units
    Nhat.by.sample <- survey.region.dht(Nhat.by.sample, samples,width,point)
    Nhat.by.sample <- Nhat.by.sample[order(Nhat.by.sample$Region.Label,
                                           Nhat.by.sample$Sample.Label),]
    if(point){
      s.area <- Nhat.by.sample$Effort.x*pi*width^2
    }else{
      s.area <- Nhat.by.sample$Effort.x*2*width
    }

    bysample.table <- with(Nhat.by.sample,
                           data.frame(Region      = Region.Label,
                                      Sample      = Sample.Label,
                                      Effort      = Effort.x,
                                      Sample.Area = s.area,
                                      Area        = Area,
                                      n           = n,
                                      Nhat        = Nhat,
                                      Nchat       = Nhat*CoveredArea/Area))

    bysample.table$Dhat <- bysample.table$Nchat/bysample.table$Sample.Area
    Nhat.by.region <- by(Nhat.by.sample$Nhat, Nhat.by.sample$Region.Label,sum)

    # Create estimate table
    numRegions <- length(unique(samples$Region.Label))
    if(numRegions > 1){
      estimate.table <- data.frame(
                          Label = c(levels(unique(samples$Region.Label)),"Total"),
                          Estimate = rep(0, numRegions + 1),
                          se = rep(NA,numRegions + 1),
                          cv = rep(NA, numRegions + 1),
                          lcl = rep(NA,numRegions + 1),
                          ucl = rep(NA, numRegions + 1))
    }else{
      estimate.table = data.frame(Label    = c("Total"),
                                  Estimate = rep(0,1),
                                  se       = rep(NA, 1),
                                  cv       = rep(NA, 1),
                                  lcl      = rep(NA, 1),
                                  ucl      = rep(NA, 1))
    }

    if(numRegions > 1){
      estimate.table$Estimate <- c(Nhat.by.region, sum(Nhat.by.region))
    }else{
      estimate.table$Estimate <- Nhat.by.region
    }

    # Create summary table
    summary.table <- Nhat.by.sample[, c("Region.Label","Area",
                                        "CoveredArea","Effort.y")]
    summary.table <- unique(summary.table)
    var.er <- sapply(split(Nhat.by.sample,Nhat.by.sample$Region.Label),
                     function(x)varn(x$Effort.x,x$n,type=options$ervar))

    if(numRegions > 1){
       var.er <- c(var.er,varn(Nhat.by.sample$Effort.x,
                               Nhat.by.sample$n,type=options$ervar))
    }

    #  jll 11_11_04; change to set missing values for nobs to 0
    #   - regions with no sightings
    nobs <- as.vector(by(bysample.table$n, bysample.table$Region, sum))
    nobs[is.na(nobs)] <- 0
    summary.table$n <- nobs
    if(group){
      summary.table$k <- tapply(Nhat.by.sample$Sample.Label,
                                Nhat.by.sample$Region.Label,length)
      colnames(summary.table) <- c("Region", "Area", "CoveredArea",
                                   "Effort", "n","k")
    }else{
      colnames(summary.table) = c("Region", "Area", "CoveredArea", "Effort", "n")
    }

    if(numRegions > 1){
      summary.table <- data.frame(Region=c(levels(summary.table$Region),"Total"),
                                  rbind(summary.table[, -1],
                                        apply(summary.table[,-1], 2, sum)))
    }

    summary.table$ER <- summary.table$n/summary.table$Effort
    summary.table$se.ER <- sqrt(var.er)
    summary.table$cv.ER <- summary.table$se.ER/summary.table$ER
    summary.table$cv.ER[summary.table$ER==0] <- 0

    # set missing values to 0
    summary.table$ER[is.nan(summary.table$ER)] <- 0
    summary.table$se.ER[is.nan(summary.table$se.ER)] <- 0
    summary.table$cv.ER[is.nan(summary.table$cv.ER)] <- 0

    # If summary of individuals for a clustered popn, add mean
    #  group size and its std error
    if(!group){
      mean.clustersize <- tapply(obs$size,obs$Region.Label,mean)
      se.clustersize <- sqrt(tapply(obs$size,obs$Region.Label,var)/
                             tapply(obs$size,obs$Region.Label,length))
      cs <- data.frame(Region=names(mean.clustersize),
                       mean.size=as.vector(mean.clustersize),
                       se.mean=as.vector(se.clustersize))

      summary.table <- merge(summary.table, cs, by.x = "Region",
                             all=TRUE,sort=FALSE)
      if(numRegions > 1){
        summary.table$mean.size[numRegions+1] <- mean(obs$size)
        summary.table$se.mean[numRegions+1] <-sqrt(var(obs$size)/length(obs$size))
      }
      #  29/05/12 lhm - moved to set missing values to 0
      summary.table$mean.size[is.na(summary.table$mean.size)] <- 0
      summary.table$se.mean[is.na(summary.table$se.mean)] <- 0
    }

    rownames(summary.table) <- 1:dim(summary.table)[1]

    # If a std error has been requested call dht.se
    if(se){
      result <- dht.se(model, summary.table, samples, obs, options, numRegions,
                       estimate.table, Nhat.by.sample)
      estimate.table <- result$estimate.table
    }

    # Create estimate table for D from same table for N
    D.estimate.table <- estimate.table
    if(numRegions > 1){
      D.estimate.table$Estimate <-D.estimate.table$Estimate/c(region.table$Area,
          sum(region.table$Area))
      D.estimate.table$se <- D.estimate.table$se/c(region.table$Area,
          sum(region.table$Area))
      D.estimate.table$lcl <- D.estimate.table$lcl/c(region.table$Area,
          sum(region.table$Area))
      D.estimate.table$ucl <- D.estimate.table$ucl/c(region.table$Area,
          sum(region.table$Area))
    }else{
      D.estimate.table$Estimate <- D.estimate.table$Estimate/region.table$Area
      D.estimate.table$se <- D.estimate.table$se/region.table$Area
      D.estimate.table$lcl <- D.estimate.table$lcl/region.table$Area
      D.estimate.table$ucl <- D.estimate.table$ucl/region.table$Area
    }

    # set missing values to 0
    D.estimate.table$Estimate[is.nan(D.estimate.table$Estimate)] <- 0
    D.estimate.table$se[is.nan(D.estimate.table$se)] <- 0
    D.estimate.table$cv[is.nan(D.estimate.table$cv)] <- 0
    D.estimate.table$lcl[is.nan(D.estimate.table$lcl)] <- 0
    D.estimate.table$ucl[is.nan(D.estimate.table$ucl)] <- 0

    # Return list depending on value of se
    # change to set missing values to 0
    # jll 6/30/06; dropped restriction that numregions > 1 on sending vc back
    if(se){
      cormat <- result$vc/(result$estimate.table$se %o% result$estimate.table$se)
      cormat[is.nan(cormat)] <- 0
      result <- list(bysample=bysample.table, summary = summary.table,
                     N=result$estimate.table, D=D.estimate.table, 
                     average.p=average.p, cormat = cormat,
                     vc=list(total=result$vc,detection=result$vc1,er=result$vc2),
                     Nhat.by.sample=Nhat.by.sample)
    }else{
      result <- list(bysample=bysample.table,summary = summary.table,
                    N = estimate.table,D = D.estimate.table, average.p=average.p,
                    Nhat.by.sample=Nhat.by.sample)
    }
    return(result)
  }

###Start of dht function

  # Code additions by jll 18-Nov-04; the following allows for a subset
  # statement to be added to create obs.table from model data rather than
  # creating obs.table separately. This only works if the data contain the
  # Sample.Label and Region.Label fields.
  point <- model$meta.data$point
  objects <- as.numeric(names(model$fitted))
  if(is.null(obs.table)){
    data <- model$data
    if("observer"%in%names(data)){
      # jll 3 Sept 2014 if dual observer I added code to use observer 1 only
      # or it was doubling sample size
      data <- data[data$observer==1,]
    }
    if("Sample.Label" %in% names(data) & "Region.Label" %in% names(data)){
      if(is.null(substitute(subset))){
         obs.table <- data[,c("object","Sample.Label","Region.Label")]
      }else{
         select <- data[eval(substitute(subset),envir=data),]
         obs.table <- select[,c("object","Sample.Label","Region.Label")]
      }
      obs.table <- obs.table[obs.table$object %in% objects,]
    }else{
      stop("Must specify obs.table because Sample.Label and/or Region.Label fields not contained in data")
    }
  }

  # Extract relevant fields from Region and Sample tables; jll 4 May 07;
  region.table <- region.table[,c("Region.Label","Area")]
  sample.table <- sample.table[,c("Region.Label","Sample.Label","Effort")]

  # Make sure input data labels are factors
  region.table$Region.Label <- factor(region.table$Region.Label)
  sample.table$Region.Label <- factor(sample.table$Region.Label,
                                      levels=levels(region.table$Region.Label))
  obs.table$Region.Label <- factor(obs.table$Region.Label,
                                   levels=levels(region.table$Region.Label))
  sample.table$Sample.Label <- factor(sample.table$Sample.Label)
  obs.table$Sample.Label <- factor(obs.table$Sample.Label,
                                   levels=levels(sample.table$Sample.Label))

  # Assign default values to options
  options <- assign.default.values(options, pdelta = 0.001, varflag = 2,
                                   convert.units = 1, ervar="R2")

  # Convert width value
  width <- model$meta.data$width * options$convert.units

  # If area is zero for all regions reset to the area of the covered region
  DensityOnly <- FALSE
  if(sum(region.table$Area)==0){
    DensityOnly <- TRUE
    # cat("Warning: Area for regions is zero. They have been set to area of covered region(strips), \nso N is for covered region.",
    #     "However, standard errors will not match \nprevious covered region SE because it includes spatial variation\n")
    Effort.by.region <- by(sample.table$Effort, sample.table$Region.Label,sum)
    region.table$Area <- ifelse(point,
                                pi*as.vector(Effort.by.region)*width^2,
                                2*as.vector(Effort.by.region)*width)
  }

  # Create obs/samples structures
  vs <- create.varstructure(model, region.table, sample.table, obs.table)
  samples <- vs$samples
  obs <- vs$obs
  region.table <- vs$region

  # handle subset feature when labels are also in data
  if(!is.null(obs$Region.Label.x)){
    obs$Region.Label <- obs$Region.Label.x
    obs$Sample.Label <- obs$Sample.Label.x
    obs$Region.Label.x <- NULL
    obs$Sample.Label.x <- NULL
    obs$Region.Label.y <- NULL
    obs$Sample.Label.y <- NULL
  }

  # Merge with fitted values
  pdot <- model$fitted
  obs <- merge(obs,data.frame(object=objects,pdot=pdot))

  # If clustered population create tables for clusters and individuals and
  # an expected S table otherwise just tables for individuals in an
  # unclustered popn
  if(!is.null(obs$size)){
    clusters <- tables.dht(TRUE)
    individuals <- tables.dht(FALSE)
    Expected.S <- individuals$N$Estimate/clusters$N$Estimate

    # This computes the se(E(s)).  It essentially uses 3.37 from Ads but in
    # place of using 3.25, 3.34 and 3.38, it uses 3.27,3.35 and an equivalent
    # cov replacement term for 3.38.  This uses line to line variability
    # whereas the other formula measure the variance of E(s) within the lines 
    # and it goes to zero as p approaches 1.
    if(se & options$varflag!=1){
      if(options$varflag==2){
        numRegions <- length(unique(samples$Region.Label))
        cov.Nc.Ncs <- rep(0,numRegions)
        scale <- clusters$summary$Area/clusters$summary$CoveredArea

        for(i in 1:numRegions){
          c.stratum.data <- clusters$Nhat.by.sample[
              as.character(clusters$Nhat.by.sample$Region.Label) ==
              as.character(region.table$Region.Label[i]), ]

          i.stratum.data <- individuals$Nhat.by.sample[
              as.character(individuals$Nhat.by.sample$Region.Label) ==
              as.character(region.table$Region.Label[i]), ]

          Li <- sum(c.stratum.data$Effort.x)
          cov.Nc.Ncs[i] <- covn(c.stratum.data$Effort.x/(scale[i]*Li),
                                c.stratum.data$Nhat/scale[i],
                                i.stratum.data$Nhat/scale[i],
                                options$ervar)
        }
      }else{
        cov.Nc.Ncs <- as.vector(by(obs$size*(1 - obs$pdot)/obs$pdot^2,
                                   obs$Region.Label, sum))
        cov.Nc.Ncs[is.na(cov.Nc.Ncs)] <- 0
      }

      cov.Nc.Ncs[is.nan(cov.Nc.Ncs)] <- 0
      cov.Nc.Ncs <- c(cov.Nc.Ncs,sum(cov.Nc.Ncs))
      cov.Nc.Ncs <- cov.Nc.Ncs+diag(t(clusters$vc$detection$partial)%*%
                      solvecov(model$hessian)$inv%*%
                      individuals$vc$detection$partial)
      se.Expected.S <- clusters$N$cv^2 + individuals$N$cv^2 -
                  2*cov.Nc.Ncs/(individuals$N$Estimate*clusters$N$Estimate)
      Expected.S[is.nan(Expected.S)] <- 0
      se.Expected.S[se.Expected.S<=0 | is.nan(se.Expected.S)] <- 0
      se.Expected.S <- Expected.S*sqrt(se.Expected.S)
      Expected.S <- data.frame(Region=clusters$N$Label,
                               Expected.S=as.vector(Expected.S),
                               se.Expected.S=as.vector(se.Expected.S))
    }else{
      Expected.S <- data.frame(Region=clusters$N$Label,
                               Expected.S=as.vector(Expected.S))
    }

    if(DensityOnly){
      clusters$N <- NULL
      individuals$N <- NULL
    }

    result <- list(clusters=clusters,
                   individuals=individuals,
                   Expected.S=as.vector(Expected.S))
  }else{
    individuals <- tables.dht(TRUE)
    if(DensityOnly){
      individuals$N <- NULL
    }
    result <- list(individuals=individuals)
  }

  class(result) <- "dht"
  return(result)
}
