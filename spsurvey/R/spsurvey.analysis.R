spsurvey.analysis <- function(sites=NULL, subpop=NULL, design=NULL,
   data.cat=NULL, data.cont=NULL, siteID=NULL, wgt=NULL, sigma=NULL,
   var.sigma=NULL, xcoord=NULL, ycoord=NULL, stratum=NULL, cluster=NULL,
   wgt1=NULL, xcoord1=NULL, ycoord1=NULL, popsize=NULL, popcorrect=FALSE,
   pcfsize=NULL, N.cluster=NULL, stage1size=NULL, support=NULL,
   sizeweight=FALSE, swgt=NULL, swgt1=NULL, vartype="Local", conf=95,
   pctval=c(5,10,25,50,75,90,95)) {

################################################################################
# Function: spsurvey.analysis
# Purpose: Create an Object of Class spsurvey.analysis
# Programmer: Tom Kincaid
# Date: September 29, 2003
# Last Revised: August 19, 2014
# Description:
#   This function creates an object of class spsurvey.analysis that contains all 
#   of the information necessary to use the analysis functions in the 
#   spsurvey package.
# Arguments:
#   sites = a data frame consisting of two variables: the first variable is site
#     IDs and the second variable is a logical vector indicating which sites to
#     use in the analysis.  If this data frame is not provided, then the data
#     frame will be created, where (1) site IDs are obtained either from the
#     design argument, the siteID argument, or both (when siteID is a formula);
#     and (2) a variable named use.sites that contains the value TRUE for all
#     sites is created.  The default is NULL.
#   subpop = a data frame describing sets of populations and subpopulations for
#     which estimates will be calculated.  The first variable is siteIDs and
#     each subsequent variable identifies a Type of population, where the
#     variable name is used to identify Type.  A Type variable identifies each
#     site with one of the subpopulations of that Type.  If this data frame is
#     not provided, then the data frame will be created, where (1) site IDs are
#     obtained either from the design argument, the siteID argument, or both
#     (when siteID is a formula); and (2) a single Type variable named all.sites
#     that contains the value "All Sites" for all sites is created.  The default
#     is NULL.
#   design = a data frame consisting of design variables.  If variable names are
#     provided as formulas in the corresponding arguments, then the formulas are
#     interpreted using this data frame.  If this data frame is not provided,
#     then the data frame will be created from inputs to the design variables in
#     the argument list.  The default is NULL.  If variable names are not
#     provided as formulas, then variables should be named as follows:
#       siteID = site IDs
#       wgt = final adjusted weights
#       xcoord = x-coordinates for location
#       ycoord = y-coordinates for location
#       stratum = stratum codes
#       cluster = stage one sampling unit codes
#       wgt1 = final adjusted stage one weights
#       xcoord1 = stage one x-coordinates for location
#       ycoord1 = stage one y-coordinates for location
#       support = support values
#       swgt = size-weights
#       swgt1 = stage one size-weights
#   data.cat = a data frame of categorical response variables.  The first
#     variable is site IDs.  Subsequent variables are response variables.
#     Missing data (NA) is allowed.  If psurvey.obj is not provided, then this
#     argument is required.  The default is NULL.
#   data.cont = a data frame of continuous response variables.  The first
#     variable is site IDs.  Subsequent variables are response variables.
#     Missing data (NA) is allowed.  The default is NULL.
#   siteID = site IDs.  This variable can be input directly or as a formula and
#     must be supplied either as this argument or in the design data frame.  The
#     default is NULL.
#   wgt = final adjusted weights, which are either the weights for a single-
#     stage sample or the stage two weights for a two-stage sample.  This
#     variable can be input directly or as a formula and must be supplied either
#     as this argument or in the design data frame.  The default is NULL.
#   sigma = measurement error variance.  This variable must be a vector
#     containing a value for each response variable and must have the names
#     attribute set to identify the response variable names.  Missing data (NA)
#     is allowed.  The default is NULL.
#   var.sigma = variance of the measurement error variance.  This variable must
#     be a vector containing a value for each response variable and must have
#     the names attribute set to identify the response variable names.  Missing
#     data (NA) is allowed.  The default is NULL.
#   xcoord = x-coordinates for location, which are either the x-coordinates for
#     a single-stage sample or the stage two x-coordinates for a two-stage
#     sample.  This variable can be input directly or as a formula and must be
#     supplied either as this argument or in the design data frame when argument
#     vartype is set to "Local".  The default is NULL.
#   ycoord = y-coordinates for location, which are either the y-coordinates for
#     a single-stage sample or the stage two y-coordinates for a two-stage
#     sample.  This variable can be input directly or as a formula and must be
#     supplied either as this argument or in the design data frame when argument
#     vartype is set to "Local".  The default is NULL.
#   stratum = the stratum codes.  This variable can be input directly or as a
#     formula.  The default is NULL.
#   cluster = the stage one sampling unit (primary sampling unit or cluster)
#     codes.  This variable can be input directly or as a formula.  The default
#     is NULL.
#   wgt1 = final adjusted stage one weights.  This variable can be input
#     directly or as a formula.  The default is NULL.
#   xcoord1 = the stage one x-coordinates for location.  This variable can be
#     input directly or as a formula.  The default is NULL.
#   ycoord1 = the stage one y-coordinates for location.  This variable can be
#     input directly or as a formula.  The default is NULL.
#   popsize = known size of the resource, which is used to perform ratio
#     adjustment to estimators expressed using measurement units for the
#     resource.  For a finite resource, this argument is either the total number
#     of sampling units or the known sum of size-weights.  For an extensive
#     resource, this argument is the measure of the resource, i.e., either known
#     total length for a linear resource or known total area for an areal
#     resource.  The argument must be in the form of a list containing an
#     element for each population Type in the subpop data frame, where NULL is a
#     valid choice for a population Type.  The list must be named using the
#     column names for the population Types in subpop. If a population Type
#     doesn't contain subpopulations, then each element of the list is either a
#     single value for an unstratified sample or a vector containing a value for
#     each stratum for a stratified sample, where elements of the vector are
#     named using the stratum codes.  If a population Type contains
#     subpopulations, then each element of the list is a list containing an
#     element for each subpopulation, where the list is named using the
#     subpopulation names.  The element for each subpopulation will be either a
#     single value for an unstratified sample or a named vector of values for a
#     stratified sample.  The default is NULL.
#     Example popsize for a stratified sample:
#       popsize = list("Pop 1"=c("Stratum 1"=750,
#                                "Stratum 2"=500,
#                                "Stratum 3"=250),
#                      "Pop 2"=list("SubPop 1"=c("Stratum 1"=350,
#                                                "Stratum 2"=250,
#                                                "Stratum 3"=150),
#                                   "SubPop 2"=c("Stratum 1"=250,
#                                                "Stratum 2"=150,
#                                                "Stratum 3"=100),
#                                   "SubPop 3"=c("Stratum 1"=150,
#                                                "Stratum 2"=150,
#                                                "Stratum 3"=75)),
#                      "Pop 3"=NULL)
#     Example popsize for an unstratified sample:
#       popsize = list("Pop 1"=1500,
#                      "Pop 2"=list("SubPop 1"=750,
#                                   "SubPop 2"=500,
#                                   "SubPop 3"=375),
#                      "Pop 3"=NULL)
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
#     a site from anextensive resource, which is required for calculation of
#     finite and continuous population correction factors.  This variable can be
#     input directly or as a formula.  The default is NULL.
#   sizeweight = a logical value that indicates whether size-weights should be
#     used in the analysis, where TRUE = use the size-weights and FALSE = do not
#     use the size-weights.  The default is FALSE.
#   swgt = the size-weight for each site, which is the stage two size-weight for
#     a two-stage sample.  This variable can be input directly or as a formula.
#     The default is NULL.
#   swgt1 = the stage one size-weight for each site.  This variable can be input
#     directly or as a formula.  The default is NULL.
#   vartype = the choice of variance estimator, where "Local" = local mean
#     estimator and "SRS" = SRS estimator.  The default is "Local".
#   conf = the confidence level.  The default is 95%.
#   pctval = the set of values at which percentiles are estimated.  The default
#     set is: {5, 10, 25, 50, 75, 90, 95}.
# Output:
#   A list of class spsurvey.analysis.  Only those sites indicated by the logical
#   variable in the sites data frame are retained in the output. The sites,
#   subpop, and design data frames will always exist in the output. At least one
#   of the data.cat and data.cont data frames will exist.  Depending upon values
#   of the input variables, other elements in the output may be NULL.  The
#   output list is composed of the following elements:
#     sites = the sites data frame.
#     subpop = the subpop data frame.
#     design = the design data frame.
#     data.cat = the data.cat data frame.
#     data.cont = the data.cont data frame.
#     sigma = measurement error variance.
#     var.sigma = variance of the estimated measurement error variance.
#     stratum.ind = a logical value that indicates whether the sample is
#       stratified, where TRUE = a stratified sample and FALSE = not a
#       stratified sample.
#     cluster.ind = a logical value that indicates whether the sample is a two-
#       stage sample, where TRUE = a two-stage sample and FALSE = not a two-
#       stage sample.
#     popsize = the known size of the resource.
#     pcfactor.ind = a logical value that indicates whether the population
#       correction factor is used during variance estimation, where TRUE = use
#       the population correction factor and FALSE = do not use the factor.
#     pcfsize = size of the resource, which is required for calculation of
#       finite and continuous population correction factors for a single-stage
#       sample.
#     N.cluster = the number of stage one sampling units in the resource 
#     stage1size = size of the stage one sampling units of a two-stage sample.
#     swgt.ind = a logical value that indicates whether the sample is a size-
#       weighted sample, where TRUE = a size-weighted sample and FALSE = not a
#       size-weighted sample.
#     vartype = the choice of variance estimator, where "Local" = local mean
#       estimator and "SRS" = SRS estimator.
#     conf = the confidence level.
#     pctval = the set of values at which percentiles are estimated.
# Other Functions Required:
#   dframe.check - check site IDs, the sites data frame, the subpop data frame,
#     and the data.cat data frame to assure valid contents and, as necessary,
#     create the sites data frame and the subpop data frame
#   vecprint - takes an input vector and outputs a character string with line
#     breaks inserted
#   uniqueID - creates unique site IDs by appending a unique number to each
#     occurrence of a site ID
#   input.check - check input values for errors, consistency, and compatibility
#     with analytical functions
# Examples:
#   Categorical variable example:
#     mysiteID <- paste("Site", 1:100, sep="")
#     mysites <- data.frame(siteID=mysiteID, Active=rep(TRUE, 100))
#     mysubpop <- data.frame(siteID=mysiteID, All.Sites=rep("All Sites", 100),
#       Resource.Class=rep(c("Good","Poor"), c(55,45)))
#     mydesign <- data.frame(siteID=mysiteID, wgt=runif(100, 10, 100),
#       xcoord=runif(100), ycoord=runif(100), stratum=rep(c("Stratum1",
#       "Stratum2"), 50))
#     mydata.cat <- data.frame(siteID=mysiteID, CatVar=rep(c("north", "south",
#       "east", "west"), 25))
#     mypopsize <- list(All.Sites=c(Stratum1=3500, Stratum2=2000),
#       Resource.Class=list(Good=c(Stratum1=2500, Stratum2=1500),
#       Poor=c(Stratum1=1000, Stratum2=500)))
#     spsurvey.analysis(sites=mysites, subpop=mysubpop, design=mydesign,
#       data.cat=mydata.cat, popsize=mypopsize)
#
#   Continuous variable example - including deconvolution estimates:
#     mydesign <- data.frame(ID=mysiteID, wgt=runif(100, 10, 100),
#       xcoord=runif(100), ycoord=runif(100), stratum=rep(c("Stratum1",
#       "Stratum2"), 50))
#     ContVar <- rnorm(100, 10, 1)
#     mydata.cont <- data.frame(siteID=mysiteID, ContVar=ContVar,
#       ContVar.1=ContVar + rnorm(100, 0, sqrt(0.25)),
#       ContVar.2=ContVar + rnorm(100, 0, sqrt(0.50)))
#     mysigma <- c(ContVar=NA, ContVar.1=0.25, ContVar.2=0.50)
#     spsurvey.analysis(sites=mysites, subpop=mysubpop[,1:2], design=mydesign,
#       data.cont=mydata.cont, siteID=~ID, sigma=mysigma,
#       popsize=mypopsize[1])
################################################################################

# Create a data frame for warning messages

   warn.ind <- FALSE
   warn.df <- NULL
   fname <- "spsurvey.analysis"

# Create the vector of design variable names

   design.names <- c("siteID", "wgt", "xcoord", "ycoord", "stratum", "cluster",
      "wgt1", "xcoord1", "ycoord1", "support", "swgt", "swgt1")

# When the design argument is supplied, check that the input value is a data
# frame

   if(!is.null(design) && !is.data.frame(design))
      stop("\nThe design argument must be a data frame.")

# Obtain site IDs

   if(is.null(siteID)) {
      if(is.null(design))
         stop("\nSince neither argument design nor argument siteID were supplied, site IDs \ncannot be obtained.")
      if(is.null(design$siteID))
         stop("\nArgument siteID was not supplied, and argument design did not contain a \nvariable named siteID.")
      siteID <- design$siteID
      nresp <- length(siteID)
   } else {
      if(inherits(siteID, "formula")) {
         if(is.null(design))
            stop("\nArgument siteID was a formula, and argument design was not supplied.")
         temp <- match(as.character(siteID)[2], names(design), nomatch=0)
         if(temp == 0)
            stop("\nArgument siteID was a formula, and argument design did not contain a column \nmatching the name used in the formula.")
         siteID <- model.frame(siteID, design)
         siteID <- siteID[,1]
         nresp <- length(siteID)
      } else {
         if(!is.factor(siteID))
            siteID <- as.factor(siteID)
         nresp <- length(siteID)
      }
   }

# Obtain final adjusted weights

   if(is.null(wgt)) {
      if(is.null(design))
         stop("\nSince neither argument design nor argument wgt were supplied, final adjusted \nweights cannot be obtained.")
      if(is.null(design$wgt))
         stop("\nArgument wgt was not supplied, and argument design did not contain a variable \nnamed wgt.")
      wgt <- design$wgt
   } else {
      if(inherits(wgt, "formula")) {
         if(is.null(design))
            stop("\nArgument wgt was a formula, and argument design was not supplied.")
         temp <- match(as.character(wgt)[2], names(design), nomatch=0)
         if(temp == 0)
            stop("\nArgument wgt was a formula, and argument design did not contain a column \nmatching the name used in the formula.")
         wgt <- model.frame(wgt, design)
         wgt <- wgt[,1]
      }
   }
   if(!is.numeric(wgt))
      stop("\nThe final adjusted weights are not numeric values.")

# Obtain x-coordinates for location

   if(is.null(xcoord)) {
      if(!is.null(design$xcoord))
         xcoord <- design$xcoord
   } else {
      if(inherits(xcoord, "formula")) {
         if(is.null(design))
            stop("\nArgument xcoord was a formula, and argument design was not supplied.")
         temp <- match(as.character(xcoord)[2], names(design), nomatch=0)
         if(temp == 0)
            stop("\nArgument xcoord was a formula, and argument design did not contain a column \nmatching the name used in the formula.")
         xcoord <- model.frame(xcoord, design)
         xcoord <- xcoord[,1]
      }
   }
   if(!is.numeric(xcoord) && !is.null(xcoord))
      stop("\nThe x-coordinates for location are not numeric values.")

# Obtain y-coordinates for location

   if(is.null(ycoord)) {
      if(!is.null(design$ycoord))
         ycoord <- design$ycoord
   } else {
      if(inherits(ycoord, "formula")) {
         if(is.null(design))
            stop("\nArgument ycoord was a formula, and argument design was not supplied.")
         temp <- match(as.character(ycoord)[2], names(design), nomatch=0)
         if(temp == 0)
            stop("\nArgument ycoord was a formula, and argument design did not contain a column \nmatching the name used in the formula.")
         ycoord <- model.frame(ycoord, design)
         ycoord <- ycoord[,1]
      }
   }
   if(!is.numeric(ycoord) && !is.null(ycoord))
      stop("\nThe y-coordinates for location are not numeric values.")

# Obtain stratum codes

   if(is.null(stratum)) {
      if(!is.null(design$stratum))
         stratum <- design$stratum
   } else if(inherits(stratum, "formula")) {
         if(is.null(design))
            stop("\nArgument stratum was a formula, and argument design was not supplied.")
         temp <- match(as.character(stratum)[2], names(design), nomatch=0)
         if(temp == 0)
            stop("\nArgument stratum was a formula, and argument design did not contain a column \nmatching the name used in the formula.")
         stratum <- model.frame(stratum, design)
         stratum <- stratum[,1]
   }

# Obtain cluster codes

   if(is.null(cluster)) {
      if(!is.null(design$cluster))
         cluster <- design$cluster
   } else if(inherits(cluster, "formula")) {
         if(is.null(design))
            stop("\nArgument cluster was a formula, and argument design was not supplied.")
         temp <- match(as.character(cluster)[2], names(design), nomatch=0)
         if(temp == 0)
            stop("\nArgument cluster was a formula, and argument design did not contain a column \nmatching the name used in the formula.")
         cluster <- model.frame(cluster, design)
         cluster <- cluster[,1]
   }

# Obtain final adjusted stage one weights

   if(is.null(wgt1)) {
      if(!is.null(design$wgt1))
         wgt1 <- design$wgt1
   } else {
      if(inherits(wgt1, "formula")) {
         if(is.null(design))
            stop("\nArgument wgt1 was a formula, and argument design was not supplied.")
         temp <- match(as.character(wgt1)[2], names(design), nomatch=0)
         if(temp == 0)
            stop("\nArgument wgt1 was a formula, and argument design did not contain a column \nmatching the name used in the formula.")
         wgt1 <- model.frame(wgt1, design)
         wgt1 <- wgt1[,1]
      }
   }
   if(!is.numeric(wgt1) && !is.null(wgt1))
      stop("\nThe final adjusted stage one weights are not numeric values.")

# Obtain stage one x-coordinates for location

   if(is.null(xcoord1)) {
      if(!is.null(design$xcoord1))
         xcoord1 <- design$xcoord1
   } else {
      if(inherits(xcoord1, "formula")) {
         if(is.null(design))
            stop("\nArgument xcoord1 was a formula, and argument design was not supplied.")
         temp <- match(as.character(xcoord1)[2], names(design), nomatch=0)
         if(temp == 0)
            stop("\nArgument xcoord1 was a formula, and argument design did not contain a column \nmatching the name used in the formula.")
         xcoord1 <- model.frame(xcoord1, design)
         xcoord1 <- xcoord1[,1]
      }
   }
   if(!is.numeric(xcoord1) && !is.null(xcoord1))
      stop("\nThe stage one x-coordinates for location are not numeric values.")

# Obtain stage one y-coordinates for location

   if(is.null(ycoord1)) {
      if(!is.null(design$ycoord1))
         ycoord1 <- design$ycoord1
   } else {
      if(inherits(ycoord1, "formula")) {
         if(is.null(design))
            stop("\nArgument ycoord1 was a formula, and argument design was not supplied.")
         temp <- match(as.character(ycoord1)[2], names(design), nomatch=0)
         if(temp == 0)
            stop("\nArgument ycoord1 was a formula, and argument design did not contain a column \nmatching the name used in the formula.")
         ycoord1 <- model.frame(ycoord1, design)
         ycoord1 <- ycoord1[,1]
      }
   }
   if(!is.numeric(ycoord1) && !is.null(ycoord1))
      stop("\nThe stage one y-coordinates for location are not numeric values.")

# Obtain support values

   if(is.null(support)) {
      if(!is.null(design$support))
         support <- design$support
   } else if(inherits(support, "formula")) {
         if(is.null(design))
            stop("\nArgument support was a formula, and argument design was not supplied.")
         temp <- match(as.character(support)[2], names(design), nomatch=0)
         if(temp == 0)
            stop("\nArgument support was a formula, and argument design did not contain a column \nmatching the name used in the formula.")
         support <- model.frame(support, design)
         support <- support[,1]
   }

# Obtain size-weights

   if(is.null(swgt)) {
      if(!is.null(design$swgt))
         swgt <- design$swgt
   } else if(inherits(swgt, "formula")) {
         if(is.null(design))
            stop("\nArgument swgt was a formula, and argument design was not supplied.")
         temp <- match(as.character(swgt)[2], names(design), nomatch=0)
         if(temp == 0)
            stop("\nArgument swgt was a formula, and argument design did not contain a column \nmatching the name used in the formula.")
         swgt <- model.frame(swgt, design)
         swgt <- swgt[,1]
   }

# Obtain stage one size-weights

   if(is.null(swgt1)) {
      if(!is.null(design$swgt1))
         swgt1 <- design$swgt1
   } else if(inherits(swgt1, "formula")) {
         if(is.null(design))
            stop("\nArgument swgt1 was a formula, and argument design was not supplied.")
         temp <- match(as.character(swgt1)[2], names(design), nomatch=0)
         if(temp == 0)
            stop("\nArgument swgt1 was a formula, and argument design did not contain a column \nmatching the name used in the formula.")
         swgt1 <- model.frame(swgt1, design)
         swgt1 <- swgt1[,1]
   }

# As necessary, assign missing values to the design variables

   if(is.null(xcoord))
      xcoord <- rep(NA, nresp)
   if(is.null(ycoord))
      ycoord <- rep(NA, nresp)
   if(is.null(stratum))
      stratum <- rep(NA, nresp)
   if(is.null(cluster))
      cluster <- rep(NA, nresp)
   if(is.null(wgt1))
      wgt1 <- rep(NA, nresp)
   if(is.null(xcoord1))
      xcoord1 <- rep(NA, nresp)
   if(is.null(ycoord1))
      ycoord1 <- rep(NA, nresp)
   if(is.null(support))
      support <- rep(NA, nresp)
   if(is.null(swgt))
      swgt <- rep(NA, nresp)
   if(is.null(swgt1))
      swgt1 <- rep(NA, nresp)

# Check the design variables for length

   if(length(wgt) != nresp)
      stop("\nThe number of final adjusted weights does not match the number of site ID \nvalues.")
   if(length(xcoord) != nresp)
      stop("\nThe number of x-coordinates for location does not match the number of site ID \nvalues.")
   if(length(ycoord) != nresp)
      stop("\nThe number of y-coordinates for location does not match the number of site ID \nvalues.")
   if(length(stratum) != nresp)
      stop("\nThe number of stratum codes does not match the number of site ID values.")
   if(length(cluster) != nresp)
      stop("\nThe number of stage one sampling unit codes does not match the number of site \nID values.")
   if(length(wgt1) != nresp)
      stop("\nThe number of final adjusted stage one weights does not match the number of \nsite ID values.")
   if(length(xcoord1) != nresp)
      stop("\nThe number of stage one x-coordinates for location does not match the number of \nsite ID values.")
   if(length(ycoord1) != nresp)
      stop("\nThe number of stage one y-coordinates for location does not match the number of \nsite ID values.")
   if(length(support) != nresp)
      stop("\nThe number of support values does not match the number of site ID values.")
   if(length(swgt) != nresp)
      stop("\nThe number of size-weight values does not match the number of site ID values.")
   if(length(swgt1) != nresp)
      stop("\nThe number of stage one size-weight values does not match the number of site ID \nvalues.")

# Create the design data frame

   design <- data.frame(siteID=siteID, wgt=wgt, xcoord=xcoord, ycoord=ycoord,
      stratum=stratum, cluster=cluster, wgt1=wgt1, xcoord1=xcoord1,
      ycoord1=ycoord1, support=support, swgt=swgt, swgt1=swgt1)

# Check the sites data frame, the design data frame, the subpop data frame, and
# the data.cont data frame to assure valid contents

   temp <- dframe.check(sites, design, subpop, data.cat, data.cont, NULL,
      design.names)
   sites <- temp$sites
   design <- temp$design
   subpop <- temp$subpop
   data.cat <- temp$data.cat
   data.cont <- temp$data.cont

# Assign variables from the design data frame

      siteID <- design$siteID
      wgt <- design$wgt
      xcoord <- design$xcoord
      ycoord <- design$ycoord
      stratum <- design$stratum
      cluster <- design$cluster
      wgt1 <- design$wgt1
      xcoord1 <- design$xcoord1
      ycoord1 <- design$ycoord1
      support <- design$support
      swgt <- design$swgt
      swgt1 <- design$swgt1
   
# Check site IDs for repeat values and, as necessary, create unique site IDs and
# output a warning message

   temp <- sapply(split(siteID, siteID), length)
   if(any(temp > 1)) {
      warn.ind <- TRUE
      temp.str <- vecprint(names(temp)[temp > 1])
      warn <- paste("The following site ID values occur more than once among the values that were \ninput to the function:\n", temp.str, sep="")
      act <- "Unique site ID values were created.\n"
      warn.df <- rbind(warn.df, data.frame(func=I(fname),
         subpoptype=NA, subpop=NA, indicator=NA, stratum=NA, warning=I(warn),
         action=I(act)))
      siteID <- uniqueID(siteID)
      sites[,1] <- siteID
      subpop[,1] <- siteID
      if(!is.null(data.cat))
         data.cat[,1] <- siteID
      if(!is.null(data.cont))
         data.cont[,1] <- siteID
   }

# Determine whether the sample is stratified

   stratum.ind <- length(unique(stratum)) > 1

# If the sample is stratified, convert stratum to a factor, determine stratum 
# levels, and calculate number of strata

   if(stratum.ind) {
      stratum <- factor(stratum)
      stratum.levels <- levels(stratum)
      nstrata <- length(stratum.levels)
   } else {
      stratum.levels <- NULL
      nstrata <- NULL
   }

# Determine whether the sample has two stages

   cluster.ind <- length(unique(cluster)) > 1

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

# Ensure that popsize is a list

   if(!is.null(popsize) && !is.list(popsize))
      stop("\nThe popsize argument must be a list")

# Determine whether the population correction factor is to be used

   if(popcorrect && is.null(support))
      stop("\nThe logical value that indicates whether finite or continuous population \ncorrection factors should be employed during variance estimation was set to \nTRUE, but support values were not provided.")
   pcfactor.ind <- popcorrect

# Determine whether the sample uses size-weights

   if(sizeweight) {
      if(is.null(swgt))
      stop("\nThe logical value that indicates whether size-weights should be employed in the analysis was set to \nTRUE, but size-weights were not provided in the design data frame.")
   if(cluster.ind && is.null(swgt1))
      stop("\nThe sample has two stages and the logical value that indicates whether size- \nweights should be employed in the analysis was set to TRUE, but stage one \nsize-weights were not provided in the design data frame.")
   }
   swgt.ind <- sizeweight

# Check the vector of measurement error variance values for correct names

   if(!is.null(sigma)) {
      temp.names <- names(data.cont)[-1]
      if(length(sigma) != length(temp.names))
         stop("\nThe vector of measurement error variance values is not the correct length.")
      if(is.null(names(sigma)))
         stop("\nThe vector of measurement error variance values must be named.")
      temp <- match(temp.names, names(sigma), nomatch=0)
      if(any(temp == 0)) {
         temp.str <- vecprint(temp.names[temp == 0])
         stop(paste("\nThe following names for the response variables do not occur among the names for \nthe vector of measurement error variance values:\n", temp.str, sep=""))
      }
      temp <- match(names(sigma), temp.names, nomatch=0)
      if(any(temp == 0)) {
         temp.str <- vecprint(names(sigma)[temp == 0])
         stop(paste("\nThe following names for the vector of measurement error variance values do not \noccur among the names for the response variables:\n", temp.str, sep=""))
      }
      sigma <- sigma[temp]
   }

# Check the vector of values for variance of the measurement error variance for
# correct names

   if(!is.null(var.sigma)) {
      if(length(var.sigma) != length(temp.names))
         stop("\nThe vector of values for variance of the measurement error variance is not the \ncorrect length.")
      if(is.null(names(var.sigma)))
         stop("\nThe vector of values for variance of the measurement error variance must be \nnamed.")
      temp <- match(temp.names, names(var.sigma), nomatch=0)
      if(any(temp == 0)) {
         temp.str <- vecprint(temp.names[temp == 0])
         stop(paste("\nThe following names for the response variables do not occur among the names for \nthe vector of values for variance of the measurement error variance:\n", temp.str, sep=""))
      }
      temp <- match(names(var.sigma), temp.names, nomatch=0)
      if(any(temp == 0)) {
         temp.str <- vecprint(names(var.sigma)[temp == 0])
         stop(paste("\nThe following names for the vector of values for variance of the measurement \nerror variance do not occur among the names for the response variables:\n", temp.str, sep=""))
      }
      var.sigma <- var.sigma[temp]
   }

# Check for compatibility of input values

   temp <- input.check(length(siteID), wgt, sigma, var.sigma, xcoord, ycoord,
      stratum.ind, stratum, stratum.levels, nstrata, cluster.ind, cluster,
      cluster.levels, ncluster, wgt1, xcoord1, ycoord1, popsize, pcfactor.ind,
      pcfsize, N.cluster, stage1size, support, swgt.ind, swgt, swgt1, vartype,
      conf, pctval=pctval, subpop=subpop)
   popsize <- temp$popsize
   pcfsize <- temp$pcfsize
   N.cluster <- temp$N.cluster
   stage1size <- temp$stage1size

# If the sample has two stages, determine whether there are a sufficient number
# of sites in each stage one sampling unit to allow variance calculation

      if(cluster.ind) {
         ntypes <- dim(subpop)[2]
         typenames <- names(subpop)
         for(itype in 2:ntypes) {
            temp <- apply(table(cluster, subpop[,itype]) == 1, 2, sum)
            ind <- tapply(cluster, subpop[,itype], function(x)
               length(unique(x)))
            if(any(temp == ind)) {
               temp.str <- vecprint(names(temp)[temp == ind])
               warn.df <<- warn.df
               stop(paste("\nA variance estimate cannot be calculated since all of the stage one sampling \nunits contain a single stage two sampling unit for the following \nsubpopulation(s) of population ", typenames[itype], ":\n", temp.str, "\nEnter the following command to view the warning messages that were generated: \nwarnprnt() \n", sep=""))
            }
            if(any(temp > 0)) {
               temp.str <- vecprint(names(temp)[temp > 0])
               warn <- paste("Since they include one or more stage one sampling units with a single site, \nthe mean of the variance estimates for stage one sampling units with two or \nmore sites will be used as the variance estimate for stage one sampling units \nwith one site for the following subpopulation(s) of population\n", typenames[itype], ":\n", temp.str, sep="")
               act <- "The mean of the variance estimates will be used.\n"
               warn.df <- rbind(warn.df, data.frame(func=I(fname),
         subpoptype=NA, subpop=NA, indicator=NA, stratum=NA, warning=I(warn),
         action=I(act)))
            }
         }
      }

# If the sample was stratified and had two stages, then reset cluster to its 
# input value

   if(stratum.ind && cluster.ind)
      cluster <- cluster.in

# Recreate the design data frame

   design <- data.frame(siteID=siteID, wgt=wgt, xcoord=xcoord, ycoord=ycoord,
      stratum=stratum, cluster=cluster, wgt1=wgt1, xcoord1=xcoord1,
      ycoord1=ycoord1, support=support, swgt=swgt, swgt1=swgt1)
   names(design) <- design.names

# Create the output list

   psurvey.obj <- list(sites=sites, subpop=subpop, design=design,
      data.cat=data.cat, data.cont=data.cont, sigma=sigma, var.sigma=var.sigma,
      stratum.ind=stratum.ind, cluster.ind=cluster.ind, popsize=popsize,
      pcfactor.ind=pcfactor.ind, pcfsize=pcfsize, N.cluster=N.cluster,
      stage1size=stage1size, swgt.ind=swgt.ind, vartype=vartype, conf=conf,
      pctval=pctval)

# Assign class "spsurvey.analysis" to the output list

   class(psurvey.obj) <- "spsurvey.analysis"

# As necessary, output a message indicating that warning messages were generated
# during execution of the program

   if(warn.ind) {
      warn.df <<- warn.df
      if(nrow(warn.df) == 1)
         cat("During execution of the program, a warning message was generated.  The warning \nmessage is stored in a data frame named 'warn.df'.  Enter the following command \nto view the warning message: warnprnt()\n")
      else
         cat(paste("During execution of the program,", nrow(warn.df), "warning messages were generated.  The warning \nmessages are stored in a data frame named 'warn.df'.  Enter the following \ncommand to view the warning messages: warnprnt() \nTo view a subset of the warning messages (say, messages number 1, 3, and 5), \nenter the following command: warnprnt(m=c(1,3,5))\n"))
   }

# Return the list

   psurvey.obj
}
