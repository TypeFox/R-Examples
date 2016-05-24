#' Mark-Recapture Distance Sampling (mrds)
#'
#' This package implements mark-recapture distance sampling
#'     methods as described in D.L. Borchers, W. Zucchini and Fewster,
#'     R.M. (1988), "Mark-recapture models for line transect surveys",
#'     Biometrics 54: 1207-1220. and Laake, J.L. (1999) "Distance sampling
#'    with independent observers: Reducing bias from heterogeneity by
#'     weakening the conditional independence assumption." in Amstrup,
#'     G.W., Garner, S.C., Laake, J.L., Manly, B.F.J., McDonald, L.L. and
#'     Robertson, D.G. (eds) "Marine mammal survey and assessment
#'     methods", Balkema, Rotterdam: 137-148 and Borchers, D.L., Laake,
#'     J.L., Southwell, C. and Paxton, C.L.G. "Accommodating unmodelled
#'     heterogeneity in double-observer distance sampling surveys". 2006.
#'     Biometrics 62:372-378.)
#'
#' Further information on distance sampling methods and example code is available at \url{http://distancesampling.org/R/}.
#'
#' For help with distance sampling and this package, there is a Google Group \url{https://groups.google.com/forum/#!forum/distance-sampling}.
#'
#' @name mrds-package
#' @aliases mrds-package mrds
#' @docType package
#' @author Jeff Laake <jeff.laake@@noaa.gov>,
#'         David Borchers <dlb@@mcs.st-and.ac.uk>,
#'         Len Thomas <len@@mcs.st-and.ac.uk>,
#'         David L. Miller <dave@@ninepointeightone.net>,
#'         Jon Bishop <jonb@@mcs.st-and.ac.uk>
#' @keywords package
#'
NULL



#' Golf tee data used in chapter 6 of Advanced Distance Sampling examples
#'
#' Double platform data collected in a line transect survey of golf tees by 2
#' observers at St. Andrews. Field sex was actually colour of the golf tee: 0 -
#' green; 1 - yellow. Exposure was either low (0) or high(1) depending on
#' height of tee above the ground. size was the number of tees in an observed
#' cluster.
#'
#'
#' @name book.tee.data
#' @docType data
#' @format The format is: List of 4 $ book.tee.dataframe:'data.frame': 324 obs.
#'   of 7 variables: ..$ object : num [1:324] 1 1 2 2 3 3 4 4 5 5 ...  ..$
#'   observer: Factor w/ 2 levels "1","2": 1 2 1 2 1 2 1 2 1 2 ...  ..$
#'   detected: num [1:324] 1 0 1 0 1 0 1 0 1 0 ...  ..$ distance: num [1:324]
#'   2.68 2.68 3.33 3.33 0.34 0.34 2.53 2.53 1.46 1.46 ...  ..$ size : num
#'   [1:324] 2 2 2 2 1 1 2 2 2 2 ...  ..$ sex : num [1:324] 1 1 1 1 0 0 1 1 1 1
#'   ...  ..$ exposure: num [1:324] 1 1 0 0 0 0 1 1 0 0 ...  $ book.tee.region
#'   :'data.frame': 2 obs. of 2 variables: ..$ Region.Label: Factor w/ 2 levels
#'   "1","2": 1 2 ..$ Area : num [1:2] 1040 640 $ book.tee.samples
#'   :'data.frame': 11 obs. of 3 variables: ..$ Sample.Label: num [1:11] 1 2 3
#'   4 5 6 7 8 9 10 ...  ..$ Region.Label: Factor w/ 2 levels "1","2": 1 1 1 1
#'   1 1 2 2 2 2 ...  ..$ Effort : num [1:11] 10 30 30 27 21 12 23 23 15 12 ...
#'   $ book.tee.obs :'data.frame': 162 obs. of 3 variables: ..$ object : int
#'   [1:162] 1 2 3 21 22 23 24 59 60 61 ...  ..$ Region.Label: int [1:162] 1 1
#'   1 1 1 1 1 1 1 1 ...  ..$ Sample.Label: int [1:162] 1 1 1 1 1 1 1 1 1 1 ...
#' @keywords datasets
NULL

#' Pronghorn aerial survey data from Wyoming
#'
#' Detections of pronghorn from fixed-wing aerial surveys in Southeastern
#' Wyoming using four angular bins defined by strut marks. Illustrates data
#' where altitude above ground level (AGL) varies during the survey.
#'
#' Each record is an observed cluster of pronghorn.  The data provide the
#' stratum for the observation, the direction of travel, the AGL at the time of
#' the observation, the angular bin which contained the center of the pronghorn
#' cluster(group), and the number of pronghorn in the group. The angular bins
#' were defined by a combination of two window and five wing strut marks to
#' define bin cutpoints for perpendicular ground distances of 0-65, 65-90,
#' 90-115, 115-165 and 165-265 meters when the plane is 300' (91.4 meters)
#' above ground level. The inner band is considered a blind region due to
#' obstruction of view beneath the plane; thus th the line is offset 65 meters
#' from underneath the plane.
#'
#' @name pronghorn
#' @docType data
#' @format A data frame with 660 observations on the following 5 variables.
#'   \describe{ \item{STRATUM}{a numeric vector}
#'   \item{direction}{a factor with levels \code{N} \code{S}
#'   representing the survey direction} \item{AGL}{height above ground
#'   level} \item{Band}{a factor with levels \code{A} \code{B} \code{C}
#'   \code{D} which represent angular bands between breaks at
#'   35.42,44.56,51.52,61.02,70.97 degrees.  These angles were set based on
#'   selected distance bins based on the target AGL.}
#'   \item{cluster}{number of pronghorn in the observed cluster} }
#' @references Laake, J., R. J. Guenzel, J. L. Bengtson, P. Boveng, M. Cameron,
#'   and M. B. Hanson. 2008.  Coping with variation in aerial survey protocol
#'   for line-transect sampling. Wildlife Research 35:289-298.
#' @source Data provided courtesy of Rich Guenzel of Wyoming Game and Fish.
#' @keywords datasets
NULL


#' Wooden stake data from 1977 survey
#'
#' Multiple surveys by different observers of a single 1km transect containing
#' 150 wooden stakes placed randomly throughout a 40 m strip (20m on either
#' side).
#'
#' @name stake77
#' @docType data
#' @format A data frame with 150 observations on the following 10 variables.
#'   \describe{ \item{StakeNo}{unique number for each stake 1-150}
#'   \item{PD}{perpendicular distance at which the stake was placed
#'   from the line} \item{Obs1}{0/1 whether missed/seen by observer 1}
#'   \item{Obs2}{0/1 whether missed/seen by observer 2}
#'   \item{Obs3}{0/1 whether missed/seen by observer 3}
#'   \item{Obs4}{0/1 whether missed/seen by observer 4}
#'   \item{Obs5}{0/1 whether missed/seen by observer 5}
#'   \item{Obs6}{0/1 whether missed/seen by observer 6}
#'   \item{Obs7}{0/1 whether missed/seen by observer 7}
#'   \item{Obs8}{0/1 whether missed/seen by observer 8} }
#' @references Burnham, K. P., D. R. Anderson, and J. L. Laake. 1980.
#'   Estimation of Density from Line Transect Sampling of Biological
#'   Populations. Wildlife Monographs:7-202.
#' @source Laake, J. 1978. Line transect estimators robust to animal movement.
#'   M.S. Thesis. Utah State University, Logan, Utah. 55p.
#' @keywords datasets
#' @examples
#' \donttest{
#' data(stake77)
#' # Extract functions for stake data and put in the mrds format
#' extract.stake <- function(stake,obs){
#'   extract.obs <- function(obs){
#'     example <- subset(stake,eval(parse(text=paste("Obs",obs,"==1",sep=""))),
#'                       select="PD")
#'     example$distance <- example$PD
#'     example$object <- 1:nrow(example)
#'     example$PD <- NULL
#'     return(example)
#'   }
#'   if(obs!="all"){
#'     return(extract.obs(obs=obs))
#'   }else{
#'     example <- NULL
#'     for(i in 1:(ncol(stake)-2)){
#'       df <- extract.obs(obs=i)
#'       df$person <- i
#'       example <- rbind(example,df)
#'     }
#'     example$person <- factor(example$person)
#'     example$object <- 1:nrow(example)
#'     return(example)
#'   }
#' }
#' extract.stake.pairs <- function(stake,obs1,obs2,removal=FALSE){
#'   obs1 <- paste("Obs",obs1,sep="")
#'   obs2 <- paste("Obs",obs2,sep="")
#'   example <- subset(stake,eval(parse(text=paste(obs1,"==1 |",obs2,"==1 ",
#'                                        sep=""))),select=c("PD",obs1,obs2))
#'   names(example) <- c("distance","obs1","obs2")
#'   detected <- c(example$obs1,example$obs2)
#'   example <- data.frame(object   = rep(1:nrow(example),2),
#'                         distance = rep(example$distance,2),
#'                         detected = detected,
#'                         observer = c(rep(1,nrow(example)),rep(2,nrow(example))))
#'   if(removal) example$detected[example$observer==2] <- 1
#'   return(example)
#' }
#' # extract data for observer 1 and fit a single observer model
#' stakes <- extract.stake(stake77,1)
#' ds.model <- ddf(dsmodel = ~mcds(key = "hn", formula = ~1), data = stakes,
#'                 method = "ds", meta.data = list(width = 20))
#' plot(ds.model,breaks=seq(0,20,2),showpoints=TRUE)
#' ddf.gof(ds.model)
#' 
#' # extract data from observers 1 and 3 and fit an io model
#' stkpairs <- extract.stake.pairs(stake77,1,3,removal=FALSE)
#' io.model <- ddf(dsmodel = ~mcds(key = "hn", formula=~1),
#'                 mrmodel=~glm(formula=~distance),
#'                 data = stkpairs, method = "io")
#' summary(io.model)
#' par(mfrow=c(3,2))
#' plot(io.model,breaks=seq(0,20,2),showpoints=TRUE,new=FALSE)
#' dev.new()
#' ddf.gof(io.model)
#' }
NULL


#' Wooden stake data from 1978 survey
#'
#' Multiple surveys by different observers of a single 1km transect containing
#' 150 wooden stakes placed based on expected uniform distribution throughout a
#' 40 m strip (20m on either side).
#'
#' The 1997 survey was based on a single realization of a uniform distribution.
#' Because it was a single transect and there was no randomization of the
#' distances for each survey, we repeated the experiment and used distances
#' that provided a uniform distribution but randomly sorted the positions along
#' the line so there was no pattern obvious to the observer.
#'
#' @name stake78
#' @docType data
#' @format A data frame with 150 observations on the following 13 variables.
#'   \describe{ \item{StakeNo}{unique number for each stake 1-150}
#'   \item{PD}{perpendicular distance at which the stake was placed
#'   from the line} \item{Obs1}{0/1 whether missed/seen by observer 1}
#'   \item{Obs2}{0/1 whether missed/seen by observer 2}
#'   \item{Obs3}{0/1 whether missed/seen by observer 3}
#'   \item{Obs4}{0/1 whether missed/seen by observer 4}
#'   \item{Obs5}{0/1 whether missed/seen by observer 5}
#'   \item{Obs6}{0/1 whether missed/seen by observer 6}
#'   \item{Obs7}{0/1 whether missed/seen by observer 7}
#'   \item{Obs8}{0/1 whether missed/seen by observer 8}
#'   \item{Obs9}{0/1 whether missed/seen by observer 9}
#'   \item{Obs10}{0/1 whether missed/seen by observer 10}
#'   \item{Obs11}{0/1 whether missed/seen by observer 11} }
#' @references Burnham, K. P., D. R. Anderson, and J. L. Laake. 1980.
#'   Estimation of Density from Line Transect Sampling of Biological
#'   Populations. Wildlife Monographs:7-202.
#' @source Laake, J. 1978. Line transect estimators robust to animal movement.
#'   M.S. Thesis. Utah State University, Logan, Utah. 55p.
#' @keywords datasets
#' @examples
#' \donttest{
#' data(stake78)
#' data(stake77)
#' # compare distribution of distances for all stakes
#' hist(stake77$PD)
#' hist(stake78$PD)
#' # Extract stake data and put in the mrds format for model fitting.
#' extract.stake <- function(stake,obs){
#'   extract.obs <- function(obs){
#'     example <- subset(stake,eval(parse(text=paste("Obs",obs,"==1",sep=""))),
#'                       select="PD")
#'     example$distance <- example$PD
#'     example$object <- 1:nrow(example)
#'     example$PD <- NULL
#'     return(example)
#'   }
#'   if(obs!="all"){
#'      return(extract.obs(obs=obs))
#'   }else{
#'     example <- NULL
#'     for(i in 1:(ncol(stake)-2)){
#'       df <- extract.obs(obs=i)
#'       df$person <- i
#'       example <- rbind(example,df)
#'     }
#'     example$person <- factor(example$person)
#'     example$object <- 1:nrow(example)
#'     return(example)
#'   }
#' }
#' extract.stake.pairs <- function(stake,obs1,obs2,removal=FALSE){
#'   obs1 <- paste("Obs",obs1,sep="")
#'   obs2 <- paste("Obs",obs2,sep="")
#'   example <- subset(stake,eval(parse(text=paste(obs1,"==1 |",obs2,"==1 ",
#'                                      sep=""))), select=c("PD",obs1,obs2))
#'   names(example) <- c("distance","obs1","obs2")
#'   detected <- c(example$obs1,example$obs2)
#'   example <- data.frame(object=rep(1:nrow(example),2),
#'                         distance=rep(example$distance,2),
#'                         detected = detected,
#'                         observer=c(rep(1,nrow(example)),
#'                                    rep(2,nrow(example))))
#'   if(removal) example$detected[example$observer==2] <- 1
#'   return(example)
#' }
#'
#' # extract data for observer 10 and fit a single observer model
#' stakes <- extract.stake(stake78,10)
#' ds.model <- ddf(dsmodel = ~mcds(key = "hn", formula = ~1), data = stakes,
#'                 method = "ds", meta.data = list(width = 20))
#' plot(ds.model,breaks=seq(0,20,2),showpoints=TRUE)
#' ddf.gof(ds.model)
#'
#' # extract data from observers 5 and 7 and fit an io model
#' stkpairs <- extract.stake.pairs(stake78,5,7,removal=FALSE)
#' io.model <- ddf(dsmodel = ~mcds(key = "hn", formula=~1),
#'                 mrmodel=~glm(formula=~distance),
#'                 data = stkpairs, method = "io")
#' summary(io.model)
#' par(mfrow=c(3,2))
#' plot(io.model,breaks=seq(0,20,2),showpoints=TRUE,new=FALSE)
#' ddf.gof(io.model)
#' }
#'
NULL

#' Single observer point count data example from Distance
#'
#' Single observer point count data example from Distance
#'
#' @name ptdata.distance
#' @docType data
#' @format The format is 144 obs of 6 variables:
#'   distance: numeric distance from center
#'   observer: Factor w/ 2 levels "1","2": 1 2 1 2 1 2 1 2 1 2 ...
#'   detected: numeric 0/1
#'   object: sequential object number
#'   Sample.Label: point label
#'   Region.Label: single region label
#' @keywords datasets
#' @examples
#' \donttest{
#' data(ptdata.distance)
#' xx <- ddf(dsmodel = ~cds(key="hn", formula = ~1), data = ptdata.distance,
#'           method = "ds", meta.data = list(point=TRUE))
#' summary(xx)
#' plot(xx,main="Distance point count data")
#' ddf.gof(xx)
#' Regions <- data.frame(Region.Label=1,Area=1)
#' Samples <- data.frame(Sample.Label=1:30,
#'                       Region.Label=rep(1,30),
#'                       Effort=rep(1,30))
#' print(dht(xx,sample.table=Samples,region.table=Regions))
#' }
NULL



#' Simulated single observer point count data
#'
#' Simulated single observer point count data with detection p(0)=1;
#' hn sigma=30; w=100
#'
#' @name ptdata.single
#' @docType data
#' @format The format is 341 obs of 4 variables: ..$
#'   distance: numeric distance from center $
#'   observer: Factor w/ 2 levels "1","2": 1 2 1 2 1 2 1 2 1 2 ...  ..$
#'   detected: numeric 0/1  $ object : sequential object number
#' @keywords datasets
#' @examples
#' \donttest{
#' data(ptdata.single)
#' xx=ddf(dsmodel = ~cds(key="hn", formula = ~1), data = ptdata.single,
#'          method = "ds", meta.data = list(point=TRUE))
#' summary(xx)
#' plot(xx,main="Simulated point count data")
#' }
NULL

#' Simulated dual observer point count data
#'
#' Simulated dual observer point count data with detection p(0)=0.8;
#' hn sigma=30; w=100 for both observers with dependency y>0, gamma=0.1
#'
#' @name ptdata.dual
#' @docType data
#' @format The format is 420 obs of 6 variables:
#' distance: numeric distance from center
#' observer: Factor w/ 2 levels "1","2": 1 2 1 2 1 2 1 2 1 2 ...
#' detected: numeric 0/1
#' person: Factor with 2 levels A,B
#' pair: Factor with 2 levels "AB" BA" $
#' object : sequential object number
#' @keywords datasets
#' @examples
#' \donttest{
#' data(ptdata.dual)
#' xx <- ddf(mrmodel=~glm(formula=~distance),
#'           dsmodel = ~cds(key="hn", formula = ~1),
#'           data = ptdata.dual, method = "io", meta.data = list(point=TRUE))
#' summary(xx)
#' plot(xx,main="Simulated point count data")
#' }
NULL


#' Simulated removal observer point count data
#'
#' Simulated removal observer point count data with detection p(0)=0.8;
#' hn sigma=30; w=100 for both observers with dependency y>0, gamma=0.1
#'
#' @name ptdata.removal
#' @docType data
#' @format The format is 408 obs of 6 variables:
#'  distance: numeric distance from center
#'  observer: Factor w/ 2 levels "1","2": 1 2 1 2 1 2 1 2 1 2 ...
#'  detected: numeric 0/1
#'  person: Factor with 2 levels A,B
#'  pair: Factor with 2 levels "AB" BA"
#'  object: sequential object number
#' @keywords datasets
#' @examples
#' \donttest{
#' data(ptdata.removal)
#' xx <- ddf(mrmodel=~glm(formula=~distance),
#'           dsmodel = ~cds(key="hn", formula = ~1),
#'           data = ptdata.removal, method = "rem",
#'           meta.data = list(point=TRUE))
#' summary(xx)
#' plot(xx,main="Simulated point count data")
#' }
NULL

#' Golden-cheeked warbler mark-recapture distance sampling analysis
#'
#' These data represent avian point count surveys conducted at 453 point sample survey locations on the 24,000 (approx) live-fire region of Fort Hood
#' in central Texas.  Surveys were conducted by independent double observers (2 per survey occasion) and as such we had a maximum of 3 paired survey histories, giving a maximum of
#' 6 sample occasions (see MacKenzie et al. 2006, MacKenzie and Royle 2005, and Laake et al. 2011 for various sample survey design details).  At each point, we surveyed for 5 minutes (technically broken into
#' 3 time intervals of 2, 2, and 1 minutes; not used here) and we noted detections by each observer and collected distance to each observation within a set of distance bins (0-50, 50-100m; Laake et al. 2011)
#' of the target species (Golden-cheeked warblers in this case) for each surveyor.  Our primary focus was to use mark-recapture distance sampling methods to estimate density of
#' Golden-cheeked warblers, and to estimate detection rates for the mark-recapture, distance, and composite model.
#'
#' @details  In addition to detailing the analysis used by Collier et al. (2013, In Review), this example documents the use of \code{mrds} for avian point count surveys and shows how density models
#' can be incorporated with occupancy models to develop spatially explicit density surface maps. For those that are interested, for the distance sampling portion of our analysis, we
#' used both conventional distance sampling (\code{cds}) and multiple covariate distance sampling (\code{mcds}) with uniform and half-normal key functions.  For the mark-recapture portion of
#' our analysis, we tended to use covariates for distance (median bin width), observer, and date of survey (days since 15 March 2011).
#'
#' We combined our \code{mrds} density estimates via a Horvitz-Thompson styled estimator with the resource selection function gradient developed
#' in Farrell et al. (2013) and estimated density on an ~3.14ha hexagonal grid across our study area, which provided a density gradient for Fort Hood.  Because there
#' was considerable data manipulation needed for each analysis to structure the data appropriately for use in \code{mrds}, rather than wrap each analysis in a single function, we have provided
#' both the Golden-cheeked warbler and Black-capped vireo analyses in their full detail.  The primary differences you will see will be changes to
#' model structures and model outputs between the two species.
#'
#' @name lfgcwa
#' @docType data
#' @format The format is a data frame with the following covariate metrics.
#' \describe{\item{PointID}{Unique identifier for each sample location; locations are the same for both species}
#' \item{VisitNumber}{Visit number to the point}
#' \item{Species}{Species designation, either Golden-cheeked warbler (GW) or Black-capped Vireo (BV)}
#' \item{Distance}{Distance measure, which is either NA (representing no detection), or the median of the binned detection distances}
#' \item{PairNumber}{ID value indicating which observers were paired for that sampling occasion}
#' \item{Observer}{Observer ID, either primary(1), or secondary (2)}
#' \item{Detected}{Detection of a bird, either 1 = detected, or 0 = not detected}
#' \item{Date}{Date of survey since 15 March 2011, numeric value}
#' \item{Pred}{Predicted occupancy value for that survey hexagon based on Farrell et al. (2013)}
#' \item{Category}{Region.Label categorization, see R package \code{mrds} help file for details on data structure}
#' \item{Effort}{Amount of survey effort at the point}
#' \item{Day}{Number of days since 15 March 2011, numeric value}
#' \item{ObjectID}{Unique ID for each paired observations} }
#'
#' @references Farrell, S.F., B.A. Collier, K.L. Skow, A.M. Long, A.J. Campomizzi, M.L. Morrison, B. Hays, and R.N. Wilkins. 2013. Using LiDAR-derived structural vegetation characteristics to develop high-resolution, small-scale, species distribution models for conservation planning. Ecosphere 43(3): 42. http://dx.doi.org/10.1890/ES12-000352.1
#' @references Laake, J.L., B.A. Collier, M.L. Morrison, and R.N. Wilkins. 2011. Point-based mark recapture distance sampling. Journal of Agricultural, Biological and Environmental Statistics 16: 389-408.
#' @references Collier, B.A., S.L. Farrell, K.L. Skow, A.M. Long, A.J. Campomizzi, K.B. Hays, J.L. Laake, M.L. Morrison, and R.N. Wilkins. 2013. Spatially explicit density of endangered avian species in a disturbed landscape.  Auk, In Review.
#' @keywords datasets
#' @author Bret Collier and Jeff Laake
#' @examples
#'
#' \dontrun{
#' data(lfgcwa)
#' xy <- cut(lfgcwa$Pred, c(-0.0001, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1),
#'  labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
#' x <- data.frame(lfgcwa, New=xy)
#'
#' # Note that I scaled the individual covariate of day-helps with
#' # convergence issues
#' bird.data <- data.frame(object=x$ObjectID, observer=x$Observer,
#'                         detected=x$Detected, distance=x$Distance,
#'                         Region.Label=x$New, Sample.Label=x$PointID,
#'                         Day=(x$Day/max(x$Day)))
#'
#' # make observer a factor variable
#' bird.data$observer=factor(bird.data$observer)
#'
#' # Jeff Laake suggested this snippet to quickly create distance medians
#' # which adds bin information to the \code{bird.data} dataframe
#'
#' bird.data$distbegin=0
#' bird.data$distend=100
#' bird.data$distend[bird.data$distance==12.5]=50
#' bird.data$distbegin[bird.data$distance==37.5]=0
#' bird.data$distend[bird.data$distance==37.5]=50
#' bird.data$distbegin[bird.data$distance==62.5]=50
#' bird.data$distend[bird.data$distance==62.5]=100
#' bird.data$distbegin[bird.data$distance==87.5]=50
#' bird.data$distend[bird.data$distance==87.5]=100
#'
#' # Removed all survey points with distance=NA for a survey event;
#' # hence no observations for use in \code{ddf()} but needed later
#' bird.data=bird.data[complete.cases(bird.data),]
#'
#' # Manipulations on full dataset for various data.frame creation
#' # for use in density estimation using \code{dht()}
#'
#' # Samples dataframe
#' xx <- x
#' x <- data.frame(PointID=x$PointID, Species=x$Species,
#'                 Category=x$New, Effort=x$Effort)
#' x <- x[!duplicated(x$PointID),]
#' point.num <- table(x$Category)
#' samples <- data.frame(PointID=x$PointID, Region.Label=x$Category,
#'                       Effort=x$Effort)
#' final.samples=data.frame(Sample.Label=samples$PointID,
#'                          Region.Label=samples$Region.Label,
#'                          Effort=samples$Effort)
#'
#' # obs dataframe
#' obs <- data.frame(ObjectID=xx$ObjectID, PointID=xx$PointID)
#' # used to get Region and Sample assigned to ObjectID
#' obs <- merge(obs, samples, by=c("PointID", "PointID"))
#' obs <- obs[!duplicated(obs$ObjectID),]
#' obs <- data.frame(object=obs$ObjectID, Region.Label=obs$Region.Label,
#'                   Sample.Label=obs$PointID)
#'
#' #Region.Label dataframe
#' region.data <- data.frame(Region.Label=c(1,2,3,4,5,6,7,8,9),
#'                           Area=c(point.num[1]*3.14,
#'                                  point.num[2]*3.14,
#'                                  point.num[3]*3.14,
#'                                  point.num[4]*3.14,
#'                                  point.num[5]*3.14,
#'                                  point.num[6]*3.14,
#'                                  point.num[7]*3.14,
#'                                  point.num[8]*3.14,
#'                                  point.num[9]*3.14))
#'
#' # Candidate Models
#'
#' GW1=ddf(
#'    dsmodel=~cds(key="unif", adj.series="cos", adj.order=1,adj.scale="width"),
#'    mrmodel=~glm(~distance),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE,point=TRUE,width=100,breaks=c(0,50,100)))
#' GW2=ddf(
#'    dsmodel=~cds(key="unif", adj.series="cos", adj.order=1,adj.scale="width"),
#'    mrmodel=~glm(~distance+observer),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE,point=TRUE,width=100,breaks=c(0,50,100)))
#' GW3=ddf(
#'    dsmodel=~cds(key="unif", adj.series="cos", adj.order=1,adj.scale="width"),
#'    mrmodel=~glm(~distance*observer),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE,point=TRUE,width=100,breaks=c(0,50,100)))
#' GW4=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE,point=TRUE,width=100,breaks=c(0,50,100)))
#' GW4FI=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance),
#'    data=bird.data,
#'    method="io.fi",
#'    meta.data=list(binned=TRUE,point=TRUE,width=100,breaks=c(0,50,100)))
#' GW5=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance+observer),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' GW5FI=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance+observer),
#'    data=bird.data,
#'    method="io.fi",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' GW6=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance*observer),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' GW6FI=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance*observer),
#'    data=bird.data,
#'    method="io.fi",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' GW7=ddf(
#'    dsmodel=~cds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance*Day),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' GW7FI=ddf(
#'    dsmodel=~cds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance*Day),
#'    data=bird.data,
#'    method="io.fi",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' GW8=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance*observer*Day),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' GW8FI=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance*observer*Day),
#'    data=bird.data,
#'    method="io.fi",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#'#GWDS=ddf(
#'#   dsmodel=~mcds(key="hn",formula=~1),
#'#   data=bird.data,
#'#   method="ds",
#'#   meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#'
#'
#'
#' #### GCWA Summary Metrics
#'
#' #AIC table building code, not exactly elegant, but I did not
#' want to add more package dependencies
#' AIC = c(GW1$criterion, GW2$criterion, GW3$criterion, GW4$criterion,
#'         GW4FI$criterion, GW5$criterion, GW5FI$criterion,
#'         GW6$criterion, GW6FI$criterion, GW7$criterion, GW7FI$criterion,
#'         GW8$criterion, GW8FI$criterion)
#'
#' #creates a set of row names for me to check my grep() call below
#' rn <- c("GW1", "GW2", "GW3", "GW4", "GW4FI", "GW5", "GW5FI", "GW6",
#'         "GW6FI", "GW7","GW7FI", "GW8", "GW8FI")
#'
#' # number of parameters for each model
#' k <- c(length(GW1$par), length(GW2$par), length(GW3$par), length(GW4$par),
#'        length(GW4FI$par), length(GW5$par), length(GW5FI$par),
#'        length(GW6$par), length(GW6FI$par), length(GW7$par),
#'        length(GW7FI$par), length(GW8$par), length(GW8FI$par))
#'
#' # build AIC table and
#' AIC.table <- data.frame(AIC = AIC, rn=rn, k=k, dAIC = abs(min(AIC)-AIC),
#'                         likg = exp(-.5*(abs(min(AIC)-AIC))))
#' # row.names(AIC.table)=grep("GW", ls(), value=TRUE)
#' AIC.table <- AIC.table[with(AIC.table, order(-likg, -dAIC, AIC, k)),]
#' AIC.table <- data.frame(AIC.table, wi=AIC.table$likg/sum(AIC.table$likg))
#' AIC.table
#'
#' # Model average N_hat_covered estimates
#' # not very clean, but I wanted to show full process, need to use
#' # collect.models and model.table here
#'
#' estimate <- c(GW1$Nhat, GW2$Nhat, GW3$Nhat, GW4$Nhat, GW4FI$Nhat,
#'               GW5$Nhat, GW5FI$Nhat, GW6$Nhat, GW6FI$Nhat, GW7$Nhat,
#'               GW7FI$Nhat, GW8$Nhat, GW8FI$Nhat)
#' AIC.values <- AIC
#'
#' # Nhat.se is calculated in mrds:::summary.io, not in ddf(), so
#' # it takes a bit to pull out
#' std.err <- c(summary(GW1)$Nhat.se, summary(GW2)$Nhat.se,
#'              summary(GW3)$Nhat.se,summary(GW4)$Nhat.se,
#'              summary(GW4FI)$Nhat.se, summary(GW5)$Nhat.se,
#'              summary(GW5FI)$Nhat.se, summary(GW6)$Nhat.se,
#'              summary(GW6FI)$Nhat.se, summary(GW7)$Nhat.se,
#'              summary(GW7FI)$Nhat.se,summary(GW8)$Nhat.se,
#'              summary(GW8FI)$Nhat.se)
#'}
#' \dontrun{
#' #Not Run
#' #requires RMark
#' library(RMark)
#' #uses model.average structure to model average real abundance estimates for
#' #covered area of the surveys
#' mmi.list=list(estimate=estimate, AIC=AIC.values, se=std.err)
#' model.average(mmi.list, revised=TRUE)
#'
#'#Not Run
#'#Best Model FI
#'#best.modelFI=AIC.table[1,]
#'#best.model
#'#Best Model PI
#'#best.modelPI=AIC.table[2,]
#'#best.modelPI
#'
#'#Not Run
#'#summary(GW7FI, se=TRUE)
#'#summary(GW7, se=TRUE)
#'
#'#Not Run
#'#GOF for models
#'#ddf.gof(GW7, breaks=c(0,50,100))
#'
#'#Not Run
#'#Density estimation across occupancy categories
#'#out.GW=dht(GW7, region.data, final.samples, obs, se=TRUE,
#'            options=list(convert.units=.01))
#'
#' #Plots--Not Run
#' #Composite Detection Function examples
#' #plot(GW7, which=3, showpoints=FALSE, angle=0, density=0,
#' #     col="black", lwd=3, main="Golden-cheeked Warbler",
#' #     xlab="Distance (m)", las=1, cex.axis=1.25, cex.lab=1.25)
#'
#' #Conditional Detection Function
#' #dd=expand.grid(distance=0:100,Day=(4:82)/82)
#' #dmat=model.matrix(~distance*Day,dd)
#' #dd$p=plogis(model.matrix(~distance*Day,dd)%*%coef(GW7$mr)$estimate)
#' #dd$Day=dd$Day*82
#' #with(dd[dd$Day==12,],plot(distance,p,ylim=c(0,1), las=1,
#' # ylab="Detection probability", xlab="Distance (m)",
#' #  type="l",lty=1, lwd=3, bty="l", cex.axis=1.5, cex.lab=1.5))
#' #with(dd[dd$Day==65,],lines(distance,p,lty=2, lwd=3))
#' #ch=paste(bird.data$detected[bird.data$observer==1],
#' #         bird.data$detected[bird.data$observer==2],
#' #         sep="")
#' #tab=table(ch,cut(82*bird.data$Day[bird.data$observer==1],c(0,45,83)),
#' # cut(bird.data$distance[bird.data$observer==1],c(0,50,100)))
#' #tabmat=cbind(colMeans(rbind(tab[3,,1]/colSums(tab[2:3,,1],
#' #                            tab[3,,1]/colSums(tab[c(1,3),,1])))),
#' #             colMeans(rbind(tab[3,,2]/colSums(tab[2:3,,2],
#' #                            tab[3,,2]/colSums(tab[c(1,3),,2])))))
#' #colnames(tabmat)=c("0-50","51-100")
#' #points(c(25,75),tabmat[1,],pch=1, cex=1.5)
#' #points(c(25,75),tabmat[2,],pch=2, cex=1.5)
#'
#' # Another alternative plot using barplot instead of points
#' # (this is one in paper)
#'
#' #ch=paste(bird.data$detected[bird.data$observer==1],
#' #         bird.data$detected[bird.data$observer==2],
#' #sep="")
#' #tab=table(ch,cut(82*bird.data$Day[bird.data$observer==1],c(0,45,83)),
#' # cut(bird.data$distance[bird.data$observer==1],c(0,50,100)))
#' #tabmat=cbind(colMeans(rbind(tab[3,,1]/colSums(tab[2:3,,1],
#' #                            tab[3,,1]/colSums(tab[c(1,3),,1])))),
#' #colMeans(rbind(tab[3,,2]/colSums(tab[2:3,,2],
#' #               tab[3,,2]/colSums(tab[c(1,3),,2])))))
#' #colnames(tabmat)=c("0-50","51-100")
#' #par(mfrow=c(2, 1), mai=c(1,1,1,1))
#' #with(dd[dd$Day==12,],
#' #     plot(distance,p,ylim=c(0,1), las=1,
#' #          ylab="Detection probability", xlab="",
#' #          type="l",lty=1, lwd=4, bty="l", cex.axis=1.5, cex.lab=1.5))
#' #segments(0, 0, .0, tabmat[1,1], lwd=3)
#' #segments(0, tabmat[1,1], 50, tabmat[1,1], lwd=4)
#' #segments(50, tabmat[1,1], 50, 0, lwd=4)
#' #segments(50, tabmat[1,2], 100, tabmat[1,2], lwd=4)
#' #segments(0, tabmat[1,1], 50, tabmat[1,1], lwd=4)
#' #segments(100, tabmat[1,2], 100, 0, lwd=4)
#' #mtext("a",line=-1, at=90)
#' #with(dd[dd$Day==65,],
#' #     plot(distance,p,ylim=c(0,1), las=1, ylab="Detection probability",
#' #          xlab="Distance", type="l",lty=1,
#' #          lwd=4, bty="l", cex.axis=1.5, cex.lab=1.5))
#' #segments(0, 0, .0, tabmat[2,1], lwd=4)
#' #segments(0, tabmat[2,1], 50, tabmat[2,1], lwd=4)
#' #segments(50, tabmat[2,1], 50, 0, lwd=4)
#' #segments(50, tabmat[2,2], 50, tabmat[2,1], lwd=4)
#' #segments(50, tabmat[2,2], 100, tabmat[2,2], lwd=4)
#' #segments(100, tabmat[2,2], 100, 0, lwd=4)
#' #mtext("b",line=-1, at=90)
#' }
#'
NULL

#' Black-capped vireo mark-recapture distance sampling analysis
#'
#' These data represent avian point count surveys conducted at 453 point sample survey locations on the 24,000 (approx) live-fire region of Fort Hood
#' in central Texas.  Surveys were conducted by independent double observers (2 per survey occasion) and as such we had a maximum of 3 paired survey histories, giving a maximum of
#' 6 sample occasions (see MacKenzie et al. 2006, MacKenzie and Royle 2005, and Laake et al. 2011 for various sample survey design details).  At each point, we surveyed for 5 minutes (technically broken into
#' 3 time intervals of 2, 2, and 1 minutes; not used here) and we noted detections by each observer and collected distance to each observation within a set of distance bins (0-25, 25-50, 50-75, 75-100m)
#' of the target species (Black-capped vireo's in this case) for each surveyor.  Our primary focus was to use mark-recapture distance sampling methods to estimate density of
#' Black-capped vireo's, and to estimate detection rates for the mark-recapture, distance, and composite model.
#'
#' @details  In addition to detailing the analysis used by Collier et al. (2013, In Review), this example documents the use of \code{mrds} for avian point count surveys and shows how density models
#' can be incorporated with occupancy models to develop spatially explicit density surface maps. For those that are interested, for the distance sampling portion of our analysis, we
#' used both conventional distance sampling (\code{cds}) and multiple covariate distance sampling (\code{mcds}) with uniform and half-normal key functions.  For the mark-recapture portion of
#' our analysis, we tended to use covariates for distance (median bin width), observer, and date of survey (days since 15 March 2011).
#'
#' We combined our \code{mrds} density estimates via a Horvitz-Thompson styled estimator with the resource selection function gradient developed
#' in Farrell et al. (2013) and estimated density on an ~3.14ha hexagonal grid across our study area, which provided a density gradient for the Fort Hood military installation.  Because there
#' was considerable data manipulation needed for each analysis to structure the data appropriately for use in \code{mrds}, rather than wrap each analysis in a single function, we have provided
#' both the Golden-cheeked warbler and Black-capped vireo analyses in their full detail.  The primary differences you will see will be changes to
#' model structures and model outputs between the two species.
#'
#' @name lfbcvi
#' @docType data
#' @format The format is a data frame with the following covariate metrics.
#' \describe{\item{PointID}{Unique identifier for each sample location; locations are the same for both species}
#' \item{VisitNumber}{Visit number to the point}
#' \item{Species}{Species designation, either Golden-cheeked warbler (GW) or Black-capped Vireo (BV)}
#' \item{Distance}{Distance measure, which is either NA (representing no detection), or the median of the binned detection distances}
#' \item{PairNumber}{ID value indicating which observers were paired for that sampling occasion}
#' \item{Observer}{Observer ID, either primary(1), or secondary (2)}
#' \item{Detected}{Detection of a bird, either 1 = detected, or 0 = not detected}
#' \item{Date}{Date of survey since 15 march 2011}
#' \item{Pred}{Predicted occupancy value for that survey hexagon based on Farrell et al. (2013)}
#' \item{Category}{Region.Label categorization, see \code{mrds} help file for details on data structure}
#' \item{Effort}{Amount of survey effort at the point}
#' \item{Day}{Number of days since 15 March 2011}
#' \item{ObjectID}{Unique ID for each paired observations} }
#'
#' @references Farrell, S.F., B.A. Collier, K.L. Skow, A.M. Long, A.J. Campomizzi, M.L. Morrison, B. Hays, and R.N. Wilkins. 2013. Using LiDAR-derived structural vegetation characteristics to develop high-resolution, small-scale, species distribution models for conservation planning. Ecosphere 43(3): 42. http://dx.doi.org/10.1890/ES12-000352.1
#' @references Laake, J.L., B.A. Collier, M.L. Morrison, and R.N. Wilkins. 2011. Point-based mark recapture distance sampling. Journal of Agricultural, Biological and Environmental Statistics 16: 389-408.
#' @references Collier, B.A., S.L. Farrell, K.L. Skow, A. M. Long, A.J. Campomizzi, K.B. Hays, J.L. Laake, M.L. Morrison, and R.N. Wilkins. 2013. Spatially explicit density of
#' endangered avian species in a disturbed landscape.  Auk, In Review.
#' @keywords datasets
#' @author Bret Collier and Jeff Laake
#' @examples
#'
#' \dontrun{
#' data(lfbcvi)
#' xy=cut(lfbcvi$Pred, c(-0.0001, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1),
#'   labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
#' x=data.frame(lfbcvi, New=xy)
#'
#' # Note that I scaled the individual covariate of day-helps with
#' # convergence issues
#' bird.data <- data.frame(object=x$ObjectID, observer=x$Observer,
#'                         detected=x$Detected, distance=x$Distance,
#'                         Region.Label=x$New, Sample.Label=x$PointID,
#'                         Day=(x$Day/max(x$Day)))
#'
#' # make observer a factor variable
#' bird.data$observer=factor(bird.data$observer)
#'
#' # Jeff Laake suggested this snippet to quickly create distance medians
#' # which adds bin information to the bird.data dataframe
#'
#' bird.data$distbegin=0
#' bird.data$distend=100
#' bird.data$distend[bird.data$distance==12.5]=25
#' bird.data$distbegin[bird.data$distance==37.5]=25
#' bird.data$distend[bird.data$distance==37.5]=50
#' bird.data$distbegin[bird.data$distance==62.5]=50
#' bird.data$distend[bird.data$distance==62.5]=75
#' bird.data$distbegin[bird.data$distance==87.5]=75
#' bird.data$distend[bird.data$distance==87.5]=100
#'
#' # Removed all survey points with distance=NA for a survey event;
#' # hence no observations for use in ddf() but needed later
#' bird.data=bird.data[complete.cases(bird.data),]
#'
#' # Manipulations on full dataset for various data.frame creation for
#' # use in density estimation using dht()
#'
#' #Samples dataframe
#' xx=x
#' x=data.frame(PointID=x$PointID, Species=x$Species,
#'              Category=x$New, Effort=x$Effort)
#' x=x[!duplicated(x$PointID),]
#' point.num=table(x$Category)
#' samples=data.frame(PointID=x$PointID, Region.Label=x$Category,
#'                    Effort=x$Effort)
#' final.samples=data.frame(Sample.Label=samples$PointID,
#'                          Region.Label=samples$Region.Label,
#'                          Effort=samples$Effort)
#'
#' #obs dataframe
#' obs=data.frame(ObjectID=xx$ObjectID, PointID=xx$PointID)
#' #used to get Region and Sample assigned to ObjectID
#' obs=merge(obs, samples, by=c("PointID", "PointID"))
#' obs=obs[!duplicated(obs$ObjectID),]
#' obs=data.frame(object=obs$ObjectID, Region.Label=obs$Region.Label,
#'                Sample.Label=obs$PointID)
#'
#' region.data=data.frame(Region.Label=c(1, 2, 3,4,5,6,7,8,9, 10),
#' Area=c(point.num[1]*3.14, point.num[2]*3.14,
#'        point.num[3]*3.14, point.num[4]*3.14,
#'        point.num[5]*3.14, point.num[6]*3.14,
#'        point.num[7]*3.14, point.num[8]*3.14,
#'        point.num[9]*3.14, point.num[10]*3.14))
#'
#' # Candidate Models
#'
#' BV1=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE,point=TRUE,width=100,breaks=c(0,50,100)))
#' BV1FI=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance),
#'    data=bird.data,
#'    method="io.fi",
#'    meta.data=list(binned=TRUE,point=TRUE,width=100,breaks=c(0,50,100)))
#' BV2=ddf(
#'    dsmodel=~mcds(key="hr",formula=~1),
#'    mrmodel=~glm(~distance),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE,point=TRUE,width=100,breaks=c(0,50,100)))
#' BV3=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance+observer),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' BV3FI=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance+observer),
#'    data=bird.data,
#'    method="io.fi",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' BV4=ddf(
#'    dsmodel=~mcds(key="hr",formula=~1),
#'    mrmodel=~glm(~distance+observer),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' BV5=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance*observer),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' BV5FI=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance*observer),
#'    data=bird.data,
#'    method="io.fi",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' BV6=ddf(
#'    dsmodel=~mcds(key="hr",formula=~1),
#'    mrmodel=~glm(~distance*observer),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' BV7=ddf(
#'    dsmodel=~cds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance*Day),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' BV7FI=ddf(
#'    dsmodel=~cds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance*Day),
#'    data=bird.data,
#'    method="io.fi",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' BV8=ddf(
#'    dsmodel=~cds(key="hr",formula=~1),
#'    mrmodel=~glm(~distance*Day),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' BV9=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance*observer*Day),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' BV9FI=ddf(
#'    dsmodel=~mcds(key="hn",formula=~1),
#'    mrmodel=~glm(~distance*observer*Day),
#'    data=bird.data,
#'    method="io.fi",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' BV10=ddf(
#'    dsmodel=~mcds(key="hr",formula=~1),
#'    mrmodel=~glm(~distance*observer*Day),
#'    data=bird.data,
#'    method="io",
#'    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#' #BV.DS=ddf(
#' #    dsmodel=~mcds(key="hn",formula=~1),
#' #    data=bird.data,
#' #    method="ds",
#' #    meta.data=list(binned=TRUE, point=TRUE, width=100,breaks=c(0,50,100)))
#'
#' #AIC table building code.
#' AIC = c(BV1$criterion, BV1FI$criterion, BV2$criterion, BV3$criterion,
#'         BV3FI$criterion, BV4$criterion,  BV5$criterion, BV5FI$criterion,
#'         BV6$criterion, BV7$criterion, BV7FI$criterion, BV8$criterion,
#'         BV9$criterion, BV9FI$criterion, BV10$criterion)
#'
#' #creates a set of row names for me to check my grep() call below
#' rn = c("BV1", "BV1FI", "BV2", "BV3",  "BV3FI", "BV4", "BV5", "BV5FI",
#'        "BV6", "BV7", "BV7FI", "BV8", "BV9", "BV9FI", "BV10")
#'
#' #Number parameters
#' k = c(length(BV1$par), length(BV1FI$par), length(BV2$par),
#'       length(BV3$par), length(BV3FI$par), length(BV4$par),
#'       length(BV5$par),length(BV5FI$par), length(BV6$par),
#'       length(BV7$par), length(BV7FI$par), length(BV8$par),
#        length(BV9$par), length(BV9FI$par), length(BV10$par))
#'
#' #build AIC table
#' AIC.table=data.frame(AIC = AIC, rn=rn, k=k, dAIC = abs(min(AIC)-AIC) ,
#'                      likg=exp(-.5*(abs(min(AIC)-AIC))))
#' #row.names(AIC.table)=grep("BV", ls(), value=TRUE)
#' AIC.table=AIC.table[with(AIC.table, order(-likg, -dAIC, AIC, k)),]
#' AIC.table=data.frame(AIC.table, wi=AIC.table$likg/sum(AIC.table$likg))
#' AIC.table
#'
#' # Model average N_hat_covered estimates
#' #  not very clean, but I wanted to show full process, need to use
#' #  collect.models and model.table here later on
#' estimate <- c(BV1$Nhat, BV1FI$Nhat, BV2$Nhat, BV3$Nhat, BV3FI$Nhat,
#'               BV4$Nhat,  BV5$Nhat, BV5FI$Nhat, BV6$Nhat, BV7$Nhat,
#'               BV7FI$Nhat, BV8$Nhat, BV9$Nhat, BV9FI$Nhat, BV10$Nhat)
#'
#' AIC.values=AIC
#'
#' # had to use str() to extract here as Nhat.se is calculated in
#' # mrds:::summary.io, not in ddf(), so it takes a bit
#' std.err <- c(summary(BV1)$Nhat.se, summary(BV1FI)$Nhat.se,
#'              summary(BV2)$Nhat.se, summary(BV3)$Nhat.se,
#'              summary(BV3FI)$Nhat.se, summary(BV4)$Nhat.se,
#'              summary(BV5)$Nhat.se, summary(BV5FI)$Nhat.se,
#'              summary(BV6)$Nhat.se, summary(BV7)$Nhat.se,
#'              summary(BV7FI)$Nhat.se,summary(BV8)$Nhat.se,
#'              summary(BV9)$Nhat.se, summary(BV9FI)$Nhat.se,
#'              summary(BV10)$Nhat.se)
#'}
#'
#' \dontrun{
#' #Not Run
#' #requires RMark
#' library(RMark)
#' #uses model.average structure to model average real abundance estimates for
#' #covered area of the surveys
#'   mmi.list=list(estimate=estimate, AIC=AIC.values, se=std.err)
#'   model.average(mmi.list, revised=TRUE)
#'
#' #Not Run
#' #Summary for the top 2 models
#'  #summary(BV5, se=TRUE)
#'  #summary(BV5FI, se=TRUE)
#'
#' #Not Run
#' #Best Model
#'  #best.model=AIC.table[1,]
#'
#' #Not Run
#' #GOF for models
#' #ddf.gof(BV5, breaks=c(0, 25, 50, 75, 100))
#'
#' #Not Run
#' #Density estimation across occupancy categories
#' #out.BV=dht(BV5, region.data, final.samples, obs, se=TRUE,
#' #           options=list(convert.units=.01))
#'
#' #Plot--Not Run
#'
#' #Composite Detection Function
#' #plot(BV5, which=3, showpoints=FALSE, angle=0, density=0, col="black", lwd=3,
#' # main="Black-capped Vireo",xlab="Distance (m)", las=1, cex.axis=1.25,
#' # cex.lab=1.25)
#'
#' }
#'
NULL

#' Tips on optimisation issues in \code{mrds} models
#'
#' Occasionally when fitting an `mrds` model one can run into optimisation issues. In general such problems can be quite complex so these "quick fixes" may not work. If you come up against problems that are not fixed by these tips, or you feel the results are dubious please go ahead and contact the package authors.
#'
#'
#' @section Debug mode:
#' One can obtain debug output at each stage of the optimisation using the \code{showit} option. This is set via \code{control}, so adding \code{control=list(showit=3)} gives the highest level of debug output (setting \code{showit} to 1 or 2 gives less output).
#'
#'
#' @section Re-scaling covariates:
#' Sometimes convergence issues in covariate (MCDS) models are caused by values of the covariate being very large, so a rescaling of that covariate is then necessary. Simply scaling by the standard deviation of the covariate can help (e.g. \code{dat$size.scaled <- dat$scale/sd(dat$scale)} for a covariate \code{size}, then including \code{size.scaled} in the model instead of \code{size}).
#'
#' It is important to note that one needs to use the original covariate (size) when computing Horvitz-Thompson estimates of population size if the group size is used in that estimate. i.e. use the unscaled size in the numerator of the H-T estimator.
#'
#'
#' @section Initial values:
#' Initial (or starting) values can be set via the \code{initial} element of the \code{control} list. \code{initial} is a list itself with elements \code{scale}, \code{shape} and \code{adjustment}, corresponding to the associated parameters. If a model has covariates then the \code{scale} or \code{shape} elements will be vectors with parameter initial values in the same order as they are specific in the model formula (using \code{showit} is a good check they are in the correct order). Adjustment starting values are in order of the order of that term (cosine order 2 is before cosine order 3 terms).
#'
#' One way of obtaining starting values is to fit a simpler model first (say with fewer covariates or adjustments) and then use the starting values from this simpler model for the corresponding parameters.
#'
#' Another alternative to obtain starting values is to fit the model (or some submodel) using Distance for Windows. Note that Distance reports the scale parameter (or intercept in a covariate model) on the exponential sclale, so one must \code{log} this before supplying it to \code{ddf}.
#'
#'
#' @section Bounds:
#' One can change the upper and lower bounds for the parameters. These specify the largest and smallest values individual parameters can be. By placing these constraints on the parameters, it is possible to "temper" the optimisation problem, making fitting possible.
#'
#' Again, one uses the \code{control} list, the elements \code{upperbounds} and \code{lowerbounds}. In this case, each of \code{upperbounds} and \code{lowerbounds} are vectors, which one can think of as each of the vectors \code{scale}, \code{shape} and \code{adjustment} from the "Initial values" section above, concatenated in that order. If one does not occur (e.g. no shape parameter) then it is simple omitted from the vector.
#'
#' @name mrds-opt
#' @docType methods
#' @author David L. Miller <dave@@ninepointeightone.net>
NULL
