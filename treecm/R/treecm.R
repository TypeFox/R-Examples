###############################################################################
#
# treecm: An R package to assess tree centre of mass
# author: Marco Bascietto <marco.bascietto@ibaf.cnr.it>
#
# This is released under a GPL license.
#
# Documentation was created using roxygen:
# package.skeleton(name="treecm", code_files="treecm.R")
# library(roxygen2); roxygenize('treecm', 'treecm', overwrite=TRUE, unlink.target=FALSE, copy.package=FALSE)
###############################################################################

#' @title Assessment of the position of the centre of mass of trees
#'
#' @description The centre of mass is a crucial data for arborists in order to consolidate a tree using steel or dynamic cables. Given field-recorded data on branchiness of a tree, the package:
#' \itemize{
#'   \item{computes and plots the centre of mass of the tree itself}
#'   \item{simulates the shift in CM position as branches are pruned}
#'   \item{computes branches slenderness coefficient in order to aid the arborist identify potentially dangerous branches}
#'   \item{computes the force acting on a ground plinth and its best position relating to the tree centre of mass, should the tree need to be stabilized by a steel cable}
#' }
#' The tree stem is ideally sectioned in logs. The weight of tree components is assessed based on
#' \itemize{
#'   \item the sum of volume of stem logs
#'   \item the sum of branches biomass
#' }
#' Field measures to be taken on logs and branches are described in \code{\link{importFieldData}} and are to be recorded on the tree itself, possibly using tree-climbing tecniques.
#' In order to help the arborist in the pruning selection process a simple plot of branch coefficient of slenderness is implemented.
#' @note \bold{Branch biomass} is computed by allometric equations relating its fresh weight (wood + leaves) to its diameter at point of insertion on the stem. \bold{Log biomass} is computed by converting its volume to weight using wood fresh density. Volume is computed using Smalian's formula (see \code{\link{logBiomass}} description).
#' A sample \code{.CSV} file is provided to guide through data filling in the field
#' @seealso \code{\link{logBiomass}}
#' \code{\link{importFieldData}}
#' @name treecm-package
#' @aliases treecm
#' @docType package
#' @import plyr
#' @author Marco Bascietto \email{marco.bascietto@@cnr.it}
#' @keywords package
#' @examples
#' data(stonePine1TreeData)
#' vectors  <- treeVectors(stonePine1TreeData)
#' CM       <- centreOfMass(vectors)
#' plot(vectors, 
#'    main = "Centre Of Mass", 
#'    col = "grey30", 
#'    txtcol = "grey30")
#' plot(CM)
#' summary(CM)
#' @references Source code is hosted at GitHub (\url{http://mbask.github.com/treecm/})
NULL




#' Green wood density data for a few tree species
#' 
#' Wood density is used to convert wood volume measures in the field to their
#' fresh weights. It is measured in \eqn{\frac{kg}{m^3}}.
#' 
#' Density is measured at humidity level 50%, very close to mean living tree
#' wood humidity. The dataset is provided as a reference only, please be
#' cautioned about using these values on your samples.
#' 
#' @name Dst
#' @docType data
#' @format A data frame with 170 observations on the following 3 variables. \code{
#' data.frame:  170 obs. of  3 variables:
#'   $ species: chr  "Abies alba" "Abies alba" "Abies balsama" "Abies grandis" ...
#'   $ group  : Factor w/ 2 levels "conifer","dicot": 1 1 1 1 1 1 1 1 1 1 ...
#'   $ density: int  545 577 529 449 465 673 689 497 673 577 ...}
#' @source Niklas, K. J. and Spatz, H.-C. (2010) Worldwide correlations of
#' mechanical properties and green wood density. American Journal of Botany,
#' 97, 1587-1594
#' @keywords datasets
#' @examples
#' 
#' data(Dst)
NULL 



#' Raw CSV file of field recorded values for a stone pine tree
#' 
#' Required data for the assessment of the centre of mass have been recorded in
#' the field for a stone pine (\emph{Pinus pinea} L.). This is an example of
#' csv file that should be fed to \code{\link{treeBiomass}} to assess tree
#' centre of mass.
#'
#' This dataset has been collected for a 17.1 metres tall stone pine whose stem was tilted approx. 20 degrees from the vertical plane (or 80 degrees from the horizontal plane). The stem has been sectioned in two logs (\code{L1} and \code{L2}), and a final branch (\code{C}).
#' 
#' The \code{.csv} file must contain all column headings listed in \code{\link{importFieldData}}, regardless of them being optional (no data in them) or mandatory.
#' 
#' @name stonePine1FieldData
#' @docType data
#' @format \code{
#' "code","azimuth","dBase","dTip","length","tipD","height","tilt","toBePruned","pathToTip"
#' "L1",275,73,41,"10.2","2.5",0,80,"FALSE","TRUE"
#' "L2",275,41,16,"3.9","2.75","10.2",80,"FALSE","TRUE"
#' "B1",190,15,0,"7.95","10.1",,,"FALSE","FALSE"
#' "B2",200,22,0,"7.95","10.4",,,"FALSE","FALSE"
#' "B3",230,15,0,"7.95","10.4",,,"FALSE","FALSE"
#' "B4",200,18,0,"7.95","11.15",,,"FALSE","FALSE"
#' "B5",180,7,0,"7.95","11.3",,,"FALSE","FALSE"
#' "B6",150,6,0,"7.95","11.3",,,"FALSE","FALSE"
#' "B7",340,16,0,"3.95","11.3",,,"FALSE","FALSE"
#' "B8",220,13,0,"7.95","11.8",,,"FALSE","FALSE"
#' "B9",165,19,0,"7.95","11.8",,,"FALSE","FALSE"
#' "B10",280,8,0,"3.95","11.9",,,"FALSE","FALSE"
#' "B11",170,9,0,"7.95","11.9",,,"FALSE","FALSE"
#' "B12",265,8,0,"7.95","12.2",,,"FALSE","FALSE"
#' "B13",75,6,0,"3.95","12.2",,,"FALSE","FALSE"
#' "B14",180,6,0,"7.95","12.2",,,"FALSE","FALSE"
#' "B15",170,6,0,"7.95","12.6",,,"FALSE","FALSE"
#' "B16",120,5,0,"7.95","12.6",,,"FALSE","FALSE"
#' "B17",10,14,0,"3.95",13,,,"FALSE","FALSE"
#' "B18",180,13,0,"7.95",13,,,"FALSE","FALSE"
#' "B19",260,13,0,"7.95","13.2",,,"FALSE","FALSE"
#' "B20",75,6,0,"3.95","13.2",,,"FALSE","FALSE"
#' "B21",75,10,0,"3.95","13.75",,,"FALSE","FALSE"
#' "B22",215,7,0,"7.95","13.75",,,"FALSE","FALSE"
#' "B23",140,7,0,"7.95","13.75",,,"FALSE","FALSE"
#' "C",275,16,0,3,3,"14.1",80,"FALSE","TRUE"
#' }
#' @source Original data collected by the author
#' @keywords datasets
#' @examples
#' library("treecm")
#' csvFileName <- system.file("data", "stonePine1FieldData.csv.gz", package = "treecm")
#' treeData <- importFieldData(
#'   csvFileName, 
#'   650, 
#'   allometryABDC
#' )
#' head(treeData$fieldData)
NULL


#' Field recorded values for a stone pine tree
#' 
#' Required data for the assessment of the centre of mass have been recorded in
#' the field for a stone pine (\emph{Pinus pinea} L.).
#' \code{\link{treeBiomass}} has already been run on the dataset, vectors have
#' yet to be computed.
#' 
#' This dataset includes a list of 4 elements:
#' \itemize{
#' \item{the \code{\link{stonePine1FieldData}} dataset}
#' \item{the density of wood}
#' \item{the allometry function to be used to compute branches biomass}
#' \item{the estimate of branches centre of mass}
#' }
#'
#' @name stonePine1TreeData
#' @docType data
#' @format The format is: 
#' \code{
#'  List of 4
#' $ fieldData   :'data.frame':  26 obs. of  10 variables:
#' ..$ azimuth   : int [1:26] 275 275 190 200 230 200 180 150 340 220 ...
#' ..$ dBase     : int [1:26] 73 41 15 22 15 18 7 6 16 13 ...
#' ..$ dTip      : num [1:26] 41 16 0 0 0 0 0 0 0 0 ...
#' ..$ length    : num [1:26] 10.2 3.9 NA NA NA NA NA NA NA NA ...
#' ..$ tipD      : num [1:26] 2.5 2.75 7.95 7.95 7.95 7.95 7.95 7.95 7.95 3.95 ...
#' ..$ height    : num [1:26] 0 10.2 10.1 10.4 10.4 ...
#' ..$ tilt      : num [1:26] 80 80 0 0 0 0 0 0 0 0 ...
#' ..$ toBePruned: logi [1:26] FALSE FALSE FALSE FALSE FALSE FALSE ...
#' ..$ pathToTip : logi [1:26] TRUE TRUE FALSE FALSE FALSE FALSE ...
#' ..$ biomass   : num [1:26] 1825 193 123 313 123 ...
#' $ density     : num 650
#' $ allometryFUN:function (x, diameter)  
#'   $ branchesCM  : num 1
#' }
#' @source Original data collected by the author
#' @keywords datasets
#' @examples
#' data(stonePine1TreeData)
#' vectors  <- treeVectors(stonePine1TreeData)
#' CM       <- centreOfMass(vectors)
#' summary(CM)
#' # The steps to recreate this dataset:
#' csvFileName <- system.file("data", "stonePine1FieldData.csv.gz", package = "treecm")
#' treeData <- importFieldData(csvFileName, 650, allometryABDC)
#' treeData <- treeBiomass(treeData)
NULL



#' Raw CSV file of field recorded values for a stone pine tree
#' 
#' Required data for the assessment of the centre of mass have been recorded in
#' the field for a stone pine (\emph{Pinus pinea} L.). This is an example of
#' csv file that should be fed to \code{\link{treeBiomass}} to assess tree
#' centre of mass.
#' 
#' This dataset has been collected for a \eqn{\approx 11} metres tall stone pine with a small number of very large branches.
#'
#' The \code{.csv} file must contain all column headings listed in \code{\link{importFieldData}}, regardless of them being optional (no data in them) or mandatory.
#' 
#' @name stonePine2FieldData
#' @docType data
#' @format \code{
#' "code","azimuth","dBase","dTip","length","tipD","height","tilt","toBePruned","pathToTip"
#' "L1",0,67,40,"6.8",0,0,90,,"TRUE"
#' "B1",250,40,,,"7.8","6.8",,,
#' "B2",240,32,,,"8.9","7.8",,,
#' "B3",55,25,,,9,9,,,
#' "B4",260,10,,,"5.5",10,,"TRUE",
#' "B5",80,36,,,"8.2","10.6",,,
#' "B6",255,27,,,"4.5","10.8",,,
#' "B7",0,40,,,"6.5","6.8",70,,"TRUE"
#' }
#' @source Original data collected by the author
#' @keywords datasets
#' @examples
#' library("treecm")
#' csvFileName <- system.file("data", "stonePine2FieldData.csv.gz", package = "treecm")
#' treeData <- importFieldData(csvFileName, 650, allometryABDC)
#' head(treeData$fieldData)
NULL



#' @title Plots a segment
#'
#' @description Plots a segmente given two set of polar coordinates (angle, distance
#' from tree base), may be used to represent buildings close to the tree on the CM plot
#' 
#' @param a0 angle of first set of coordinates
#' @param d0 distance of first set of coordinates
#' @param a1 angle of second set of coordinates
#' @param d1 distance of second set of coordinates
#' @return \code{NULL}
#' @export
#' @importFrom graphics points
#' @importFrom graphics segments
plotPolarSegment <- function(a0, d0, a1, d1) {
  xy0 <- toCartesianXY(a0, d0)
  xy1 <- toCartesianXY(a1, d1)
  points(xy0[1], xy0[2])
  points(xy1[1], xy1[2])
  segments(xy0[1], xy0[2], xy1[1], xy1[2])
}

#' @title Imports field data from csv file
#'
#' @description Imports \code{csv} file holding field recorded data returning a list holding field and other key data provided as arguments. 
#'
#' Field measures to be taken on logs and branches include:
#' \itemize{
#'   \item{\bold{Azimuth} (\code{azimuth}): mean angle of orientation of the branch or log measured from the base of the tree (usually with magnetic north as reference, measured clockwise), \emph{mandatory}}
#'   \item{\bold{Diameter at base} (\code{dBase}): diameter at insertion point on the stem, for branches, diameter of the lower section for logs, \emph{mandatory}}
#'   \item{\bold{Diameter at top} (\code{dTip}): always 0 for branches, diameter of the higher section for logs, defaults to 0, \emph{mandatory} only for logs}
#'   \item{\bold{Length} (\code{length}): length of logs or branches, \emph{mandatory} for logs, \emph{mandatory} for branches only if their slenderness coefficient is to be computed}
#'   \item{\bold{Distance} (\code{tipD}): length of branch or log projection on the ground, starting from tree base to tip of branch or log, \emph{mandatory}}
#'   \item{\bold{Height} (\code{height}): height of branch insertion on the stem or height of lower section of the log, to be used to compute \code{z} coordinate of CM, defaults to NA, \emph{mandatory} only for \code{z} determination of CM}
#'   \item{\bold{Tilt} (\code{tilt}): mean branch tilt or log tilt from the horizontal plane (eg a vertical branch is 90 degrees, an horizontal branch is 0 degrees), to be used to compute \code{z} coordinate of CM, defaults to 0, \emph{mandatory} only for \code{z} determination of CM. Note however that the tree tip should be considered as a branch, not a log, in order to account for foliage biomass. In this case the tilt value should be recorded otherwise it would default to 0, \emph{ie} an horizontal branch}
#'   \item{\bold{To be pruned} (\code{toBePruned}): a boolean value, defaults to FALSE, \emph{optional}}
#'   \item{\bold{Path to tip} (\code{pathToTip}): a boolean value TRUE on logs and branches that make part of the ``main stem'' of the tree, defaults to FALSE, \emph{optional}}
#' }
#' 
#' @param fileName Name of csv file holding field data
#' @param dst Fresh density of wood of the tree
#' @param branchesAllometryFUN the function that should compute branch biomass from its diameter
#' @param bCM Estimated position of the centre of mass of branches, a real number from 0.01 (CM at branch base) to 1.00 (CM at branch tip). As a rule of thumb, average live branches, with an average amount of foliage, have CM approx. from \eqn{1/3} to \eqn{2/3} of their length. bCM = 1.0 (default value)
#' @seealso \code{\link{getCoordinatesAndMoment}}
#' @return a list of 4 elements: field data, wood fresh density, allometryFUN function and branches CM
#' @export
#' @importFrom utils read.csv
importFieldData <- function(fileName, dst, branchesAllometryFUN, bCM = 1) {
  tree <- read.csv(fileName, row.names = 1)
  
  tree <- within(tree, {
    dTip[is.na(dTip)] <- 0
    tilt[is.na(tilt)] <- 0
    toBePruned[is.na(toBePruned)] <- FALSE
    pathToTip[is.na(pathToTip)] <- FALSE
  })
  
  list(
    fieldData    = tree, 
    density      = dst, 
    allometryFUN = branchesAllometryFUN,
    branchesCM   = bCM
  )
  
}

#' @title Stores branches CM in an object
#'
#' @description Stores branches CM in the object provided as an argument. branchesCM has to be in the range 0.01 - 1
#'
#' @note Method \code{\link{treeVectors}} must be invoked to take changes into effect
#'
#' @param object the object of class \code{treeData}
#' @param value the new branch CM
#' @return the updated list
#' @export
setBranchesCM <- function(object, value) {
  if (is.list(object) && value > 0 && value <=1) 
    object$branchesCM <- value
  return(object)
}
