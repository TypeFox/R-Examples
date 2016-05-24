#' @title Computes the modulus of the CM vector
#'
#' @description The Centre of Mass vector starts from the tree base and points towards the CM.
#' Its modulus is the distance between the CM (x, y, z) and tree base (0, 0, 0).
#'
#' @note This function is mainly needed to compute the moment of the tree. The CM modulus is the tree moment arm.
#'
#' @param object an object of \code{CM} class
#' @return a real number
#' @export
#' @family Stabilization
#' @examples 
#' library(treecm)
#' data(stonePine1TreeData)
#' vectors <- treeVectors(stonePine1TreeData)
#' CM <- centreOfMass(object=vectors)
#' print(centreOfMassModulus(CM))
#' treeMoment <- buildTreeMomentObject(
#'  centreOfMassModulus(CM)
#'  , treeTotalBiomass(stonePine1TreeData)
#'  , centreOfMassAngle(CM)
#'  )
centreOfMassModulus <- function(object) {
  with(object, sqrt(x^2 + y^2 + z^2))
}

#' @title Returns the angle between CM modulus and the tree weight vector
#'
#' @description This function is mainly needed to compute the moment of the tree. The angle is need to compute the projection of the tree weight normal to the CM modulus
#'
#' @param object an object of \code{CM} class
#' @return a real number in radians
#' @export
#' @family Stabilization
#' @examples 
#' library(treecm)
#' data(stonePine1TreeData)
#' vectors <- treeVectors(stonePine1TreeData)
#' CM <- centreOfMass(vectors)
#' print(centreOfMassAngle(CM))
#' treeMoment <- buildTreeMomentObject(
#'  centreOfMassModulus(CM)
#'  , treeTotalBiomass(stonePine1TreeData)
#'  , centreOfMassAngle(CM)
#'  )
centreOfMassAngle <- function(object) {
  with(object, atan(sqrt(x^2 + y^2) / z))
}

#' @title Returns tree logs and branches being part of the ``main stem''
#'
#' @description The ``main stem'' is not very clearly defined. Most softwood species (\emph{ie} \emph{Picea}, \emph{Abies}, \emph{Pseudotsuga} etc.) have only one stem given by their dominant apex. Hardwood species (\emph{Quercus}, \emph{Tilia} etc) and some softwood species (\emph{eg} \emph{Pinus pinea}) do not exhibit a dominant apex and branches often enlarge and grow taller than the apex. In these latter cases one has to select an appropriate path from tree base to tree tip, according to what may be considered the ``main stem'' of the tree. Both in former and latter cases the path selection has to be laid down in the \code{fieldData} data frame.
#'
#' @note Selected branches and logs have a \code{TRUE} value in the \code{pathToTip} column. This is a necessary step
#' towards anchor force detemrination, as the returned data frame has to be submitted to \code{\link{getPlinthForce}}
#'
#' @param treeData A named list that includes a \code{fieldData} data frame element, holding \code{pathToTip}-, \code{azimuth}-, \code{length}-, \code{tilt}-named columns
#' @return a data frame subsetted from the \code{fieldData} data frame having \code{TRUE} selected branches and logs, with three columns: \code{azimuth}, \code{length}, \code{tilt}. The first row if filled with zeros.
#' @export
#' @family Stabilization
#' @seealso \code{\link{importFieldData}}
#' @examples 
#' library(treecm)
#' data(stonePine1TreeData)
#' logs <- logPathSelection(stonePine1TreeData)
logPathSelection <- function(treeData) {
  rbind(0, with(treeData, fieldData[
    fieldData$pathToTip, c("azimuth", "length", "tilt")
    ]
  ))
}

#' @title Returns vector cartesian coordinates
#'
#' @description Given a modulus, a tilt angle and an azimuth angle it returns the vector cartesian coordinates
#'
#' @param x a named vector of three elements (z, x, y)
#' @return a list holding z, x, y coordinates
#' @export
#' @family Stabilization
toCartesianXYZ <- function(x) {
  mat <- matrix(
    c(
      sin(x["tiltRad"]), rep(cos(x["tiltRad"]), 2),
      1, sin(x["azRad"]), cos(x["azRad"]),
      rep(x["length"], 3)
      ),
    nrow = 3, ncol = 3
    )
  apply(mat, 1, prod)
}

#' @title Min/max values for the anchor position along the stem
#'
#' @description The anchor level must not be lower than the tree CM (for obvious static reasons) and higher than ``main stem'' height
#'
#' @param logs a data frame holding the selected logs (see \code{\link{logPathSelection}})
#' @param CM an object of \code{CM} class
#' @return a named vector of 2 elements: \item{z}{the height of the CM} \item{hMax}{the height of the ``main stem'' tip}
#' @export
#' @family Stabilization
#' @examples 
#' library(treecm)
#' data(stonePine1TreeData)
#' vectors <- treeVectors(stonePine1TreeData)
#' CM <- centreOfMass(vectors)
#' logs <- logPathSelection(stonePine1TreeData)
#' anchorRange(logs, CM)
anchorRange <- function(logs, CM) {
  c(CM["z"], hMax = cumsum(logs$length * sin(logs$tilt*pi/180))[nrow(logs)])
}

#' @title Computes the force on the plinth on the ground
#'
#' @description To stabilize the tree a steel cable is connected from an anchor point on the tree to a plinth on the ground.
#' The function computes the force on the plinth (needed to choose the appropriate steel cable and to build the plinth itself) and the maximum security azimuth (the angle relative to the North from the tree base).
#' Force is computed by comparing the moment of the tree and the moment of the anchor, whose arm is the vector from tree base to the anchor point.
#' The anchor point is defined as the distance from tree base, along the stem. Note that this distance equals anchor height only when the stem is perfectly vertical and straight.
#'
#' @note The function is vectorized both for anchor distance from tree base (\code{l.stem} parameter) and for cable length (\code{d}). It is not possible to pass invalid \code{l.stem} values, see \code{\link{anchorRange}}.
#'
#' @param l.stem the distance from tree base to anchor point, along the stem
#' @param d the length of the cable (in metres)
#' @param logs a data frame holding the selected logs (see \code{\link{logPathSelection}})
#' @param treeMoment the moment of the tree as computed by \code{\link{calcMoment}}
#' @param CM an object of \code{CM} class
#' @return a named list of 6 elements: \item{force}{the force in Newton on the plinth} \item{distanceOnGround}{the distance from tree base to the plinth} \item{anchorAlongStem}{the distance from tree base to the anchor (\emph{ie} \code{l.stem})} \item{cableLength}{the length of the cable (\emph{ie} \code{d})} \item{anchorHeight}{true height of the anchor over ground} \item{azimuth}{the azimuth of the plint relative to the tree base}
#' @export
#' @family Stabilization
#' @examples 
#' library(treecm)
#' data(stonePine1TreeData)
#' vectors <- treeVectors(stonePine1TreeData)
#' CM <- centreOfMass(vectors)
#' treeMoment <- buildTreeMomentObject(
#'   centreOfMassModulus(CM)
#'   , treeTotalBiomass(stonePine1TreeData)
#'   , centreOfMassAngle(CM)
#'   )
#' treeMoment <- calcMoment(treeMoment)
#' logs <- logPathSelection(stonePine1TreeData)
#' plinth <- data.frame(getPlinthForce(10, 20, logs, getMoment(treeMoment), CM))
getPlinthForce <- function(l.stem, d, logs, treeMoment, CM) {
  ## Controllo congruenza dei dati definiti dall'utente
  ## Ferma se: l'altezza dell'ancoraggio sul fusto e' sotto il baricentro o se e' maggiore dell'altezza del fusto
  aR <- anchorRange(logs, CM)
  stopifnot(l.stem >= aR[["z"]], l.stem <= aR[["hMax"]])
  rm(aR)
  
  ## trasformo gli angoli in radianti
  logs$tiltRad <- logs$tilt * pi / 180
  logs$azRad   <- logs$azimuth * pi / 180
  
  ## costruisce una matrice con le lunghezze cumulate dei toppi e le coordinate z,x,y 
  ## delle cime dei toppi
  anchorData <- as.data.frame(cbind(l.sum = cumsum(logs$length), t(apply(logs, 1, toCartesianXYZ))))
  colnames(anchorData)[2:4] <- c("l.z", "l.x", "l.y")

  ## calcolo le risultanti del vettore piede-toppi
  anchorData$lz.sum <- cumsum(anchorData$l.z)
  anchorData$lx.sum <- cumsum(anchorData$l.x)
  anchorData$ly.sum <- cumsum(anchorData$l.y)

  f.anchor <- vector(); d.x <- vector(); lz.sumV <- vector(); 
  for (i in seq_along(l.stem)) {
    ## numero del toppo sottostante a quello ove viene posto l'ancoraggio (paragrafo 3.6)
    logIndex <- which.max(anchorData$l.sum >= l.stem[i]) - 1
    
    ## modulo del vettore dalla cima del toppo precedente al punto di ancoraggio (par. 3.6, punto 2)
    anchor <- l.stem  - anchorData$l.sum[logIndex]
  
    ## coordinate del punto di applicazione del vettore anchor (par. 3.6 punti 3, 4, 5)
    priorLog <- unlist(anchorData[logIndex, c("lz.sum", "lx.sum", "ly.sum")])
    
    ## Passiamo al toppo successivo
    logIndex <- logIndex + 1
    
    ## Coordinate del punto di ancoraggio (par. 3.6 punto 6)
    anchorPoint <- with(logs, toCartesianXYZ(
      c(
        tiltRad = tiltRad[logIndex]
        , azRad = azRad[logIndex]
        , length = anchor[i]
      ))
    + priorLog)
    
    ## Modulo del vettore lanchor dal piede dell'albero al punto di ancoraggio (par. 3.7)
    anchorPoint[["l.anchor"]] <- sqrt(sum(anchorPoint^2))
    
    ## Angolo beta (par. 3.8)
    betaAngle <- acos(anchorPoint[["lz.sum"]] / d) - acos(anchorPoint[["lz.sum"]] / anchorPoint[["l.anchor"]])
  
    ## Calcolo di f.anchor (par. 4.1)
    f.anchor <- append(f.anchor, treeMoment / (anchorPoint[["l.anchor"]] * sin(betaAngle))) # [N]
    
    ## Distanza orizzontale tra piede dell'albero e plinto
    ## assumendo che il plinto venga messo in asse con la proiezione
    ## a terra dell'ancoraggio
    d.x <- append(d.x, sqrt(d^2 - anchorPoint[["lz.sum"]]^2) - sqrt(anchorPoint[["lx.sum"]]^2 + anchorPoint[["ly.sum"]]^2))
    
    lz.sumV <- append(lz.sumV, anchorPoint[["lz.sum"]])
  }
  ## Angolo dal piede dell'albero (opposto a quello baricentro-piede)
  polarAngle <- (toPolar(CM[["x"]], CM[["y"]])[1] + 180) %% 360

  return(list(
    "force" = f.anchor
    , "distanceOnGround" = d.x
    , "anchorAlongStem" = rep(l.stem, each = length(d))
    , "cableLength" = rep(d, length(l.stem))
    , "anchorHeight" = rep(lz.sumV, each = length(d))
    , "azimuth" = polarAngle
  ))
}

#' @title Constructor for the generic class moment
#'
#' @description Create an instance of a moment class holding mass, arm length and angle between the arm and the weight vector
#'
#' @param length the length of the arm of the moment
#' @param mass the mass of the moment (in kg)
#' @param angle the angle between the moment arm and the vector pointing towards the ground (the weight vector)
#' @return an instance of \code{moment} and list classes
#' @export
#' @family Stabilization momentClass
buildMomentObject <- function(length, mass, angle) {
  moment <- list(length = length, mass = mass, angle = angle)
  class(moment) <- c("moment", class(moment))
  return(moment)
}

#' @title Computes moment and returns the moment object
#'
#' @description Moment is computed as \eqn{M=l \cdot F}, where l is moment arm, F is the component of the force (\emph{mass times g}) normal to moment arm
#'
#' @param object an instance of moment class
#' @param g the standard gravity
#' @return the updated moment object
#' @export
#' @family Stabilization momentClass
#' @seealso \code{\link{getPlinthForce}}
calcMoment <- function(object, g = 9.81) {
  object$moment <- with(object, length * mass * sin(angle) * g)
  return(object)
}

#' @title Returns the moment
#'
#' @description For a description of how moment is computes see \code{\link{calcMoment}}. The function computes moment, if not yet done in a previous user call to \code{\link{calcMoment}}, and only returns it without updating the moment object
#'
#' @param object an instance of moment class
#' @return the moment figure
#' @export
#' @family Stabilization momentClass
#' @seealso \code{\link{getPlinthForce}}
getMoment <- function(object) {
  if (is.null(object$moment))
    object <- calcMoment(object)
  return(object$moment)
}

#' @title Constructor for the generic class treeMoment
#'
#' @description The class inherits from \code{moment} without adding any more properties
#'
#' @param length the distance from tree base to CM, as computed by \code{\link{centreOfMassModulus}}
#' @param mass the tree mass (in kg), not its weight, as computed by \code{\link{treeTotalBiomass}}
#' @param angle the angle between the moment arm (from tree base to CM) and the vector pointing towards the ground (the tree weight vector), as computed by \code{\link{centreOfMassAngle}}
#' @return an instance of \code{treeMoment}, \code{moment} and list classes
#' @export
#' @family Stabilization momentClass
#' @seealso \code{\link{calcMoment}}
buildTreeMomentObject <- function(length, mass, angle) {
  treeMoment <- buildMomentObject(length, mass, angle)
  class(treeMoment) <- c("treeMoment", class(treeMoment))
  return(treeMoment)
}
