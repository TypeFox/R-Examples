#' @title Returns the biomass of a stone pine tree
#'
#' @description Returns total biomass of a stone pine tree (wood and leaves, dry state) in kg given the 
#' diameter at breast height, using an allometric equation
#'
#' @note Use this function at you own risk, it has been validated for trees (ie: >24 cm diameters). 
#' @note The allometric equation takes the form of a pure quadratic equation
#' @seealso \code{\link{pureQuadraticEquation}}
#' @references Cutini, A. and Hajny, M. and Gugliotta, O. and Manetti, M. and Amorini, E. 2009, Effetti della struttura del popolamento sui modelli di stima del volume e della biomassa epigea (Pineta di Castelfusano - Roma) \emph{Forest@@}, \bold{6}, 75--84 
#'   Tipo B
#' @export
#' @param x a data frame holding diameters of branches
#' @param diameter the name of the column holding diameter of the x data frame, diameters should be in cm 
#' @return the total biomass of the branch of a stone pine (in kg, dry state)
#' @family Biomass
allometryCutini2009 <- function(x, diameter) {
  a <- -198.236
  b <- 0.620
  pureQuadraticEquation(a, b, as.double(x[diameter]))
}

#' @title Returns the biomass of a maritime pine branch
#'
#' @description Returns the woody biomass of a maritime pine branch (dry state, no leaves!) in kg given the 
#' diameter, using an allometric equation
#'
#' @note The allometric equation has been validated for <10 cm diameter branches, extrapolation on larger branches my yield unreasonable results.
#' @note The allometric equation takes the form of a power equation
#' @seealso \code{\link{powerEquation}}
#' @references Port\'{e}, A. and Trichet, P. and Bert, D. and Loustau, D. 2002, Allometric relationships for branch and tree woody biomass of Maritime pine (\emph{Pinus pinaster} Ait.) \emph{Forest Ecology and Management}, \bold{158}, 71--83
#' @export
#' @param x a data frame holding diameters of branches
#' @param diameter the name of the column holding diameter of the x data frame, diameters should be in cm 
#' @return the woody biomass (dry state, no leaves!) of the branch of a maritime pine (in kg)
#' @family Biomass
allometryPorte2002 <- function(x, diameter) {
  a <- 21.228
  b <- 2.818
  powerEquation(a, b, as.double(x[diameter])) / 1000
}

#' @title Returns the fresh weight of a stone pine branch
#'
#' @description Returns the fresh biomass of a stone pine branch in kg given the 
#' diameter, using an allometric equation
#'
#' @note The allometric equation has been validated for 8-16 cm diameter branches. 
#' @note The allometric equation takes the form of a power equation
#' @seealso \code{\link{powerEquation}}
#' @references Data collected by A. Ascarelli, non linear regression by M. Bascietto
#' @export
#' @param x a data frame holding diameters of branches
#' @param diameter the name of the column holding diameter of the x data frame, diameters should be in cm 
#' @return the fresh biomass of the branch of a stone pine (in kg)
#' @family Biomass
allometryAsca2011 <- function(x, diameter) {
  a <- 0.7201
  b <- 1.8882
  powerEquation(a, b, as.double(x[diameter]))
}

#' @title Returns the fresh weight of a stone pine branch
#'
#' @description Returns the fresh biomass of a stone pine branch in kg given the 
#' diameter, using an allometric equation
#'
#' @note The allometric equation has been validated for 5-16 cm diameter branches. 
#' @note The allometric equation takes the form of a power equation. This equation 
#' yields more correct results than \code{\link{allometryAsca2011}} since it has been
#' built on a wider range of branch diameters and it superseeds it. 
#' @seealso \code{\link{powerEquation}}
#' @references Data collected by A. Ascarelli and integrated by small diameter branches by 
#' M. Bascietto and B. De Cinti, non linear regression by M. Bascietto
#' @export
#' @param x a data frame holding diameters of branches
#' @param diameter the name of the column holding diameter of the x data frame, diameters should be in cm 
#' @return the fresh biomass of the branch of a stone pine (in kg)
#' @family Biomass
allometryABDC <- function(x, diameter) {
  a <- 0.16843
  b <- 2.43523
  powerEquation(a, b, as.double(x[diameter]))
}

#' @title Returns the result of a pure quadratic equation
#'
#' @description Returns the result of the pure quadratic equation \eqn{Y = a + bX^2}
#' given \eqn{a}, \eqn{b} and \eqn{X}
#'
#' @param a the parameter \eqn{a} in the pure quadratic equation
#' @param b the parameter \eqn{a} in the pure quadratic equation
#' @param x the dependent variable
#' @export
#' @return the dependent variable (\eqn{Y})
#' @family Biomass
pureQuadraticEquation <- function(a, b, x) {
  a + b * x^2
}

#' @title Returns the result of an exponential equation
#'
#' @description Returns the result of the exponential equation \eqn{Y = a * X^b}
#' given \eqn{a}, \eqn{b} and \eqn{X}
#'
#' @param a the parameter \eqn{a} in the exponential equation
#' @param b the parameter \eqn{b} in the exponential equation
#' @param x the independent variable
#' @export
#' @return the dependent variable (\eqn{Y})
#' @family Biomass
powerEquation <- function(a, b, x) {
  a * x^b
}

#' @title Estimates the wood biomass of logs and truncated branches
#'
#' @description Estimates the wood biomass of logs and truncated branches by
#' computing their volume (using Smalian's formula) and converting it
#' to fresh weight using wood fresh density.
#'
#' Smalian's formula: \eqn{V=\frac{Sb+Sd}{2}l} where \eqn{V} is the log volume, 
#' \eqn{Sb} is the aerea of the basal (lower) section, \eqn{Sd} is the 
#' area of the higher section and \eqn{l} is the length of the log.
#'
#' @note Diameters used to compute section areas should be measured under the bark layer! When this is not the case (scarcely ever!) and diameters include bark thickness the log biomass is somewhat over-estimated!
#' @param x the data frame holding the measures needed to perform the estimation
#' @param lowerD The name of the data frame column holding diameter of the lower section in cm
#' @param higherD The name of the data frame column holding the diameter of the higher section (usually smaller!) in cm
#' @param logLength The name of the data frame column holding the length of the log or branch in m
#' @param density The name of the data frame column holding the fresh density of the wood, defined as \eqn{D=\frac{V_f}{W_f}} where \eqn{V_f} is wood volume measured in the field (i.e. satured with water) in \eqn{m^3} and \eqn{W_f} is wood fresh weight in kg. Fresh density is measured in \eqn{\frac{kg}{m^3}}
#' @family Biomass
#' @references la Marca, O. \emph{Elementi di dendrometria} 2004, Patron Editore (Bologna), p. 119
logBiomass <- function(x, lowerD, higherD, logLength, density) {
  lowerS   <- pi * (as.double(x[lowerD]) / 200)^2 
  higherS  <- pi * (as.double(x[higherD]) / 200)^2
  l        <- as.double(x[logLength])
  volume   <- (lowerS + higherS) / 2 * l
  volume * density
}

#' @title Computes masses of branches and logs
#'
#' @description Computes branches biomass using an allometric function provided in \code{object$allometryFUN} and logs weight using Smalian's formula.
#'
#' Branches are telled apart from logs in the raw data frame (\code{object$fieldData}) because their final diameter is 0 (ie they have a tip) whereas logs have a final diameter > 0.
#'
#' @param object an object of \code{treeData} class
#' @return an object of \code{treeData} class
#' @seealso \code{\link{logBiomass}}
#' @family Biomass
#' @export
treeBiomass <- function(object) {
  ## gets stem and cut branches (ie. diameter at tip > 0) biomass, by converting its fresh volume to dry weight
  object$fieldData$biomass[(object$fieldData$dTip > 0)] <- as.vector(
      apply(
        object$fieldData[(object$fieldData$dTip > 0),], 
        1, 
        logBiomass, 
        lowerD    = "dBase", 
        higherD   = "dTip", 
        logLength = "length",
        object$density
      )
    )
    
    ## gets branches (ie. diameter at tip = 0) total biomass (wood + leaves) 
  object$fieldData$biomass[(object$fieldData$dTip == 0)] <- as.vector(
      apply(
        object$fieldData[(object$fieldData$dTip == 0),], 
        1, 
        object$allometryFUN, 
        diameter = "dBase"
      )
    )
  return(object)
}


#' @title Returns the total biomass of the tree
#'
#' @description This is just a helper function, it sums the biomass of all logs and branches, as previously computed by \code{\link{treeBiomass}}
#'
#' @note This function may be used to compute the moment of the tree. Tree biomass (multiplied by standard gravity) is the tree force applied to its CM.
#'
#' @param treeData A named list that includes a \code{fieldData} data frame element, holding a \code{biomass}-named column. Note that the \code{biomass} column is added to the data frame by a previous call to \code{\link{treeBiomass}} function
#' @return a real number or \code{FALSE} if the \code{biomass} column is \code{NA}
#' @export
#' @family Biomass
#' @examples 
#' library(treecm)
#' data(stonePine1TreeData)
#' print(treeTotalBiomass(stonePine1TreeData))
treeTotalBiomass <- function(treeData) {
  if (any(grepl("biomass", colnames(treeData$fieldData)))) 
    return(sum(treeData$fieldData$biomass))
  else
    return(FALSE)
}
