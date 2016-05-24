#' Simulate components of continuous variables of population data
#'
#' Simulate components of continuous variables of population data by resampling
#' fractions from survey data. The continuous variable to be split and any
#' categorical conditioning variables need to be simulated beforehand.
#'
#' @name simComponents
#' @param simPopObj a \code{\linkS4class{simPopObj}}-object.
#' @param total a character string specifying the continuous variable of dataP
#' that should be split into components. Currently, only one variable can be
#' split at a time.
#' @param components a character vector specifying the components in
#' \code{dataS} that should be simulated for the population data.
#' @param conditional an optional character vector specifying categorical
#' conditioning variables for resampling. The fractions occurring in
#' \code{dataS} are then drawn from the respective subsets defined by these
#' variables.
#' @param replaceEmpty a character string; if \code{conditional} specifies at
#' least two conditioning variables, this determines how replacement cells for
#' empty subsets in the sample are obtained. If \code{"sequential"}, the
#' conditioning variables are browsed sequentially such that replacement cells
#' have the same value in one conditioning variable and minimum Manhattan
#' distance in the other conditioning variables. If no such cells exist,
#' replacement cells with minimum overall Manhattan distance are selected. The
#' latter is always done if this is \code{"min"} or only one conditioning
#' variable is used.
#' @param seed optional; an integer value to be used as the seed of the random
#' number generator, or an integer vector containing the state of the random
#' number generator to be restored.
#' @return An object of class \code{\linkS4class{simPopObj}} containing survey
#' data as well as the simulated population data including the components of
#' the continuous variable specified by \code{total} and \code{components}.
#' @note The basic household structure, any categorical conditioning variables
#' and the continuous variable to be split need to be simulated beforehand with
#' the functions \code{\link{simStructure}}, \code{\link{simCategorical}} and
#' \code{\link{simContinuous}}.
#' @author Stefan Kraft and Andreas Alfons and Bernhard Meindl
#' @export
#' @seealso \code{\link{simStructure}}, \code{\link{simCategorical}},
#' \code{\link{simContinuous}}, \code{\link{simEUSILC}}
#' @keywords datagen
#' @examples
#' data(eusilcS)
#' \dontrun{
#' ## approx. 20 seconds computation time
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize",
#'   strata="db040", weight="db090")
#' simPopObj <- simStructure(data=inp, method="direct",
#'   basicHHvars=c("age", "rb090", "hsize", "pl030", "pb220a"))
#' simPopObj <- simContinuous(simPopObj, additional = "netIncome",
#'   regModel = ~rb090+hsize+pl030+pb220a+hsize,
#'   method="multinom", upper=200000, equidist=FALSE, nr_cpus=1)
#'
#' # categorize net income for use as conditioning variable
#' sIncome <- manageSimPopObj(simPopObj, var="netIncome", sample=TRUE, set=FALSE)
#' sWeight <- manageSimPopObj(simPopObj, var="rb050", sample=TRUE, set=FALSE)
#' pIncome <- manageSimPopObj(simPopObj, var="netIncome", sample=FALSE, set=FALSE)
#'
#' breaks <- getBreaks(x=unlist(sIncome), w=unlist(sWeight), upper=Inf, equidist=FALSE)
#' simPopObj <- manageSimPopObj(simPopObj, var="netIncomeCat", sample=TRUE,
#'   set=TRUE, values=getCat(x=unlist(sIncome), breaks))
#' simPopObj <- manageSimPopObj(simPopObj, var="netIncomeCat", sample=FALSE,
#'   set=TRUE, values=getCat(x=unlist(pIncome), breaks))
#'
#' # simulate net income components
#' simPopObj <- simComponents(simPopObj=simPopObj, total="netIncome",
#'   components=c("py010n","py050n","py090n","py100n","py110n","py120n","py130n","py140n"),
#'   conditional = c("netIncomeCat", "pl030"), replaceEmpty = "sequential", seed=1 )
#'
#' class(simPopObj)
#' }
simComponents <- function(simPopObj, total="netIncome",
  components = c("py010n", "py050n", "py090n", "py100n", "py110n", "py120n", "py130n", "py140n"),
  conditional=c(getCatName(total), "pl030"), replaceEmpty=c("sequential", "min"), seed) {

  ##### initializations
  if ( !missing(seed) ) {
    set.seed(seed)
  }

  dataP <- simPopObj@pop@data
  dataS <- simPopObj@sample@data
  w <- simPopObj@pop@strata

  if ( length(total) != 1 ) {
    stop("currently only one continuous variable can be split at a time\n")
  }
  if ( !isTRUE(length(components) >= 2) ) {
    stop("at least two components are required for this procedure!\n")
  }
  varNames <- c(w=w, total=total, components, conditional)
  replaceEmpty <- match.arg(replaceEmpty)
  N <- nrow(dataP)

  # check data
  if ( all(varNames %in% colnames(dataS)) ) {
    dataS <- dataS[, varNames, with=F]
  } else {
    stop("undefined variables in the sample data!\n")
  }
  if ( !all(c(total, conditional) %in% colnames(dataP)) ) {
    stop("undefined variables in the population data!\n")
  }

  # exclude observations
  include <- function(x) !(is.na(x) | x == 0)
  dataS <- dataS[include(dataS[[total]]),]
  indPop <- (1:N)[include(dataP[[total]])]
  dataPop <- dataP[indPop,]

  # data frame dataFrac contains the fractions of the components
  #dataFrac <- prop.table(as.matrix(dataS[, components, with=FALSE]), 1)
  dataFrac <- dataS[,.SD, .SDcols=c(components, total)]
  dataFrac[[total]] <- abs(dataFrac[[total]])
  dataFrac <- dataFrac[,lapply(.SD, function(x) { x/dataFrac[[total]]}), .SDcols=components]
  dataFrac <- as.matrix(dataFrac)
  # matrix simFrac stores the simulated fractions
  simFrac <- matrix(ifelse(is.na(dataP[[total]]), NA, 0), nrow=N, ncol=length(components))
  colnames(simFrac) <- components

  if ( length(conditional) > 0 ) {
    ## class tables
    tabS <- table(dataS[, conditional, with=FALSE])
    tabP <- table(dataPop[, conditional, with=FALSE])
    if ( !identical(dimnames(tabS), dimnames(tabP)) ) {
      stop("incompatible factor levels in conditioning variables!\n")
    }

    ## replace empty combinations in sample by nearest neighbours
    indS <- 1:length(tabS)
    empty <- which(tabS == 0)  # empty cells
    if ( length(empty) ) {
      donors <- indS[-empty]  # donors
      indTabS <- expand.grid(lapply(dim(tabS), function(k) 1:k)) # indices
      indEmpty <- indTabS[empty, , drop=FALSE]  # indices of empty cells
      indDonors <- indTabS[donors, , drop=FALSE]  # indices of donors

      # reorder donors such that last variable varies fastest
      ord <- do.call("order", indDonors[, names(indDonors), drop=FALSE])
      indDonors <- indDonors[ord, , drop=FALSE]
      donors <- donors[ord]

      # compute replacement cells
      if ( replaceEmpty == "sequential" && length(conditional) == 1 ) {
        replaceEmpty <- "min"
        warning("sequential search of replacement cells not applicable for only one conditioning variable!\n")
      }
      fun <- switch(replaceEmpty, sequential=seqMinDist, min=minDist)
      replace <- apply(indEmpty, 1, fun, indDonors, donors)

      # replace with donor cell
      indS[empty] <- replace
    }

    # skip empty combinations in population
    indP <- which(tabP > 0)
    indS <- indS[indP]

    # split the indices of population data by conditioning variables
    indSplitP <- split(indPop, dataPop[, conditional, with=FALSE])
    # split the indices of sample data by conditioning variables
    indSplitS <- split(1:nrow(dataS), dataS[, conditional, with=FALSE])
    # split weights by conditioning variables
    weights <- split(dataS[[w]], dataS[, conditional, with=FALSE])

    # sample indices
    indices <- lapply(1:length(indP), function(i) {
      indicesS <- indSplitS[[indS[i]]]
      if ( length(indicesS) > 1 ) {
        sample(indicesS, size=length(indSplitP[[indP[i]]]), replace=TRUE, prob=weights[[indS[i]]])
      } else {
        rep.int(indicesS, length(indSplitP[[indP[i]]]))
      }
    })

    # replicate the fractions of the components
    simFrac[unlist(indSplitP),] <- dataFrac[unlist(indices),]
  } else {
    if ( nrow(dataS) > 1 ) {
      indices <- spSample(length(indPop), dataS[[w]])
    } else {
      indices <- rep.int(1, length(indPop))
    }
    simFrac[indPop,] <- dataFrac[indices,]
  }

  out <- abs(dataP[[total]]) * simFrac
  for ( i in 1:ncol(out) ) {
    dataP[,colnames(out)[i]] <- out[,i]
  }
  simPopObj@pop@data <- dataP
  invisible(simPopObj)
}
