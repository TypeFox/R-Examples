select_cont <- function(dat, nCont, type, areaVar, fixed) {
  datList <- makeDataList(dat, areaVar)
  obs <- if (type == "unit") {
    makeObsUnit(nCont, sapply(datList, nrow))
  } else {
    makeObsArea(nCont, sapply(datList, nrow), fixed)
  }
  mapply(selectObsFromData, datList, obs, fixed, SIMPLIFY = FALSE) %>% rbind_all
}

makeDataList <- function(dat, areaVar) {
  if (is.null(areaVar)) {
    list(dat)
  } else {
    split(dat, dat[areaVar])
  }
}

makeObsUnit <- function(nCont, nrowsOfDatList) {
  if (any(nCont %% 1 != 0)) { # nCont is a vector of probs
    mapply(rbinom, nrowsOfDatList, 1, nCont, SIMPLIFY = FALSE) %>% sapply(sum)
  } else {
    rep(nCont, length.out = length(nrowsOfDatList))
  }
}

makeObsArea <- function(nCont, nrowsOfDatList, fixed) {
  
  # if nCont is not a scalar, then expect a vector with positions
  isScalar <- length(nCont) == 1 
  
  numberOfAreas <- if (isScalar && (nCont %% 1) != 0) { # nCont is a prob
    rbinom(length(nrowsOfDatList), 1, nCont) %>% sum
  } else if (isScalar) {
    nCont
  } else {
    # Then numberOfAreas is never used.
    NULL
  }
  
  # Selection of non contaminated obs:
  selection <- if (fixed && isScalar) {
    1:max(length(nrowsOfDatList) - numberOfAreas, 1)
  } else if (isScalar) {
    sample.int(length(nrowsOfDatList), max(0, length(nrowsOfDatList) - numberOfAreas))
  } else {
    seq_along(nrowsOfDatList)[-nCont]
  }
  
  nrowsOfDatList[selection] <- 0
  nrowsOfDatList
}

selectObsFromData <- function(dat, obs, fixed) {
  selection <- if (fixed) {
    if (obs == 0) numeric(0) else (max(nrow(dat) - obs + 1, 1)):nrow(dat)
  } else {
    sample.int(nrow(dat), min(obs, nrow(dat)))
  }
  dat[!((1:nrow(dat)) %in% selection), ] <- 0
  dat$idC <- (1:nrow(dat)) %in% selection
  dat
}
