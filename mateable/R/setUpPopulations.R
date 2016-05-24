##' simulateScene generates a matingScene object -- a simulated population
##' in a standard format with individuals randomly assigned a mating schedule,
##' a location, and S-alleles
##'
##' @title Simulate a Mating Scene
##' @param size integer number of plants
##' @param meanSD date mean start date
##' @param sdSD date standard deviation of start date
##' @param skSD skew of the start date of the population
##' @param meanDur numeric duration in days
##' @param sdDur standard deviation of duration in days
##' @param xRange range of spatial extent of individuals along x-axis
##' @param yRange range of spatial extent of individuals along y-axis
##' @param distro unimplemented
##' @param sAlleles integer count of S-Alleles that could be in the population
##'
##' @return matingScene data frame -- see \code{\link{makeScene}}
##' @seealso \code{\link{makeScene}}
##' @author Stuart Wagenius
##' @examples
##' simulateScene()
##' \dontrun{simulateScene(NULL)}
simulateScene <- function(size = 30, meanSD = "2012-07-12", sdSD = 6, meanDur = 11,
                          sdDur = 3, skSD = 0 ,xRange = c(0, 100), yRange = c(0, 100),
                          distro = "unif", sAlleles = 10) {
  md <- as.integer(as.Date(meanSD, "%Y-%m-%d"))
  sd <- as.integer(md + round(sn::rsn(n = size, 0, omega = sdSD, alpha = skSD), 0))
  ed <- as.integer(sd + abs(round(rnorm(size, meanDur, sdDur), 0)))

  if (distro != "unif")
    warning("distro must be unif")
  xv <- runif(size, min = xRange[1], max = xRange[2])
  yv <- runif(size, min = yRange[1], max = yRange[2])
  sM <- sample(x = 1:sAlleles, size = size, replace = TRUE)
  if (sAlleles == 2) {
    sP <- 3 - sM
  } else {
    sP <- sapply(sM, FUN = function(x) sample((1:sAlleles)[-x], 1))
  }
  df <- data.frame(id = 1:size, start = sd, end = ed, x = xv,
                   y = yv, s1 = sM, s2 = sP)
  makeScene(df, startCol = "start", endCol = "end", xCol = "x", yCol = "y",
            idCol = "pla", dateFormat = "1970-01-01")
}

##' Turns a data frame with information about temporal, spatial, or
##' genetic mating data into a matingScene object using a standard format.
##'
##' @title Create a matingScene object from a data frame
##' @param df a data frame containing information about a mating scene,
##' namely coordinate of individuals in space, time, and mating type.
##' @param multiYear logical indicating whether or not to split the result into
##' a list by year
##' @param startCol character name of column with start dates
##' @param endCol character name of column with end dates
##' @param xCol character name of column with x or E coordinates
##' @param yCol character name of column with y or N coordinates
##' @param s1Col character name of one column with S-allele
##' @param s2Col character name of another column with S-alleles
##' @param idCol character name for column with unique identifier
##' @param dateFormat character indicating either (1) the format of the start and end
##' date columns if those columns are characters or (2) the origin for the start
##' and end date columns if those columns are numeric. It is used in as.Date
##'
##' @return a matingScene object, either a single dataframe in standard format
##' or a list of dataframes. Attributes of the matingScene object indicate the type of
##' information in the data frame, including the original column names,
##' and the origin of the date columns. If multiYear = TRUE,
##' the return value will be a list of matingScene data frames where each
##' element in the list represents one year. See details for more information
##' on attributes and how to work with multi-year data.
##' @details The input dataframe can contain information about locations of
##' individuals in 1, 2, or 3 dimensions of a mating scenes.
##' The function currently allows two spatial coordinates. The user specifies
##' the names of the columns and they will be saved xCol and yCol in the
##' matingScene object. MatingScene objects currently save temporal
##' coordinates for individuals as start and end date of mating activity
##' within a year. Mating type coordinates are saved as mating type alleles.
##' Columns are named id, start, end, x, y, s1, and s2 for
##' idCol, startCol, endCol, xCol, yCol, s1Col, and s2Col respectively.
##' The attributes "t", "s", and "mt" will be set to TRUE if the data frame
##' has temporal, spatial, or mating type data, respectively and
##' will be FALSE otherwise. The attribute originalNames contains all the
##' names of the columns in the original data frame.\cr
##' The start and end columns will be changed to integers relative to the start
##' day of the population. So the first day of the first individual to become
##' receptive will be 1 and so on. The attribute origin contains the
##' origin that can be used when converting the columns start and end
##' from integers to dates.\cr
##' If no temporal data are available except the year in which it was
##' collected and df is a multi-year data set, put the collection year into the
##' column labelled as startCol and set dateFormat = "%Y" and that will split
##' the data appropriately.
##' @author Danny Hanson
makeScene <- function (df, multiYear = FALSE, startCol = "start", endCol = "end", xCol = "x",
                       yCol = "y", s1Col = "s1", s2Col = "s2", idCol = "id",
                       dateFormat = "%Y-%m-%d") {
  if (multiYear) {
    if (dateFormat == "%Y") {
      dates <- as.Date(as.character(df[, startCol]), dateFormat)
    } else {
      dates <- as.Date(df[, startCol], dateFormat)
    }
    df$year <- as.numeric(format(dates, "%Y"))
    years <- levels(as.factor(df$year))
    newScene <- list()
    for (i in 1:length(years)) {
      newScene[[as.character(years[i])]] <-
        makeScene(df[df$year %in% years[i],], F, startCol, endCol, xCol, yCol,
                  s1Col, s2Col, idCol, dateFormat)
    }
  } else {
    newScene <- data.frame(id = character(nrow(df)))

    if (idCol %in% names(df)) {
      newScene$id <- df[, idCol]
    } else {
      newScene$id <- 1:nrow(df)
    }

    attr(newScene, "t") <- FALSE
    attr(newScene, "s") <- FALSE
    attr(newScene, "mt") <- FALSE
    attr(newScene, "originalNames") <- names(df)

    if (all(c(startCol, endCol) %in% names(df))) {
      attr(newScene, "t") <- TRUE
      newScene$start <- as.integer(as.Date(df[, startCol], dateFormat))
      firstDay <- min(newScene$start)
      newScene$start <- newScene$start - firstDay + 1
      newScene$end <- as.integer(as.Date(df[, endCol], dateFormat)) - firstDay + 1
      newScene$duration <- newScene$end - newScene$start + 1
      origin <- as.Date(firstDay-1, "1970-01-01")

      attr(newScene, "origin") <- origin
    }

    if (all(c(xCol, yCol) %in% names(df))) {
      attr(newScene, "s") <- TRUE
      newScene$x <- df[, xCol]
      newScene$y <- df[, yCol]
    }

    if (all(c(s1Col, s2Col) %in% names(df))) {
      attr(newScene, "mt") <- TRUE
      newScene$s1 <- as.factor(df[, s1Col])
      newScene$s2 <- as.factor(df[, s2Col])
    }

    # not going to add this for now because it's unlikely we'll make our
    # own generics or use oop
    # class(newScene) <- "matingScene"
  }
  newScene
}
