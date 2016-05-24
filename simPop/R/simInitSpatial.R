#' Generation of smaller regions given an existing spatial variable and a
#' table.
#' 
#' This function allows to manipulate an object of class
#' \code{\linkS4class{simPopObj}} in a way that a new variable containing
#' smaller regions within an already existing broader region is generated. The
#' distribution of the smaller region within the broader region is respected.
#' 
#' The distributional information must be contained in an input table that
#' holds combinations of characteristics of the broader region and the smaller
#' regions as well as population counts (which may be available from a census).
#' 
#' @name simInitSpatial
#' @param simPopObj an object of class \code{\linkS4class{simPopObj}}.
#' @param additional a character vector of length one holding the variable name
#' of the variable containing smaller geographical units. This variable name
#' must be available as a column in input argument \code{tspatial}.
#' @param region a character vector of length one holding the variable name of
#' the broader region. This variable must be available in the input
#' \code{tspatial} as well as in the sample and population slots of input
#' \code{simPopObj}.
#' @param tspatial a data.frame containing three columns. The broader region
#' (with the variable name being the same as in input \code{region}, the
#' smaller geographical units (with the variable name being the same as in
#' input \code{additional} and a third column containing a numeric vector
#' holding counts.))
#' @return An object of class \code{\linkS4class{simPopObj}} with an additional
#' variable in the synthetic population slot.
#' @author Bernhard Meindl
#' @keywords manip
#' @export
#' @examples
#' data(eusilcS)
#' data(eusilcP)
#' 
#' # no districts are available in the population, so we have to generate those
#' # we randomly assign districts within "region" in the eusilc population data
#' # each hh has the same district
#' simulate_districts <- function(inp) {
#'   hhid <- "hid"
#'   region <- "region"
#' 
#'   a <- inp[!duplicated(inp[,hhid]),c(hhid, region)]
#'   spl <- split(a, a[,region])
#'   regions <- unique(inp[,region])
#' 
#'   tmpres <- lapply(1:length(spl), function(x) {
#'     codes <- paste(x, 1:sample(3:9,1), sep="")
#'     spl[[x]]$district <- sample(codes, nrow(spl[[x]]), replace=TRUE)
#'     spl[[x]]
#'   })
#'   tmpres <- do.call("rbind", tmpres)
#'   tmpres <- tmpres[,-c(2)]
#'   out <- merge(inp, tmpres, by.x=c(hhid), by.y=hhid, all.x=TRUE)
#'   invisible(out)
#' }
#' 
#' eusilcP <- simulate_districts(eusilcP)
#' table(eusilcP$district)
#' 
#' # we generate a synthetic population
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' simPopObj <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#' 
#' # we generate the input table using the broad region (variable 'region')
#' # and the districts, we have generated before.
#' # we
#' tab <- as.data.frame(xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$region + eusilcP$district))
#' colnames(tab) <- c("db040", "district", "Freq")
#' 
#' simPopObj <- simInitSpatial(simPopObj, additional="district", region="db040", tspatial=tab)
#' 
simInitSpatial <- function(simPopObj, additional, region, tspatial) {
  # simplified version of generateValues_distribution
  generateValues_spatial <- function(dataTable, dataPop, params) {
    if( !nrow(dataTable) ) {
      return(character())
    }
    grid <- expand.grid(dataTable[,params$additional])
    
    # percentages
    perc <- dataTable$freq / sum(dataTable$freq)
    
    # draw and exapand
    sizes <- dataPop[,.N, key(dataPop)]
    sim <- rep(sample(grid[,1], nrow(sizes), prob=perc, replace=TRUE), sizes$N)
    sim
  }
  
  x <- NULL
  dataP <- simPopObj@pop
  dataS <- simPopObj@sample
  data_pop <- dataP@data
  data_sample <- dataS@data
  basic <- simPopObj@basicHHvars

  if ( length(additional) != 1 ) {
    stop("currently exactly one additional spatial variable can be generated!\n")
  }
  if ( length(region) != 1 ) {
    stop("argument 'region' must be a character vector of length 1!\n")
  }  
  if ( additional %in% colnames(data_pop) ) {
    stop("variable specified in argument 'additional' already exists in the population!\n")
  }
  if ( ncol(tspatial) != 3 ) {
    stop("please check input 'tspatial'! It must have exactly 3 columns!\n")
  }
  
  freqs <- tspatial[,ncol(tspatial)]
  if ( !is.numeric(freqs) ) {
    stop("last column of input table must contain numeric values!\n")
  }
  tspatial <- tspatial[,-ncol(tspatial), drop=F]
  
  m <- match(additional, colnames(tspatial))
  if ( is.na(m) ) {
    stop("variable specified in argument 'additional' (",additional,") is not available in input table 'tspatial'!\n")
  }
  add <- tspatial[,m]
  
  m <- match(region, colnames(tspatial))
  if ( is.na(m) ) {
    stop("variable specified in argument 'additional' (",region,") is not available in input table 'tspatial'!\n")
  }
  reg <- tspatial[,m]
  
  # check other variables levels
  m <- match(region, colnames(data_pop))
  if ( is.na(m) ) {
    stop("variable listed in argument 'region' is not available in the synthetic population data of of input 'simPopObj'!\n")
  }
  m <- match(region, colnames(data_sample))
  if ( is.na(m) ) {
    stop("variable listed in argument 'region' is not available in the sample dataset of input 'simPopObj'!\n")
  }  
  
  # generation of our table
  tab <- data.frame(reg, add, freqs)
  colnames(tab) <- c(region, additional, "freq")  
  
  for ( i in 1:2 ) {
    a <- sort(unique(as.character(tab[,i])))
    m <- match(colnames(tab)[i], colnames(data_pop))
    b <- sort(unique(as.character(data_pop[[m]])))    
    if ( any(a!=b) ) {
      stop("We fould a problem in variable ", colnames(tspatial)[i],". Values in input table and synthetic population do not match!\n")
    }
  }

  # list indStrata contains the indices of dataP split by region
  N <- nrow(data_pop)
  indStrata <- split(1:N, data_pop[[region]])
  
  # predictor variables
  predNames <- dataP@hhid  # in spatial case, it can only be the hhid
  
  params <- list()
  params$additional <- additional
  
  values <- lapply(levels(data_sample[[region]]), function(x) {
    generateValues_spatial(
      dataTable=subset(tab, tab[,region]==x),
      dataPop=data_pop[indStrata[[x]], predNames, with=F], params)
  })
  
  ## add new categorical variables to data set and return
  data_pop[[additional]] <- unlist(values)
  
  # check
  simPopObj@pop@data <- data_pop
  invisible(simPopObj)
}
