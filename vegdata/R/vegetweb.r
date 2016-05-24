vw.survey <- function(searchstring, ...) {
#  require(httr)
  surveys <- fromJSON('http://botanik3.botanik.uni-greifswald.de/floradb-rs/service/v1/surveys?portalId=3')$survey
  if(is.numeric(searchstring)) {
    df <- data.frame(Projekt_ID = surveys[grep(searchstring, surveys$id),'id'], Projekttitel = surveys[grep(searchstring, surveys$id),'title'], Kustode= paste(surveys[grepl(searchstring, surveys$title),'owner']$firstName, surveys[grepl(searchstring, surveys$title),'owner']$lastName))
  } else {
    df <- data.frame(Projekt_ID = surveys[grep(searchstring, surveys$title),'id'], Projekttitel = surveys[grep(searchstring, surveys$title),'title'], Kustoden = paste(surveys[grepl(searchstring, surveys$title),'owner']$firstName, surveys[grepl(searchstring, surveys$title),'owner']$lastName) )
  }
  return(df)
}


vw.site <- function(user, password, survey, basket, ...) {
#  require(httr)
  if(!missing(basket)) stop('Webservices for vegetation plot baskets are not yet (Nov 2015) implemented in vegetweb.de')
  message('This is a test implementation. Passwords will be send through http, i.e. unsecured.')
  r <- GET(paste('http://botanik3.botanik.uni-greifswald.de/floradb-rs/service/v1/snapshots', survey[1], '?occurrenceAttribute=COVERAGE_MEAN', sep='/'),  authenticate(user, password), add_headers("Accept : application/json")) 
  stop_for_status(x=r, task = 'access the dataset.')
  data <- content(r, "parsed", "application/json")
  nbplots <- length(unique(sapply(data$data, '[[', 'sampleUUID')))
  message('Number of plots: ', nbplots)
  if(nbplots > 0) {
  fields <- names(data$header[[1]])
  message('Available information: ', paste(fields, collapse=' '))
  site <- do.call(rbind.data.frame, data$header) 
if(length(survey) > 1)
   for(i in 2:length(survey)) {
     r <- GET(paste('http://botanik3.botanik.uni-greifswald.de/floradb-rs/service/v1/snapshots', survey[1],  '?occurrenceAttribute=COVERAGE_MEAN', sep='/'),  authenticate(user, password), add_headers("Accept : application/json"))
	  stop_for_status(x=r, 'access the dataset.')
    data <- content(r, "parsed", "application/json")
    assign(paste(fields,i, sep='.'), names(data$header[[1]]))
    message('Available information: ', paste(fields, collapse=' '))
    site.tmp <- do.call(rbind.data.frame, data$header) 
  	cols1 <- names(site)
  	cols2 <- names(site.tmp)
   	All <- union(cols1, cols2)
  	miss1 <- setdiff(All, cols1)
  	miss2 <- setdiff(All, cols2)
  	site[, c(as.character(miss1))] <- NA
  	site.tmp[,c(as.character(miss2))] <- NA
  	site <- rbind(site, site.tmp)
  }

  ### Survey Area
  if(!'SURF_AREA' %in% names(site)) site$SURF_AREA <- NA
  n <- sum(site$SURF_AREA == 0 | is.na(site$SURF_AREA))
  if(n>0) message(paste(n, ' releves without survey area information.'))
  return(site)
}}


vw.veg <- function(user, password, survey, basket, taxeval = TRUE, ...) {
#  require(httr)
  message('This is a test implementation. Passwords will be send through http, i.e. unsecured.')
  if(is.character(survey)) {
    surveyid <- vw.survey(survey)
    if(nrow(surveyid) > 1) {stop("More than one survey found, please restrict.")}
      surveyid <- surveyid$id
  } else surveyid <- survey
  if(!missing(basket)) stop('Webservices for vegetation plot baskets are not yet (Nov 2015) implemented in vegetweb.de')
  r <- GET(paste('http://botanik3.botanik.uni-greifswald.de/floradb-rs/service/v1/snapshots', surveyid, '?occurrenceAttribute=COVERAGE_MEAN', sep='/'),  authenticate(user, password), add_headers("Accept : application/json"))
  stop_for_status(x=r, 'access the dataset.')
  data <- content(r, "parsed", "application/json")
  nbplots <- length(unique(sapply(data$data, '[[', 'sampleUUID')))
  message('Number of plots: ', nbplots)
  obs <- data.frame(RELEVE_NR = sapply(data$data, '[[', 'sampleUUID'), TaxonUsageID = sapply(data$data, '[[', 'germanSLNo'), COVER_PERC = sapply(data$data, '[[', 'value'), LAYER = sapply(data$data, '[[', 'layer'))
  if(taxeval)  
    obs <- taxval(obs, refl = 'GermanSL 1.3', check.critical = FALSE, ...)
  lc = c("layer"); values = "COVER_PERC"; dec=1
  collab <- as.vector(obs$TaxonUsageID)
  rowlab <- as.vector(obs$RELEVE_NR)
  cat('combining occurrences using type', toupper(lc), 'and creating vegetation matrix ... \n')
  layerfun <- function(x) round((1 - prod(1 - x/100)) * 100, dec)
  results <- tapply(as.numeric(obs[, values]), list(rowlab, collab), layerfun)
  results[is.na(results)] <- 0
  veg <- as.data.frame(results)
  class(veg) <- c("veg", "data.frame")
  attr(veg, 'taxreflist') <- 'GermanSL 1.3'
  return(veg)
}

