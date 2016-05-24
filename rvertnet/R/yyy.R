# Wrapper for vertsummary function
vertsumwrapper <- function(input = NULL, verbose = TRUE){
  
  if (!class(input) %in% c("list", "data.frame")) {
    stop("Input must be of class list or data.frame", call. = FALSE)
  }
  if (is(input, "list"))  input <- input$data
  
  # recs <- number of records in the data frame
  recs <- nrow(input)
  
  # coords <- number of records with viable lat and long data
  # errest <- number of "coords" records with viable coordinate uncertainty estimate
  if(is.null(input$decimallatitude) & is.null(input$decimallongitude)){
    coords <- 0
  } else{ 
    coords <- NULL
  }
  if(inherits(input$coordinateuncertaintyinmeters, "NULL")){
    errest <- 0
  } else {
    errest <- NULL
  }
  if (is.null(coords)) {
    coords <- sum(complete.cases(input[,c('decimallatitude','decimallongitude')]))
    # checking for good lat/long data (if not, use only the above line)
    input$decimallatitude <- as.numeric(as.character(input$decimallatitude))
    input$decimallongitude <- as.numeric(as.character(input$decimallongitude))
    if(is.null(errest)){
      input$coordinateuncertaintyinmeters <- as.numeric(as.character(input$coordinateuncertaintyinmeters))
    }
    mappable <- input[complete.cases(input[,c('decimallatitude','decimallongitude')]),]
    mappable <- subset(mappable, input$decimallatitude < 90 & input$decimallatitude > -90)
    mappable <- subset(mappable, input$decimallongitude < 180 & input$decimallongitude > -180)
    if(nrow(mappable) < coords){
      bad <- coords - nrow(mappable)
      mssg(verbose, paste(bad, " record(s) with bad coordinates"))
      coords <- coords - bad
    }
    if(is.null(errest)){
      mappable <- subset(mappable, input$coordinateuncertaintyinmeters > 0 &
                           input$coordinateuncertaintyinmeters < 20020000)
      if((errest <- nrow(mappable)) < coords){
        bad <- coords - errest
      }
    }
  }
  
  # instcoll <- number of records from each institution+collection
  removeDups <- function(x) {
    paste(unique(unlist(strsplit(x, split=" "))), collapse = " ")
  }
  if(inherits(input$institutioncode, "NULL") & inherits(input$collectioncode, "NULL")){
    instcoll <- NA
  } else {
    instcoll <- as.matrix(paste(input[,'institutioncode'], 
                                input[,'collectioncode'], sep = " "))
    instcoll <- table(apply(instcoll, 1, removeDups))
  }
  
  # country <- number of records from each country
  if(inherits(input$country, "NULL")){country <- NA} else{country <- table(input$country)}
  
  # year <- number of records by year
  if(inherits(input$year, "NULL")){year <- NA} else{year <- table(input$year)}
  
  # taxon <- number of records by taxonomic name
  taxon <- as.matrix(paste(input[,'genus'], input[,'specificepithet'], sep = " "))
  if(!inherits(input$infraspecificepithet, "NULL")){
    taxon <- as.matrix(paste(taxon, input[,'infraspecificepithet'], sep = " "))
  }
  taxon <- gsub(" NA", "", taxon) # remove unknowns - usually infrasp.ep
  taxon <- table(apply(taxon, 1, removeDups))
  
  # return summary
  structure(list("recs" = recs, "coords" = coords, "errest" = errest, "instcoll" = instcoll,
                 "country" = country, "year" = year, "taxon" = taxon), class="vertsummary")
}

#' @export
print.vertsummary <- function(x, ...){
  cat(paste0("Number of records ($recs): ", x$recs), sep = "\n") 
  cat(paste("Records with decimal lat/long (-90<lat<90, -180<long<180) ($coords): ", x$coords, sep = ""), sep = "\n")
  cat(paste("Records with lat/long and coordinate uncertainty estimate (0<errest<20020000) ($errest): ", x$errest), sep = "\n")
  cat("Record count by institution/collection ($instcoll): ", sep = "\n")
  print(x$instcoll)
  cat("\nRecord count by country ($country): ", sep = "\n")
  print(x$country)
  cat("\nRecord count by year ($year): ", sep = "\n")
  print(x$year)
  cat("\nRecord count by taxon ($taxon): ", sep = "\n")
  print(x$taxon) 
}
