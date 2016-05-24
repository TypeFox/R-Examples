selectHMDdata <-
function(country=NULL, 
                          data=c("Population",
                            "Deaths",
                            "Exposures",
                            "Rates"), 
                          sex=c("Females",
                            "Males",
                            "Total"), 
                          ages=NULL, 
                          years=NULL){
  
  ## Input: 
  ## country = Character string for the country name
  ## data = Character string showing type of
  ##        demographic series
  ## sex = Character string showing sex of
  ##       demographic series
  ## ages = Vector of ages to extract from data
  ## years = Vector of years to extract from data
  
  ## Output: 
  ## result = a matrix with a subset of country,
  ##          type of data, sex, ages and years
  ##          and country-data-sex attributes 
  ##          and selected ages and years
  ##          as row/colnames
  HMDdata <- NULL
  data(HMDdata, envir = environment())
  ## about country:
  ## if the user doesn't provide the country name
  ## stop the function and gives an error message
  if(is.null(country)) stop("Select a country")
  ## position of the choosen country
  pos.country <- charmatch(tolower(country),
                           tolower(names(HMDdata)))
  ## if the user unclear country-name
  ## (e.g. po for Poland or Portugal)
  if(pos.country==0) stop("Provide distinguishable country name")
  ## if the user provides wrong country name
  if(is.na(pos.country)) stop("Country not found")
  ## name of the country itself
  country.hat <- names(HMDdata)[pos.country] 
  
  ## about data:
  ## if the user doesn't provide the kind of data
  ## 1 - give a warning and automatically select the first
  ##     (i.e. Population)
  if(missing(data)) warning("Population data automatically selected")
  ## choosing the type of data
  data1 <- match.arg(data)
  ## position of the choosen type of data
  pos.data <- switch(data1,Population=5,
                     Deaths=6,Exposures=7,Rates=8)
  ## kind of data
  data.hat <- c("Population", "Deaths",
                "Exposures", "Rates")[pos.data-4]
    
  ## about sex:
  ## if the user doesn't provide the sex
  ## 1 - give a warning and automatically select the first
  ##     (i.e. Females)
  if(missing(sex)) warning("Female data automatically selected")
  ## choosing the sex
  sex1 <- match.arg(sex)
  ## position of the choosen sex
  pos.sex <- switch(sex1,Females=1,Males=2,Total=3)
  ## chosen sex
  sex.hat <- c("Females", "Males", "Total")[pos.sex]
  
  ## selected data with all ages and years
  result0 <- HMDdata[[pos.country]][[pos.data]][[pos.sex]]
  
  ## checking for ages and years out of the available ranges
  ## all possible ages
  all.ages <- as.numeric(rownames(result0))
  ## all available years
  all.years <- as.numeric(colnames(result0))
  ## if provided ages are out of the range
  age.check <- any(is.na(match(ages, all.ages)))
  ## if provided years are out of the range
  year.check <- any(is.na(match(years, all.years)))
  ## printing errors
  stop.ages <- "ages out of the available ranges (0-110)"
  stop.years <- paste("years out of the available ranges (",
                      all.years[1], "-",
                      all.years[length(all.years)], ")",
                      sep="")
  ## possible errors and massociate error messages
  if(age.check & year.check) stop(paste(stop.ages, stop.years, sep="\n       "))
  if(age.check & !year.check) stop(stop.ages)
  if(!age.check & year.check) stop(stop.years)
    
  ## about ages:
  ## if the user doesn't provide ages:
  ## 1 - give a warning
  if(missing(ages)) warning("All ages automatically selected (0-110)")
  ## 2 - take all the ages
  if(missing(ages)){ages1 <- all.ages}else{ages1 <- ages}
  ## which ages have been selected
  whi.ages <- which(!is.na(match(rownames(result0), ages1)))
  
  ## about years:    
  ## if the user doesn't provide years
  ## 1 - give a warning
  if(missing(years)) warning(paste("All available years automatically selected (", 
                                   all.years[1],  "-",
               all.years[length(all.years)], ")", sep=""))
  ## 2 - take all the available years
  if(missing(years)){
    years1 <- all.years}
  else{
    years1 <- years
  }
  ## which years have been selected
  whi.years <- which(!is.na(match(colnames(result0),
                                  years1)))
  
  ############## final matrix ...
  result <- result0[whi.ages, whi.years]
  ## ... and its attributes
  attr(result, "country-data-sex") <- paste(country.hat,
                                            data.hat,
                                            sex.hat,
                                            sep="-")
  class(result) <- "HMDdata"
  result
}
