#' Gapminder data source.
#'
#' Gapminder describes itself as a "fact tank" that promotes a fact based world view. 
#' On their website they provide a service that allows to create animated charts for various indicators, differentiated by country.
#' They also provide the underlying datasets for download. This S3 class serves as a wrapper for easy access to a subset of these data.
#'
#' Please note that neither Gapminder nor the package developer/maintainer are the data provider, except for a few cases.
#' Therefore you will have to go to the source to find out the terms of use for the specific indicator.
#'
#' This class defines some resources of the Gapminder Project.
#' See \code{queries(gapminder())} for a list of resources.
#'
#' This is a proof of concept for the \code{urldata} function.
#'
#' @seealso \code{\link{urldata}}
#'
#' @examples
#' \dontrun{
#'   gm <- gapminder()
#'   queries(gm)
#'   query(gm, "ReligionAndBabies")
#' }
#' @return Mashup
#' @references 
#' \url{http://www.gapminder.org}
#' @export
gapminder <- function() datamart(
  res_gapminder("Population", "phAwcNAVuyj0XOoBL_n5tAQ&gid=0"),
  res_gapminder("MainReligion", "0ArtujvvFrPjVdHUzTGVicFJZQ1NjaFhqelV5SDNxMVE&gid=1"),
  res_gapminder("TotalFertilityRate", "phAwcNAVuyj0TAlJeCEzcGQ&gid=0"),
  res_gapminder("PerCapitaCO2Emissions", "phAwcNAVuyj1gkNuUEXOGag&gid=0"),
  res_gapminder("IncomePerCapita", "phAwcNAVuyj1jiMAkmq1iMg&gid=0"),
  res_gapminder("InfantMortalityRate", "phAwcNAVuyj0NpF2PTov2Cw&gid=0"),
  res_gapminder("LifeExpectancyAtBirth", "phAwcNAVuyj2tPLxKvvnNPA&gid=0"),
  res_gapminder("AdolescentFertilityRate", "pyj6tScZqmEdIphYUHxcdLg&gid=0"),
  res_gapminder("BirthsAttendedBySkilledHealthStaff", "pyj6tScZqmEfKY9bk02DBYA&gid=0"),
  res_gapminder("ContraceptiveUse", "pyj6tScZqmEewsQOoKrtYJQ&gid=0"),
  res_gapminder("CrudeBirthRate", "tUSeGJOQhafugwUvHvY-wLA&gid=0"),
  res_gapminder("MaternalMortalityRate", "pyj6tScZqmEcVezxiMlWaRw&gid=0"),
  res_gapminder("Under5MortalityRate", "phAwcNAVuyj05ZR69usyQIg&gid=0"),
  res_gapminder("CrudeDeathRate", "tHyj-2jRvK3CCNJOc5Vm-HQ&gid=0"),
  res_gapminder("PopulationGrowth", "pyj6tScZqmEcl2xDWbuJ8fg&gid=0"),
  res_gapminder("SugarConsumption", "phAwcNAVuyj2sdmdhX9zuKg&gid=0"),
  res_gapminder("GDP", "pyj6tScZqmEfI4sLVvEQtHw&gid=0"),
  res_gapminder("ConsumerPricesIndex", "pyj6tScZqmEc3xNIyXiZ6EA&gid=0"),
  res_gapminder("GDPImplicitDeflator", "pyj6tScZqmEcaHt8Y6cxXQg&gid=0"),
  res_gapminder("CoalConsumption", "pyj6tScZqmEc1TmMiFdmOVg&gid=0"),
  res_gapminder("HydroelectricityConsumption", "pyj6tScZqmEdNbX4qj9QLTA&gid=0"),
  res_gapminder("NaturalGasConsumption", "pyj6tScZqmEcx9pD804Q0Aw&gid=0"),
  res_gapminder("NuclearConsumption", "pyj6tScZqmEfiy57wnt-tEA&gid=0"),
  res_gapminder("OilConsumption", "pyj6tScZqmEcm0fIa0IVtKw&gid=0"),
  res_gapminder("CoalProduction", "pyj6tScZqmEdDid2ts7KvHg&gid=0"),
  res_gapminder("ElectricityGeneration", "pyj6tScZqmEehRG-9mMHYdg&gid=0"),
  res_gapminder("NaturalGasProduction", "pyj6tScZqmEfv2K6dZmskWg&gid=0"),
  res_gapminder("OilProduction", "pyj6tScZqmEdNIa3ckVXaCQ&gid=0"),
  res_gapminder("PrimaryEnergyConsumption", "pyj6tScZqmEeTCOezV8a3HA&gid=0"),
  res_gapminder("CO2Emissions", "phAwcNAVuyj1NHPC9MyZ9SQ&gid=0"),
  res_gapminder("SulfurEmissions", "t9SYWh7siLJDzyZYN1R4HfQ&gid=0"),
  res_gapminder("TotalForestArea", "pp59adS3CHWeB1N1HlpFQVQ&gid=0"),
  res_gapminder("PrimaryForestArea", "pp59adS3CHWeECA6Gf__BNQ&gid=0"),
  res_gapminder("PlantedForestArea", "pp59adS3CHWc4aJd9fV8zZg&gid=0"),
  res_gapminder("WoodRemoval", "pp59adS3CHWe8O-N9RgxzDw&gid=0"),
  res_gapminder("BiomassStockInForest", "pp59adS3CHWcsSl830EklJA&gid=0"),
  res_gapminder("TotalWaterWithdrawal", "rIG3ZWxv381t2bIL2BNaIVw&gid=0"),
  res_gapminder("SurfaceArea", "pyj6tScZqmEeiMy8j86qDTg&gid=0"),
  res_gapminder("BadTeethPerChild", "phAwcNAVuyj3Os9LVO_pRDA&gid=0"),
  res_gapminder("PeopleLivingWithHIV", "pyj6tScZqmEe1GaiYJX2qGA&gid=0"),
  res_gapminder("MalariaReportedCases", "pp59adS3CHWczfPHQMiqxCg&gid=0"),
  res_gapminder("MalariaReportedDeaths", "pp59adS3CHWfZGL9qouvTbQ&gid=0"),
  res_gapminder("WorkingHoursPerWeek", "rIMebcn9Eo2jSIm09HBLihg&gid=0"),
  res_gapminder("UrbanPopulation", "pyj6tScZqmEfH89V6UQhpZA&gid=0"),
  res_gapminder("WomensAgeAtFirstMarriage", "t4eF8H_jq_xyKCUHAX6VT1g&gid=0"),
  res_gapminder("NumberOfBillionaires", "tNWhbu-1UIPPxtmRHtnINOQ&gid=0"),
  res_gapminder("GiniIndex", "pyj6tScZqmEcjeKHnZq6RIg&gid=0"),
  res_gapminder("BroadbandSubscribers", "pyj6tScZqmEcuy6dYkzGhfw&gid=0"),
  res_gapminder("CellPhones", "pyj6tScZqmEcKuNdFCUo6TQ&gid=0"),
  res_gapminder("PersonalComputers", "pyj6tScZqmEfUXdC83YSzfw&gid=0"),
  res_gapminder("PatentApplications", "pyj6tScZqmEd5FA9xlfO9eA&gid=0"),
  res_gapminder("PatentsGranted", "pyj6tScZqmEdMioz5VJKXHw&gid=0"),
  res_gapminder("PatentsInForce", "pyj6tScZqmEe371ZVZl73eA&gid=0"),
  res_gapminder("ArmsExports", "pyj6tScZqmEeTIhjRrVQtQA&gid=0"),
  res_gapminder("ArmsImports", "pyj6tScZqmEfnPl7VRfT9WA&gid=0"),
  res_gapminder("HumanDevelopmentIndex", "tyadrylIpQ1K_iHP407374Q&gid=0"),
  
  # example graphs
  resfunc("ReligionAndBabies", depends=c("TotalFertilityRate", "IncomePerCapita", "MainReligion"), fun=gapm_religion_babies)
)

res_gapminder <- function(resource, code) urldata(
  resource=resource, 
  template=sprintf("https://docs.google.com/spreadsheet/pub?key=%s&output=csv", code),
  extract.fct=function(uri) RCurl::getURL(uri, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")),
  transform.fct=function(x) {
    res <- read.csv(textConnection(x), na.strings=c("..", "-"), stringsAsFactor=FALSE)
    # as date, if possible
    tm <- try(as.Date(paste(substring(tail(colnames(res),-1),2),1,1,sep="-0")), silent=TRUE)
    if(!inherits(tm, "try-error")) { 
      colnames(res)[1] <- "country"
      res <- reshape(res,  idvar="country", varying=list(2:ncol(res)), times=strtail(colnames(res)[2:ncol(res)],-1), v.names=resource, direction="long")
      res <- res[!is.na(res[,resource]),]
    }
    return(res)
  }
)

gapm_religion_babies <- function(TotalFertilityRate, IncomePerCapita, MainReligion, verbose=FALSE, ...) {
  # babies per woman
  idx <- TotalFertilityRate$time=="2008"
  babies <- TotalFertilityRate[idx,"TotalFertilityRate"]
  names(babies) <- TotalFertilityRate[idx,"country"]
  countries <- names(babies)
  
  # income per capita, PPP adjusted
  idx <- IncomePerCapita$time=="2008"
  income <- IncomePerCapita[idx, "IncomePerCapita"]
  names(income) <- IncomePerCapita[idx, "country"]
  countries <- intersect(countries, names(income))
  
  # religion
  religion <- MainReligion[,"Group"]
  names(religion) <- MainReligion[,"Entity"]
  religion[religion==""] <- "unknown"
  colcodes <- c(
    Christian="blue", 
    "Eastern religions"="red", 
    Muslim="green", "unknown"="grey"
  )
  countries <- intersect(countries, names(religion))
  
  # plot
  par(mar=c(4,4,0,0)+0.1)
  plot(
    x=income[countries], 
    y=babies[countries], 
    col=colcodes[religion[countries]], 
    log="x",
    xlab="Income per Person, PPP-adjusted", 
    ylab="Babies per Woman"
  )
  legend(
    "topright", 
    legend=names(colcodes), 
    fill=colcodes, 
    border=colcodes
  )
  return(invisible(NULL))
}

