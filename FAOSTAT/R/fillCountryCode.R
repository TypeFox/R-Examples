##' A function to get country code when not available in data.
##'
##' This function can be useful when a dataset provided does not have
##' a country code available.
##'
##' @param country The column name of the data which contains
##' the country name
##' @param data The data frame to be matched
##' @param outCode The output country code system, defaulted to FAO standard.
##' @export

fillCountryCode = function(country, data, outCode = "FAOST_CODE"){
  unqCountry = unique(data[, country])
  n = length(unqCountry)
  countryCODE = rep(NA, n)
  for(i in 1:n){
    ind = which(as.matrix(FAOcountryProfile[,
      c("OFFICIAL_FAO_NAME", "SHORT_NAME", "FAO_TABLE_NAME",
        "UNOFFICIAL1_NAME", "UNOFFICIAL2_NAME", "UNOFFICIAL3_NAME")]) ==
      unqCountry[i], arr.ind = TRUE)
    which.row = ind[, 1]
    if(length(unique(which.row)) == 1)
      countryCODE[i] = FAOcountryProfile[unique(which.row), outCode]
  }
  if(anyDuplicated(na.omit(countryCODE)))
    warning(paste0("Duplicated ", outCode, " matched, double check the data"))
  if(any(is.na(countryCODE)))
    warning(paste0("Certain ", outCode , " were not matched."))
  warning("Please check the correct China has been specified.")
  def = data.frame(unqCountry, countryCODE)
  colnames(def) = c(country, outCode)
  merge(x = data, y = def, by = country, all.x = TRUE)
}

utils::globalVariables(names = c("FAOcountryProfile"))