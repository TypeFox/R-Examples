##' Construct year to year change
##'
##' A function for constructing year to year change
##'
##' @param origVar The variable in which the year to year change is to
##' be calculated
##' @param country The column representing the index of country.
##' @param year The column represing the index of year.
##' @param data The data frame containing the data
##' @param newVarName The name assigned to the new variable, if
##' missing then .CH will be appended.
##' @param n The period for the change rate to be calculated.
##' @return A data frame containing the computed year to year change
##' rate.
##' @export

chConstruct = function(data, origVar, country = "FAOST_CODE", year = "Year",
    newVarName = NA, n = 1){
  tmp = arrange(subset(data, select = c(country, year, origVar)),
                 get(country), get(year))
  unqCountry = unique(tmp[, country])
  chVar = double()
  for(i in 1:length(unqCountry)){
    tmpgr = chgr(as.numeric(unlist(subset(tmp, select = origVar,
        subset = get(country) == unqCountry[i]))), n)
    chVar = c(chVar, tmpgr)
  }
  ch.df = data.frame(tmp[c(country, year)], chVar)
  colnames(ch.df)[3] <- ifelse(is.na(newVarName),
                               paste(origVar, ".CH", n, sep = ""),
                               newVarName)
  ch.df
}
