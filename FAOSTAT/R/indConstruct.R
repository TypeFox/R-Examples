##' Construct indices
##'
##' A function for constructing indices
##'
##' @param data The data frame containing the data
##' @param origVar The variable in which the indices is to be computed
##' @param newVarName The name assigned to the new variable, if missing
##' then .SC/.SH/.GR/.CH/.IND will be appended depending on the type of
##' construction.
##' @param baseYear The year which will serve as the base
##' @return The indice
##' @export
##' @examples
##' test.df = data.frame(FAOST_CODE = rep(1, 100), Year = 1901:2000,
##'                       test = 1:100)
##' indConstruct(test.df, origVar = "test", baseYear = 1950)
##' 
indConstruct = function(data, origVar, newVarName = NA, baseYear = 2000){
  tmp = arrange(subset(data, select = c("FAOST_CODE", "Year", origVar)),
                 FAOST_CODE, Year)
  unqCountry = unique(tmp$FAOST_CODE)
  indVar = double()
  for(i in 1:length(unqCountry)){
      allValues = unlist(subset(tmp, select = origVar,
        subset = FAOST_CODE == unqCountry[i]))
      baseYearValue = unlist(subset(tmp, select = origVar,
        subset = FAOST_CODE == unqCountry[i] &
        Year == baseYear))
      ## Make the index NA if the base year value is not available
      if(length(baseYearValue) == 0)
        baseYearValue = NA
      tmpind = allValues/baseYearValue * 100
    indVar = c(indVar, tmpind)
  }
  ind.df = data.frame(tmp[c("FAOST_CODE", "Year")], indVar)
  colnames(ind.df)[3] = ifelse(is.na(newVarName),
                                paste(origVar, ".IND", sep = ""),
                                newVarName)
  ind.df
}

utils::globalVariables(names = "FAOST_CODE")
