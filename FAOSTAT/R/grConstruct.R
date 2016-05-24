##' Construct Growth rate
##'
##' A function for constructing growth rate variables.
##'
##' @param data The data frame containing the data
##' @param origVar The variable in which the growth is to be calculated
##' @param newVarName The name assigned to the new variable, if missing
##' then .SC/.SH/.GR will be appended depending on the type of
##' construction.
##' @param type The type of growth rate, can be least squares or
##' geometric
##' @param n The period for the growth rate to be calculated (Refer to
##' the lsgr or the geogr functions.)
##' @return A data frame containing the computed growth rate.
##' @export
##' @examples
##' test.df2 = data.frame(FAOST_CODE = rep(c(1, 5000), each = 5),
##'                        Year = rep(1990:1994, 2),
##'                        a = rep(1:5, 2), b = rep(1:5, 2))
##' grConstruct(test.df2, origVar = "a", type = "geo", n = 1)
##' grConstruct(test.df2, origVar = "a", type = "geo", n = 3)
##' grConstruct(test.df2, origVar = "a", type = "geo", n = 5)
grConstruct = function(data, origVar, newVarName = NA,
                        type = c("geo", "ls", "ch"), n = 1){
  type = match.arg(type)
  tmp = arrange(subset(data, select = c("FAOST_CODE", "Year", origVar)),
                 FAOST_CODE, Year)
  unqCountry = unique(tmp$FAOST_CODE)
  grVar = double()
  if(type == "ls"){
    for(i in 1:length(unqCountry)){
      tmpgr =
        lsgr(as.numeric(unlist(subset(tmp,
                                      select = origVar,
                                      subset = FAOST_CODE == unqCountry[i]))), n)
      grVar = c(grVar, tmpgr)
    }
  } else if(type == "geo"){
    for(i in 1:length(unqCountry)){
      tmpgr =
        geogr(as.numeric(unlist(subset(tmp,
                                       select = origVar,
                                       subset = FAOST_CODE == unqCountry[i]))),
              n)
      grVar = c(grVar, tmpgr)
    }
  } else if(type == "ch"){
    for(i in 1:length(unqCountry)){
      tmpgr =
        chgr(as.numeric(unlist(subset(tmp,
                                      select = origVar,
                                      subset = FAOST_CODE == unqCountry[i]))),
             n)
      grVar = c(grVar, tmpgr)
    }
  }
  gr.df = data.frame(tmp[c("FAOST_CODE", "Year")], grVar)
  colnames(gr.df)[3] = ifelse(is.na(newVarName),
                               paste(origVar, ".GR", n, sep = ""),
                               newVarName)
  gr.df
}

utils::globalVariables(names = "FAOST_CODE")
