##' Construct share variable
##' 
##' A function for constructing the share of a variable of an aggregated
##' variable.
##'
##' The share of a variable can be share of the World (if additional
##' variable were not supplied) or share of another variable (per Capita
##' if population was supplied).
##'
##' @param data The data frame containing both the share variable and
##' the aggregated variable
##' @param totVar The aggregated variable.
##' @param shareVar The subset of the aggregated variable which to be
##' divided by.
##' @param newVarName The name assigned to the new variable, if missing
##' then .SC/.SH/.GR will be appended depending on the type of
##' construction
##' @return A data frame with the new constructed variable
##' @export
##' @examples
##' ## Total variables provided, scale by totVar
##' test.df = data.frame(FAOST_CODE = 1, Year = 1990:1994, a = 1:5, b = 1:5)
##' shConstruct(data = test.df, totVar = "a", shareVar = "b")
##'
##' ## Total variables not provided, scale by world aggregate.
##' test.df2 = data.frame(FAOST_CODE = rep(c(1, 5000), each = 5),
##'                        Year = rep(1990:1994, 2),
##'                        a = rep(1:5, 2), b = rep(1:5, 2))
##' shConstruct(data = test.df2, totVar = NA, shareVar = "b")
shConstruct = function(data, totVar, shareVar, newVarName = NA){
  data = arrange(data, FAOST_CODE, Year)
  if(is.na(totVar)){
    shareData = subset(data, select = c("FAOST_CODE", "Year", shareVar))
    ## world.df = subset(data, subset = FAOST_CODE == 5000,
    ##                    select = c("Year", shareVar))
    ## colnames(world.df) = c("Year", "World")
    ## data = arrange(merge(data, world.df), FAOST_CODE, Year)
    ## totVar = "World"
    
    world.df = ddply(.data = shareData, .variables = .(Year),
      .fun = function(x) sum(x[[shareVar]], na.rm = TRUE))
    colnames(world.df) = c("Year", "World")
    fd = merge(shareData, world.df, by = "Year")
    newVarName = ifelse(is.na(newVarName), paste(shareVar, ".SC", sep = ""),
      newVarName)
    fd[[newVarName]] = fd[[shareVar]]/fd[["World"]]
    shVar.df = subset(fd, select = c("FAOST_CODE", "Year", newVarName))
    
    ## shVar = data[shareVar]/data[totVar]
    ## colnames(shVar) = ifelse(is.na(newVarName),
    ##                           paste(shareVar, ".SC", sep = ""),
    ##                           newVarName)
  } else {
    ## I don't think we should multiply by 100
    shareData = subset(data, select = c("FAOST_CODE", "Year", shareVar, totVar))
    newVarName = ifelse(is.na(newVarName), paste(shareVar, ".SH", sep = ""),
      newVarName)    
    shareData[[newVarName]] = data[[shareVar]]/data[[totVar]]
    shVar.df = subset(shareData, select = c("FAOST_CODE", "Year", newVarName))

  }
  arrange(shVar.df, FAOST_CODE, Year)
  ## merge(data, shVar.df, by = c("FAOST_CODE", "Year"))
  ## data.frame(data[c("FAOST_CODE", "Year")], shVar)
}

utils::globalVariables(names = "FAOST_CODE")
