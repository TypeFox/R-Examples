remind2magpie <- function(x) {
  #MAgPIE regions
  reg <- c("AFR.1", "CPA.2", "EUR.3", "FSU.4", "LAM.5", "MEA.6", "NAM.7", "PAO.8", "PAS.9", "SAS.10")
  if(!is.magpie(x)) stop("input is not a MAgPIE object!")
  x2 <- new.magpie(cells_and_regions=reg,years=getYears(x),names=getNames(x))
  x2["AFR",,] <- x["AFR",,]
  x2["CPA",,] <- x["CHN",,]
  x2["EUR",,] <- x["EUR",,]
  x2["FSU",,] <- x["RUS",,]
  x2["LAM",,] <- x["LAM",,]
  x2["MEA",,] <- x["MEA",,]
  x2["NAM",,] <- x["USA",,]
  x2["PAO",,] <- x["ROW",,]
  x2["PAS",,] <- x["OAS",,]
  x2["SAS",,] <- x["IND",,]
  return(x2) 
}