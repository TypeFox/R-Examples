#' Simulated Metal Contents
#'
#' A data set with two time-series of simulated metal contents.  These data
#' have mainly been used to test the package. Users may find this dataset
#' convenient as an example to construct their own data sets.
#' The columns represent the following information:
#'  \itemize{
#'      \item{OBJECTID: water body code, e.g., NL89_os;}
#'      \item{PAR: parameter, e.g., Cadmium;}
#'      \item{color: colors in density function;}
#'      \item{DATE: date according to ISO 8601 (YYYY-mm-dd) for point values
#'          or year YYYY for annual means;}
#'      \item{VALUE: numerical value.}
#'      \item{TARGET: e.g., the target value for the European Water Framework 
#'          Directive;}
#'      \item{UNIT: measurement unit of PAR. This unit should be the same for all
#'          records with the same PAR and is the same for both VALUE 
#'          and TARGET;}
#'      \item{transfrom: data transformation, i.e., log, logit, NA.}
#'  }
#'
#'  @docType data
#'  @keywords datasets
#'  @name metal
#'  @usage data(metal)
NULL






#' Annual Average 1,2-dichloroethane Concentration
#'
#' Annual arithmetic average concentration of 1,2-dichloroethane (DCA) in
#' a specific water body (\eqn{\mu g/l}{ug/l}), based on Baggelaar et al., (2010)
#'
#' The columns represent the following information:
#'  \itemize{
#'      \item{OBJECTID: water body code}
#'      \item{PAR: parameter, in this case 1,2-dichloroethane;}
#'      \item{color: colors in density function;}
#'      \item{DATE: year;}
#'      \item{VALUE: annual arithmetic average concentration}
#'      \item{TARGET: target according to the European Water Framework Directive;}
#'      \item{UNIT: measurement unit (\eqn{\mu g/l}{ug/l}).}
#'  }
#'
#'  @docType data
#'  @keywords datasets
#'  @name DCA
#'  @references Baggelaar, P., O. van Tongeren, R. Knoben, & W. van Loon, 2010. 
#'      Rapporteren van de betrouwbaarheid van KRW-beoordelingen. H2O 16, p.21--25
#'  @usage data(DCA)
NULL




#' Annual Average Environmental Quality Ratio for Macrofauna.
#'
#' Annual artithmetic average environmental quality ratio's (EQR) for Macrofauna
#' in a specific water body, based on Baggelaar et al., (2010)
#'
#' The columns represent the following information:
#'  \itemize{
#'      \item{OBJECTID: water body code}
#'      \item{PAR: parameter, in this case EQR;}
#'      \item{color: colors in density function;}
#'      \item{DATE: year;}
#'      \item{VALUE: annual arithmetic average EQR}
#'      \item{TARGET: target EQR;}
#'      \item{transfrom: applied transform}
#'  }
#'
#'  @docType data
#'  @keywords datasets
#'  @name EQR
#'  @references Baggelaar, P., O. van Tongeren, R. Knoben, & W. van Loon, 2010. 
#'      Rapporteren van de betrouwbaarheid van KRW-beoordelingen. H2O 16, p.21--25
#'  @usage data(EQR)
NULL