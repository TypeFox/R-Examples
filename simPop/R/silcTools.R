################################################################
### utility functions for the sga project on EU-SILC
#' @name silcTools
#' @rdname silcTools
#' @aliases getAge getGender getEcoStat getCitizenship getHsize restructureHHid factorNA
#' @title Utility functions for socio-economic data sets
#' @description Various utility functions mainly used for simulating EU-SILC data
#' @param birth year of birth
#' @param year current year
#' @param data data.frame (for some functions optional)
#' @author Andreas Alfons, Matthias Templ
NULL

#' @rdname silcTools
#' @name getAge
#' @examples
#' birth <- sample(1950:2000, 20)
#' getAge(birth, 2013)
#' @export
getAge <- function(birth, year, data = NULL) {
  # this function is only applicable if the income
  # reference period is the previous calender year
  if(!is.null(data)) {
    if(missing(birth)) birth <- "rb080"
    birth <- data[, birth]
    if(missing(year)) year <- "rb010"
    year <- data[, year]
  }
  if(is.factor(year)) year <- as.numeric(as.character(year))
  if(is.factor(birth)) year <- as.numeric(as.character(birth))
  year - 1 - birth
}

#' @rdname silcTools
#' @name getGender
#' @param gender variable including information on gender
#' @param labels labels of a factor variable
#' @export
#' @keywords internal
#' @examples
#' data(eusilcS)
#' head(getGender("rb090", labels = c("ma","fe"), data=eusilcS))
getGender <- function(gender, labels = c("male","female"), data = NULL) {
  if(!is.null(data)) {
    if(missing(gender)) gender <- "rb090"
    gender <- c(data[, gender])
  }
  if(is.list(gender)) gender <- gender[[1]]
  factor(c(gender), labels=labels)
}


#getHsize <- function(data)
#{
#  tab <- table(data$rb040)
#  hsize <- rep(tab, tab)
#}

#' @rdname silcTools
#' @name getEcoStat
#' @param ecoStat variable holding information on the economic status
#' @examples
#' lev <- c("Employee working full-time", "Employee working part-time",
#'          "Self-employed working full-time", "Self-employed working part-time",
#'          "Unemployed", "Pupil, student, further training, unpaid work experience",
#'          "In retirement", "Permanently disabled", "In compulsory military community or service",
#'          "Fulfilling domestic tasks", "Other inactive person")
#' g <- getEcoStat("pl030", eusilcS, lev)
#' table(g)
#' @export
getEcoStat <- function(ecoStat , data , levels) {  ## variable pl030 (economic status)
  if(missing(ecoStat)) ecoStat <- "pl031"
  ecoStat <- factor(data[, ecoStat])
  levels(ecoStat) <- levels
  return(ecoStat)
}

#' @rdname silcTools
#' @name getCitizenship
#' @examples
#' data(eusilcS)
#' ## destroy info on pb220a to show afterwards the usage of the function
#' owncountry <- "AT"
#' EU <- c("DE","BE","BG","CY","CZ","DK","EE","EL","ES","FI","FR","GR","HU","IE",
#'         "IT","LT","LU","LV","MT","NL","PL","PT","RO","SI","SE","SK","UK")
#' other <- c("CAN","CH","CSA","HR","IS","ME","MK","NAF","NME","NO",
#'            "OAF","OAS","OCE","OEU","OT","OTH","TR","USA","WAF")
#' eusilcS$fakepb220a <-  factor(sample(c(owncountry, EU, other), nrow(eusilcS), replace = TRUE))
#' table(eusilcS$fakepb220a)
#' eusilcS$fakepb220a <- getCitizenship(citizenship = "fakepb220a",
#'                                   data=eusilcS, owncountry=owncountry,
#'                                   EU=EU, other=other)
#' table(eusilcS$fakepb220a)
#' @export
getCitizenship <- function(citizenship, data, owncountry, EU, other) {
  if(missing(citizenship)) citizenship <- "pb220a"
  citizenship <- data[, citizenship]
  indNA <- which(levels(citizenship) == "")
  indOC <- which(levels(citizenship) == owncountry)
  indEU <- which(levels(citizenship) %in% EU)
  indOther <- which(levels(citizenship) %in% other)
  levels <- character(nlevels(citizenship))
  levels[indNA] <- NA
  levels[indOC] <- owncountry
  levels[indEU] <- "EU"
  levels[indOther] <- "Other"
  levels(citizenship) <- levels
  return(citizenship)
}


#' @rdname silcTools
#' @name getHsize
#' @param hhid name or index of variable holding the information on household ID
#' @examples
#' data(eusilcS)
#' hsize <- getHsize(data=eusilcS)
#' table(hsize)
#' @export
getHsize <- function(data, hhid)
{
  if(missing(hhid)) hhid <- "db030"
  tab <- table(data[,hhid]) #table(data$rb040)
  hsize <- rep(tab, tab)
  hsize <- as.numeric((hsize))
  return(hsize)
}

#' @rdname silcTools
#' @name restructureHHid
#' @examples
#' hhid <- c(6,6,3,3,3,2,1,1,8,9,9,9,9,7,7)
#' hhid
#' df <- data.frame("hhid"=hhid)
#' restructureHHid(df, "hhid")
#' @export
restructureHHid <- function(data, hhid){
  if(missing(hhid)) hhid <- "db030"
  tab <- table(data[, hhid])
  hsize <- rep(tab, tab)
  db030 <- as.numeric(names(hsize))
  return(db030)
}

#Function factorNA from package simPop: includes NAs as an extra level in the factor
#' @rdname silcTools
#' @name factorNA
#' @examples
#' hhid <- factor(c(6,6,3,3,3,2,1,1,NA,9,9,9,9,7,7))
#' hhid
#' factorNA(hhid)
#' @export
factorNA <- function(x, always = FALSE, newval = NA) {
  always <- isTRUE(always)
  if(is.na(newval)){
    if(is.factor(x)) {
      l <- levels(x)
      if(NA %in% l || !(always || any(is.na(x)))) x
      else {
        l <- c(l, NA)
        factor(x, levels=c(levels(x), NA), exclude=c())
      }
    } else {
      if(always) {
        factor(c(NA, x), exclude=c())[-1] # little trick
      } else factor(x, exclude=c())
    }
  } else { # newval is a character
    if(!is.character(newval)) stop("newval must be NA or a character string")
    if(is.factor(x)) {
      l <- levels(x)
      l <- c(l, newval)
      x <- as.character(x)
      x[is.na(x)] <- newval
      factor(x)
    }
  }
}


#
# # Function uni.distribution: random draws from the weighted univariate distribution of
# # the original data (maybe better from the SUF, but then the SUF always has to be used as well)
# univariate.dis <- function(puf,data,additional,w){
#   if (sum(is.na(data[,additional]))>0 & sum(is.na(data[,additional])) != dim(data)[1]) {
#     var <- factorNA(data[,additional],always=TRUE)
#   } else if (sum(is.na(data[,additional])) == dim(data)[1]) {
#     var <- factor(c(NA, data[,additional]), exclude=c())[-1]
#   } else {
#     var <- as.factor(data[,additional])
#   }
#   tab <- wtd.table(var,weights=data[,w],type="table")
#   p <- tab/sum(data[,w])
#   puf[,additional] <- sample(x=levels(var)[levels(var) %in% names(tab)],size=dim(puf)[1],prob=p,replace=T)
#   return(puf)
# }
#
# # Function con.distribution: random draws from the weighted conditional distribution
# # (conditioned on a factor variable)
# conditional.dis <- function(puf,data,additional,conditional,w){
#   if (sum(is.na(data[,additional]))>0 & sum(is.na(data[,additional])) != dim(data)[1]) {
#     var <- factorNA(data[,additional],always=TRUE)
#   } else if (sum(is.na(data[,additional])) == dim(data)[1]) {
#     var <- factor(c(NA, data[,additional]), exclude=c())[-1]
#   } else {
#     var <- as.factor(data[,additional])
#   }
#   puf[,additional] <- NA
#   for (i in 1:length(levels(puf[,conditional]))) {
#     tab <- wtd.table(var[data[,conditional]==levels(data[,conditional])[i]],weights=data[data[,conditional]==levels(data[,conditional])[i],w],type="table")
#     p <- tab/sum(tab)
#     puf[which(puf[,conditional]==levels(puf[,conditional])[i]),additional] <- sample(x=levels(var)[levels(var) %in% names(tab)],size=dim(puf[which(puf[,conditional]==levels(data[,conditional])[i]),])[1],prob=p,replace=T)
#   }
#   return(puf)
# }

