


#' @title Expected / Population Hazard Data Sets Usage in \pkg{popEpi}
#' @author Joonas Miettinen
#' @name pophaz
#' @description 
#' 
#' Several functions in \pkg{popEpi} make use of population or expected
#' hazards in computing the intended estimates (e.g. \code{\link{survtab}}).
#' This document explains using such data sets in this package.
#' 
#' @details 
#' 
#' Population hazard data sets (pophaz for short) in \pkg{popEpi} should
#' be \code{data.frame}s in the "long" format where one of the columns must be
#' named \code{haz} (for hazard), and other columns define the values or 
#' levels in variables relating to subjects in your data. For example,
#' \code{\link{popmort}} contains Finnish population mortality hazards
#' by sex, calendar year, and 1-year age group.
#' 
#' \tabular{rrrr}{
#' sex \tab year \tab agegroup \tab haz \cr
#' 0 \tab 1951 \tab 0 \tab 0.036363176\cr  
#' 0 \tab 1951 \tab 1 \tab 0.003616547\cr
#' 0 \tab 1951 \tab 2 \tab 0.002172384\cr
#' 0 \tab 1951 \tab 3 \tab 0.001581249\cr
#' 0 \tab 1951 \tab 4 \tab 0.001180690\cr
#' 0 \tab 1951 \tab 5 \tab 0.001070595
#' }
#' 
#' The names of the columns should match to the names of the variables
#' that you have in your subject-level data. Time variables in your pophaz
#' may also correspond to \code{Lexis} time scales; see 
#' \code{\link{survtab}}.
#' 
#' Any time variables (as they usually have) should be coded consistently:
#' When using fractional years in your data, the time variables in your pophaz
#' must also be coded in fractional years. When using e.g. \code{Date}s in your
#' data, ensure that the pophaz time variables are coded at the level of days
#' (or \code{Date}s for calendar time). 
#' 
#' The \code{haz} variable in your pophaz should also be coded consistently
#' with the used time variables. E.g. \code{haz} values in life-tables
#' reported as deaths per person-year should be multipleid by 365.25 when
#' using day-level time variables.
#' 
#' If you have your population hazards in a \code{ratetable} object
#' usable by functions in \pkg{survival} and \pkg{relsurv}, you may
#' transform them to long-format \code{data.frame}s using
#' \code{\link{as.data.frame.ratetable}}. Ensure, however, that the
#' created \code{haz} column is coded at the right level (events per
#' days or years typically).
#' 
#' National statistical institutions, the WHO, and e.g. the Human
#' Life-Table Database supply life-table data.
#' 

NULL









