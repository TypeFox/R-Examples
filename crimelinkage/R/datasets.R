#' Ficticious dataset of crime events
#'
#' Some realistic, but fictious, crime incident data.
#' @usage data(crimes)
#' @format 490 crime events
#' \describe{
#'    \item{crimeID}{The crime ID number}
#'    \item{X, Y}{Spatial coordinates}
#'    \item{MO1}{A categorical MO variable that takes values {1,\dots,31}}
#'    \item{MO2}{A categorical MO variable that takes values {a,\dots,h}}
#'    \item{MO3}{A categorical MO variable that takes values {A,\dots,O}}
#'    \item{DT.FROM}{The earliest possible Date-time of the crime.}
#'    \item{DT.TO}{The latest possible Date-time of the crime}
#'  }
#' @examples
#' head(crimes)
#' @source Ficticious data, but hopefully realistic
 "crimes"


#' Ficticious offender data
#'
#' Offender table relating crimes (\code{crimeID}) to offenders (\code{offenderID})
#' @usage data(offenders)
#' @format 1357 offenders committed 1377 crimes
#' \describe{
#'  \item{offenderID}{ID number of offender}
#'  \item{crimeID}{ID number of crime}
#' }
#' @examples
#' head(offenders)
#' @source Ficticious data, but hopefully realistic
 "offenders"
