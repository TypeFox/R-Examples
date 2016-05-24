if (getRversion() >= "2.15.1") globalVariables(c("Units"))
#' @encoding UTF-8
#' @title Unit Converter
#'
#' @description Converts data from a numeric value from one measurement system to another. For instance, distances in miles to kilometers.
#'
#' @param x A numeric value or vector of data values to be converted.
#' @param from A character defining the original unit.
#' @param to A character defining the target unit.
#'
#' @details \code{NA} is returned if a conversion cannot be found.
#' \tabular{lll}{
#'  \cr\cr
#'  \bold{Weight and mass}\tab \tab\cr
#'  Gram  \tab g \tab metric\cr
#'  Slug  \tab sg \cr
#'  Pound mass (avoirdupois)  \tab lbm \cr
#'  U (atomic mass unit)  \tab u \cr
#'  Ounce mass (avoirdupois)  \tab ozm \cr
#'  \cr\bold{Distance}\tab \cr
#'  Meter  \tab  m \tab metric \cr
#'  Statute mile  \tab mi \cr
#'  Nautical mile  \tab Nmi \cr
#'  Inch  \tab in \cr
#'  Foot  \tab ft \cr
#'  Yard  \tab yd \cr
#'  Angstrom  \tab ang \tab metric \cr
#'  Pica  \tab pica \cr
#'  \tab  \cr
#'  \cr\bold{Time}\tab \cr
#'  Year  \tab yr \cr
#'  Day  \tab day \cr
#'  Hour  \tab hr \cr
#'  Minute  \tab mn \cr
#'  Second  \tab sec \cr
#'  \cr\bold{Pressure}\tab \cr
#'  Pascal  \tab Pa (or p) \cr
#'  Atmosphere  \tab atm (or at) \cr
#'  mm of Mercury  \tab mmHg \cr
#'  \tab  \cr
#'  \cr\bold{Force}\tab \cr
#'  Newton  \tab N \tab metric \cr
#'  Dyne  \tab dyn (or dy) \cr
#'  Pound force  \tab lbf \cr
#'  \cr\bold{Energy}\tab \cr
#'  Joule  \tab J \tab metric \cr
#'  Erg  \tab e \cr
#'  Thermodynamic calorie  \tab c \cr
#'  IT calorie  \tab cal \tab metric \cr
#'  Electron volt  \tab eV (or ev) \tab metric \cr
#'  Horsepower-hour  \tab HPh (or hh) \cr
#'  Watt-hour  \tab Wh (or wh) \tab metric \cr
#'  Foot-pound  \tab flb \cr
#'  BTU  \tab BTU (or btu) \cr
#'  \tab  \cr
#'  \cr\bold{Power}\tab \cr
#'  Horsepower  \tab HP (or h) \cr
#'  Watt  \tab W (or w) \tab metric \cr
#'  \cr\bold{Magnetism}\tab \cr
#'  Tesla  \tab T \tab metric \cr
#'  Gauss  \tab ga \tab metric \cr
#'  \tab  \cr
#'  \cr\bold{Temperature}\tab \cr
#'  Degree Celsius  \tab C (or cel) \cr
#'  Degree Fahrenheit  \tab F (or fah) \cr
#'  Kelvin  \tab K (or kel) \tab metric \cr
#'  \cr\bold{Liquid measure}\tab \cr
#'  Teaspoon  \tab tsp \cr
#'  Tablespoon  \tab tbs \cr
#'  Fluid ounce  \tab oz \cr
#'  Cup  \tab cup \cr
#'  U.S. pint  \tab pt (or us_pt) \cr
#'  U.K. pint  \tab uk_pt \cr
#'  Quart  \tab qt \cr
#'  Gallon  \tab gal \cr
#'  Liter  \tab l (or lt) \tab metric \cr
#'}
#' @examples
#' converter(c(5.6, 6.7), "in", "m")
#'
#' @export
`converter` <- function(x, from, to){
  if(from == "C") {
    if(to=="F") return(x *1.8+32)
  }
  if(from == "F") {
    if(to=="C") return((x -32) *5/9)
  }

  factor <- Units[Units$from == from & Units$to==to, "factor"]
  if(length(factor)==0) factor <- NA

  return(x * factor)

}
NULL
