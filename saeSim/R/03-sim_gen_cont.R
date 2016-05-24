#' Generation Component for contamination
#' 
#' One of the components which can be added to a \code{sim_setup}. It is applied after functions added with \code{\link{sim_gen}}.
#' 
#' @inheritParams sim_agg
#' @inheritParams sim_gen
#' 
#' @param nCont gives the number of contaminated observations. Values between 0 
#'   and 1 will be treated as probability. If type is 'unit' and length is
#'   larger than 1, the expected length is the number of areas. If type is
#'   'area' and length is larger than 1 the values are interpreted as area
#'   positions; i.e. \code{c(1, 3)} is interpreted as the first and 3rd area in
#'   the data is contaminated.
#' @param type "unit" or "area" - unit- or area-level contamination.
#' @param areaVar character with variable name(s) identifying areas.
#' @param fixed TRUE fixes the observations which will be contaminated. FALSE
#'   will result in a random selection of observations or areas.
#' 
#' @seealso \code{\link{sim_gen}}
#' @export
#' 
#' @examples
#' sim_base_lm() %>% 
#'   sim_gen_cont(gen_norm(name = "e"), nCont = 0.05, type = "unit", areaVar = "idD") %>%
#'   as.data.frame
sim_gen_cont <- function(simSetup, generator, nCont, type, areaVar = NULL, fixed = TRUE) {
  
  generator <- gen_cont(generator, nCont, type, areaVar, fixed)
  
  sim_setup(
    simSetup, 
    new("sim_fun", order = 2, call = match.call(), generator)
  )
}

gen_cont <- function(generator, nCont, type, areaVar, fixed) {
  force(generator); force(nCont); force(type); force(areaVar); force(fixed)
  check_cont_input(nCont, type, fixed)
  
  genFun <- function(dat) {
    contData <- generator(dat)
    contData <- select_cont(contData, nCont, type, areaVar, fixed)
    replace_contData(contData, dat)
  }
  
  preserve_attributes(genFun)
  
}

check_cont_input <- function(nCont, type, fixed) {
  
  if (length(nCont) == 1 && nCont == 1) {
    warning("nCont is equal to 1 - one area is contaminated.")
  }
  
  if (!(type %in% c("area", "unit"))) {
    stop("Supported types are area and unit!")
  }
  
  if (length(nCont) > 1 & (type == "area")) {
    message("A vector of nCont with type 'area' is interpreted as vector with positions!")
    if (!fixed) stop("With this settings, fixed should be 'TRUE'!")
    stopifnot(all(nCont %% 1 == 0))
  }
}

replace_contData <- function(contData, dat) {
  vars <- names(contData)[names(contData) %in% names(dat)]
  for (var in vars) contData[var] <- replace_cont(contData[[var]], dat[[var]])
  contData
}

replace_cont <- function(var1, var2) ifelse(var1 == 0, var2, var1)


