#' @import rootSolve
require(rootSolve)
#' @rdname AQSysOthmer
#' @title Othmer's Equation - Tieline's correlation
#' @description Othmer's equation to correlate tieline's data applying the lever's rule.
#' @references Othmer, D.F. and P.E. Tobias, Liquid -Liquid Extraction Data -Toluene and
#' Acetaldehyde Systems. Industrial & Engineering Chemistry, 1942. 34(6): p. 690-692.
#' @param ... Additional optional arguments. None are used at present.
#' @param TLdt - Tieline Experimental data that will be used in the nonlinear fit
#' @export AQSysOthmer
#' @return Parameters K, n and Statistical data
#' @examples
#' # TLdt is a data.frame which contains series of Tieline's mass fraction
#' # (upper-rich component, bottom-rich component and water)
#' # Each column in the data.frame represents a series of one component mass fraction
#' # For example, an empty data.frame for four tielines can be obtaining using:
#' TLdt<-matrix(NA,nrow=4,ncol=6)
#' # Variables order must follows the sequence presented below:
#' # "mfXt","mfYt","mfXb","mfYb","mfWt","mfWb"
#' # In which: mf stands for mass fraction; X and Y for the component
#' # rich in bottom and upper phase, respectively; t or b for top and
#' # bottom phases and W for water.
#' # Then you just need to load the data.frame in the function:
#' \dontrun{
#' AQSysOthmer(TLdt)
#'}
AQSysOthmer <- function(TLdt,...) {
  # store tieline data into a dataframe variable. It might be a better approach check if
  # user stored it in a dataframe and if not trigger an error.
  TLdt <- as.data.frame(TLdt)
  # tieline data is a set of mass fractions of all systems components obtained
  # experimentally for the system's upper and bottom phase.
  # the line bellow set the dataset header
  names(TLdt) <- c("mfXt","mfYt","mfXb","mfYb","mfWt","mfWb")
  # the system below will calculate n and K for a given set of tielines
  suppressWarnings(
    FFn <- nls(
      log((1 - mfYt) / mfYt) ~ log(K * (((
        1 - mfXb
      ) / mfXb)) ^ n),
      start = list(n = 1,K = 1),
      algorithm = "port",
      lower = 10 ^ -10,
      data = TLdt,na.exclude
    )
  )
  # return all calculated parameters
  FFn
}
#' @name AQSysBancroft
#' @title Bancroft's Potential Equation - Tieline's correlation
#' @description Bancroft's equation to correlate tieline's data.
#' @references Othmer, D.F. and P.E. Tobias, Liquid-Liquid Extraction Data -Toluene and
#' Acetaldehyde Systems. Industrial & Engineering Chemistry, 1942. 34(6): p. 690-692.
#' @export AQSysBancroft
#' @param ... Additional optional arguments. None are used at present.
#' @param TLdt - Tieline Experimental data that will be used in the nonlinear fit
#' @return Parameters K1, r and Statistical data
#' @examples
#' # TLdt is a data.frame which contains series of Tieline's mass fraction
#' # (upper-rich component, bottom-rich component and water)
#' # Each column in the data.frame represents a series of one component mass fraction
#' # For example, an empty data.frame for four tielines can be obtaining using:
#' TLdt<-matrix(NA,nrow=4,ncol=6)
#' # Variables order must follows the sequence presented below:
#' # "mfXt","mfYt","mfXb","mfYb","mfWt","mfWb"
#' # In which: mf stands for mass fraction; X and Y for the component
#' # rich in bottom and upper phase, respectively; t or b for top and
#' # bottom phases and W for water.
#' # Then you just need to load the data.frame in the function:
#' \dontrun{
#' AQSysBancroft(TLdt)
#'}
AQSysBancroft <- function(TLdt,...) {
  # store tieline data into a dataframe variable. It might be a better approach check if
  # user stored it in a dataframe and if not trigger an error.
  TLdt <- as.data.frame(TLdt)
  # tieline data is a set of mass fractions of all systems components obtained
  # experimentally for the system's upper and bottom phase.
  # the line bellow set the dataset header
  names(TLdt) <- c("mfXt","mfYt","mfXb","mfYb","mfWt","mfWb")
  # the system below will calculate r and K1 for a given set of tielines
  #suppressWarnings(
  FFn <- nls(
    log((mfWb / mfXb)) ~ log(K1 * ((mfWt / mfYt) ^ r)),
    start = list(r = 1,K1 = 1),
    algorithm = "port",
    lower = 10 ^ -2,
    data = TLdt,na.exclude
  )
  #)
  # return all calculated parameters
  FFn
}
