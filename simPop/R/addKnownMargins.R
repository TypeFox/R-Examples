#' add known margins/totals
#'
#' add known margins/totals for a combination of variables for the population
#' to an object of class \code{\linkS4class{simPopObj}}.
#'
#' @name addKnownMargins
#' @param inp a \code{simPopObj} containing population and household survey
#' data as well as optionally margins in standardized format.
#' @param margins a \code{data.frame} containing for a combination of unique
#' variable levels for n-variables the number of known occurences in the
#' population. The numbers must be listed in the last column of data.frame
#' 'margins' while the characteristics must be listed in the first 'n' columns.
#' @details The function takes a data.frame containing known marginals/totals for a some
#' variables that must exist in the population (stored in slot 'pop' of input
#' object 'inp') and updates slot 'table' of the input object. This slot
#' finally contains the known totals.
#'
#' households are drawn from the data and new ID's are generated for the new
#' households.
#' @return an object of class \code{\linkS4class{simPopObj}} with updated slot
#' 'table'.
#' @author Bernhard Meindl
#' @keywords manip
#' @export
#' @examples
#' data(eusilcS)
#' data(eusilcP)
#' \dontrun{
#' ## approx. 20 seconds computation time
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' inp <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#' inp <- simCategorical(inp, additional=c("pl030", "pb220a"), method="multinom",nr_cpus=1)
#'
#' margins <- as.data.frame(
#'   xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$region + eusilcP$gender + eusilcP$citizenship))
#' colnames(margins) <- c("db040", "rb090", "pb220a", "freq")
#' inp <- addKnownMargins(inp, margins)
#' str(inp)
#' }
addKnownMargins <- function(inp, margins) {
  dataP <- inp@pop@data
  margins <- as.data.frame(margins)
  if ( !class(margins) == "data.frame" ) {
    stop("input argument 'margins' must inherit class 'data.frame'!\n")
  }
  if ( any(duplicated(margins)) ) {
    stop("'margins' must not contain duplicated rows!\n")
  }
  if ( !class(margins[,ncol(margins)]) == "numeric" ) {
    stop("last column of input 'margins' must contain the numbers (must be numeric)!\n")
  } else {
    vals <- margins[,ncol(margins)]
    margins <- margins[,-ncol(margins), drop=FALSE]
  }
  if ( !class(inp) == "simPopObj" ) {
    stop("input argument 'inp' must be of class 'simPopObj'!\n")
  }
  if ( !all(colnames(margins) %in% colnames(dataP)) ) {
    stop("all variables existing in input 'margins' must also be existing in the
      synthetic population existing in slot 'pop' of input object 'inp'!\n")
  }

  # order: do all levels exist?
  ind <- match(colnames(margins), colnames(dataP))

  frame <- expand.grid(lapply(ind, function(x) {
    if ( is.factor(dataP[[x]]) ) {
      levels(dataP[[x]])
    } else {
      unique(dataP[[x]])
    }
  }))
  colnames(frame) <- colnames(margins)
  frame <- as.data.table(frame)

  margins$N <- vals
  margins <- as.data.table(margins)
  setkeyv(frame, colnames(frame))
  setkeyv(margins, colnames(margins))
  frame <- merge(frame, margins, all.x=TRUE)
  frame <- frame[,lapply(.SD, as.character)]
  frame$N <- as.numeric(frame$N)

  ind <- which(is.na(frame$N))
  if ( length(ind) > 0 ) {
    frame$N[ind] <- 0
  }

  if ( !is.null(inp@table) ) {
    message("Note: currently stored marginals/totals are going to be overwritten!\n")
  }
  inp@table <- frame
  invisible(inp)
}

