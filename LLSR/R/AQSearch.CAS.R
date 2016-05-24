####################################################################################################################
#' @rdname AQSearch.CAS
#' @name AQSearch.CAS
#' @description Search CAS codes in LLSR's Chemical database.
#' @title Search CAS codes in Chemical database.
#' @method AQSearch CAS
#' @export AQSearch.CAS
#' @export
#' @param db A list of three data.frame variables within all systems data are stored.
#' @param ChemString A String with a text fragment of the chemical name the user want to search the correspondent CAS.
#' @param ... Additional optional arguments. None are used at present.
#' @examples
#' \dontrun{
#' AQSearch.CAS(db, "NaOH")
#'}
####################################################################################################################
AQSearch.CAS <- function(db, ChemString, ...) {
  #
  # names(db) <- c("CAS.CODE", "CHEM.NAME", "CHEM.COMMON")
  # Checks if variable db is corretly formatted. See AQSearchUtils.R for more details.
  db.check(db)
  # Initialization of ans, which will be used to return results.
  ans <- data.frame()
  # Check if search parameter is not null.
  if (is.null(ChemString) == FALSE) {
    # convert both input String and db string to uppercase and compare them.
    # function fetch data from chem name and common name.
    db.grep <-
      unique(c(grep(
        toupper(ChemString), toupper(db$db.cas[,2])
      ),
      grep(
        toupper(ChemString), toupper(db$db.cas[,3])
      )))
    # check if search returned any result
    if ((length(db.grep) != 0) & (db.grep[1] != 0)) {
      # if yes, store in the output variable
      ans <- db$db.cas[db.grep,]
    }else{
      # if no result was found, it triggers an error (check AQSys.err.R for details)
      AQSys.err("8")
    }
  } else {
    # if string isn't valid, it triggers an error (check AQSys.err.R for details)
    AQSys.err("4")
  }
  # return search results
  invisible(ans)
}
