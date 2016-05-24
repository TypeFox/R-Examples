####################################################################################################################
#' @rdname AQSearch
#' @name AQSearch
#' @description This function allow the user to search the LLSR database to find any ATPS that matches the used criteria.
#' @export
####################################################################################################################
AQSearch <- function(db, ...)
  UseMethod("AQSearch")
####################################################################################################################
#' @rdname AQSearch
#' @title Search function for ATPS Systems data
#' @description This function allow the user to search the package database to find any ATPS that matches the used criteria.
#' @details The function return the systems that matches the criteria submit by the user. Beware that using it with subset = FALSE will return the system's entry reference as a MD5 encoded string.
#' @param db A list of three data.frame variables within all systems data are stored.
#' @param db.ph A numeric variable containing the pH to be searched within DB.
#' @param db.upper A String variable containing the CAS of the lower's phase component to be searched within DB.
#' @param db.lower A String variable containing the CAS of the lower's phase component to be searched within DB.
#' @param db.temp A numeric variable containing the Temperature (in Kelvin) to be searched within DB.
#' @param db.additive A String variable containing the CAS of the additive component to be searched within DB.
#' @param subset If FALSE, the function returns the index of the systems that matches the search. If TRUE, the function returns a subset with all system data.
#' @param ... Additional optional arguments. None are used at present.
#' @method AQSearch default
#' @export
#' @return System's index (when subset = FALSE) or system's parameters (subset = TRUE)
#' @examples
#' \dontrun{
#' AQSearch(db,db.upper="85100-78-3", subset = TRUE)
#'}
####################################################################################################################
AQSearch.default <-
  function(db = NULL, db.ph = NULL, db.upper = NULL, db.lower = NULL,
           db.temp = NULL, db.additive = NULL, subset = FALSE, ...) {
    #
    #inialize db.ans
    db.ans <- data.frame()
    # create and initialise a list using the function's parameters
    db.params <- c(db.ph, db.upper, db.lower, db.temp, db.additive)
    # Checks if variable db is corretly formatted. See AQSearchUtils.R for more details.
    db.check(db)
    # if all parameters are null, the search is not valid and it triggers an error (check AQSys.err.R for details)
    if (all(unlist(lapply(db.params, is.null))))
      AQSys.err("6")
    # output variable is initialised with data from db.
    db.grep <- db$db.sys
    # each conditional will return a result set that excludes no match data.
    # search a system that matchs the upper-phase component, if search parameter is not null.
    if (is.null(db.upper) == FALSE) {
      db.grep <- db.grep[grep(db.upper, db.grep[,2]),]
    }
    # search a system that matchs the lower-phase component, if search parameter is not null.
    if (is.null(db.lower) == FALSE) {
      db.grep <- db.grep[grep(db.lower, db.grep[,3]),]
    }
    # search a system that matchs the system's pH, if search parameter is not null.
    if (is.null(db.ph) == FALSE) {
      db.grep <- db.grep[grep(db.ph, db.grep[,4]),]
    }
    # search a system that matchs the additive component, if search parameter is not null.
    if (is.null(db.additive) == FALSE) {
      db.grep <- db.grep[grep(db.additive, db.grep[,5]),]
    }
    # search a system that matchs the system's temperature, if search parameter is not null.
    if (is.null(db.temp) == FALSE) {
      db.grep <- db.grep[grep(db.temp, db.grep[,7]),]
    }
    # If search isn't null, evaluate data
    if (length(rownames(db.grep)) != 0) {
      # If subset parameter is true, return a db with the matched data
      if (subset == TRUE) {
        # use result's indexes to populate output datalist
        db.ans <- db.grep[rownames(db.grep),]
        for (k in 1:nrow(db.ans)) {
          # undigest reference entry
          db.ans[k,1] <- undigest(db, db.ans[k,1])
        }
        #invisible(db.grep[rownames(db.grep),])
        #Add bibliographic reference of system's data
        names(db.ans)[1] <- "Bibliographyc Reference"
        invisible(db.ans)
      }
      # otherwise return indexes of matched data
      else
        invisible(rownames(db.grep))
      #
    }else{
      # Triggers an "no results" error
      AQSys.err("5")
    }
    #
  }
