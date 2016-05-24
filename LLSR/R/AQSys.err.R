# Each number used as argument triggers an error with its respective message
AQSys.err <- function (err, ...) {
  switch(
    err,
    "0" = {
      errmsg <- "The selected Equation doesn't exist."
      stop(errmsg, call. = FALSE)
    },
    "1" = {
      errmsg <-
        "The path must point to a '.xls' or '.xlsx' Microsoft Excel worksheet."
      stop(errmsg, call. = FALSE)
    },
    "2" = {
      errmsg <-
        "There was a problem when calculating the number of systems to be analysed. Please make sure the data is formatted as the standard provided."
      stop(errmsg, call. = FALSE)
    },
    "3" = {
      errmsg <- " must be a data.frame."
      stop(..., errmsg, call. = FALSE)
    },
    "4" = {
      errmsg <- "Argument 'db' is missing. Parameter can not be NULL."
      stop(errmsg, call. = FALSE)
    },
    "5" = {
      errmsg <-
        "Your search had no results. Try removing a few parameters during your next search."
      stop(errmsg, call. = FALSE)
    },
    "6" = {
      errmsg <-
        "At least one of the parameters (pH, Temp, additive, UP.Rich or LP.Rich) must be not NULL."
      stop(errmsg, call. = FALSE)
    },
    "7" = {
      errmsg <-
        "Input variable db must be a list coontaining three data.frame variables (db.cas, db.ref and db.sys)."
      stop(errmsg, call. = FALSE)
    },
    "8" = {
      errmsg <- "Your search had no results."
      stop(errmsg, call. = FALSE)
    },
    {
      errmsg <- "An Unknown error ocourred."
      stop(errmsg, call. = FALSE)
    }
  )
  
}
