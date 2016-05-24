

#' @title Cast \code{data.table}/\code{data.frame} from long format to wide format
#' @author Matti Rantanen, Joonas Miettinen
#'
#' @description 
#' Convenience function for using \code{\link[data.table]{dcast.data.table}}
#' and \code{\link[reshape2]{dcast}};
#' inputs are character strings (names of variables) instead of a formula.
#' 
#' @param data a \code{data.table} or \code{data.frame}
#' @param columns a character string vector; the (unique combinations of the) 
#' levels of these variable will be different rows
#' @param rows a character string vector; the (unique combinations of the) 
#' levels of these variable will be different columns
#' @param values a character string; the variable which will be represented
#' on rows and columns as specified by \code{columns} and \code{rows}
#' @import data.table
#' @import stats
#' @export cast_simple
#' @details This function is just a small interface for \code{dcast} / 
#' \code{dcast.data.table} and less flexible than the originals.
#' 
#' Note that all \code{data.table} objects are also \code{data.frame} 
#' objects, but that each have their own \code{dcast} method.
#' \code{\link[data.table]{dcast.data.table}} is faster.
#' 
#' If any values in \code{value.vars} need to be 
#' aggregated, they are aggregated using \code{sum}.
#' See \code{?dcast}.
#' 
#' @examples 
#' \dontrun{
#' ## e.g. silly counts from a long-format table to a wide format
#' test <- copy(sire)
#' test$dg_y <- year(test$dg_date)
#' test$ex_y <- year(test$ex_date)
#' tab <- ltable(test, c("dg_y","ex_y"))
#' cast_simple(tab, columns='dg_y', rows="ex_y", values="obs")
#' }


cast_simple <- function(data=NULL, columns=NULL, rows=NULL, values=NULL) {
  if (!is.data.frame(data)) stop("data needs be a data.frame or data.table")
  if (is.null(data) || nrow(data) == 0L) stop("data NULL or has no rows")
  
  if (is.null(columns)) stop("columns cannot be NULL")
  
  msg <- paste0("Missing 'columns' variables: %%VARS%%")
  all_names_present(data, columns, msg = msg)
  msg <- paste0("Missing 'rows' variables: %%VARS%%")
  all_names_present(data, rows, msg = msg)
  msg <- paste0("Missing 'values' variables: %%VARS%%")
  all_names_present(data, values, msg = msg)
  
  ## allow rows = NULL 
  rowsNULL <- FALSE
  if (is.null(rows)) rowsNULL <- TRUE
  if (rowsNULL) rows <- "1"
  
  ## sometimes rows names appear to be like expressions, e.g. 'factor(V1)'
  ## (and this function only uses string-column-names, so that's fine.)
  actualRows <- rows
  if (length(rows) > 1L || rows != "1") {
    rows <- makeTempVarName(names = c(names(data), columns), 
                            pre = paste0("RN", 1:length(rows)))
    on.exit(setnames(data, rows, actualRows), add = TRUE)
    setnames(data, actualRows, rows)
  }
  ## same for cols
  actualCols <- columns
  columns <- makeTempVarName(names = c(names(data), rows), 
                             pre = paste0("CN", 1:length(columns)))
  on.exit(setnames(data, columns, actualCols), add = TRUE)
  setnames(data, actualCols, columns)
  
  form <- paste0(paste0(rows, collapse = " + "), " ~ ", 
                 paste0(columns, collapse = " + "))
  form <- as.formula(form)
  
  ## note: dcast probably usually finds the methods for data.frame / data.table,
  ## but this method is more certain
  if (is.data.table(data)) {
    d <- dcast.data.table(data, formula = form, value.var=values, 
                          drop=FALSE, fun.aggregate=sum)[]
  } else {
    d <- dcast(data, formula = form, value.var = values, 
               drop = FALSE, fun.aggregate = sum)[]
  }
  if (rowsNULL) set(d, j = names(d)[1L], value = NULL)
  wh_rows <- which(rows %in% names(d))
  if (sum(wh_rows, na.rm = TRUE)) setnames(d, rows[wh_rows], actualRows[wh_rows])
  
  d
}


#' @title Convert NA's to zero in data.table
#' @author Joonas Miettinen
#' @description Given a \code{data.table DT}, replaces any \code{NA} values
#' in the variables given in \code{vars} in \code{DT}. Takes a copy of the 
#' original data and returns the modified copy.
#' @import data.table
#' @param DT \code{data.table} object
#' @param vars a character string vector of variables names in \code{DT};
#' if \code{NULL}, uses all variable names in \code{DT}
#' @export na2zero
#' @details Given a \code{data.table} object, converts \code{NA} values
#' to numeric (double) zeroes for all variables named in \code{vars} or
#' all variables if \code{vars = NULL}.
na2zero = function(DT, vars = NULL) { 
  if (!is.data.table(DT)) stop("DT must be a data.table")
  DT <- copy(DT)
  
  navars <- vars
  if (is.null(navars)) navars <- names(DT)
  all_names_present(DT, navars)
  for (k in navars) {
    DT[is.na(get(k)), (k) := 0]
  }
  
  return(DT[])
}


#' @title Convert factor variable to numeric 
#' @description Convert factor variable with numbers as levels into a numeric variable
#' @param x a factor variable with numbers as levels
#' @export fac2num
#' @details
#' For example, a factor with levels \code{c("5","7")} is converted into 
#' a numeric variable with values \code{c(5,7)}.
#' @seealso
#' \code{\link{robust_values}}
#' @source 
#' \href{http://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-an-integer-numeric-without-a-loss-of-information}{Stackoverflow thread}
#' @examples
#' ## this is often not intended
#' as.numeric(factor(c(5,7))) ## result: c(1,2)
#' ## but this
#' fac2num(factor(c(5,7))) ## result: c(5,7)
#' 
#' ## however
#' as.numeric(factor(c("5","7","a"))) ## 1:3
#' 
#' fac2num(factor(c("5","7","a"))) ## result: c(5,7,NA) with warning
#' 
#' 
fac2num <- function(x) {
  as.numeric(levels(x))[x]
}



#' @title Convert date objects to fractional years
#' @author Joonas Miettinen
#' @description Using Date objects, calculates given 
#' dates as fractional years.
#' @param dates a vector or column of Date objects or right kind of character strings, see Details
#' @param format a character string; if \code{dates} is a character vector, 
#' specifies the format; see \code{\link{as.Date}}
#' @param year.length character string, either \code{'actual'} or 
#' \code{'approx'}; can be abbreviated; see Details
#' @import data.table
#' @export get.yrs
#' @details
#' 
#' \code{dates} should preferably be a \code{date}, \code{Date} or \code{IDate} 
#' object, 
#' although they can also be character strings in a format
#' specified by \code{format} (passed to \code{\link{as.Date}}).
#' 
#' When \code{ year.length = 'actual' }, fractional years are calculated as 
#' \code{ year + day_in_year/365 } for non-leap-years
#' and as \code{ year + day_in_year/366 } for leap years. 
#' If \code{ year.length = 'approx' }, fractional years are always
#' calculated as in \code{ year + day_in_year/365.242199 }. 
#' 
#' There is a slight difference, then, between the two methods
#' when calculating durations between fractional years. For
#' meticulous accuracy one might instead want to calculate durations using
#' dates (days) and convert the results to fractional years.
#'  
#' Note that dates are effectively converted to fractional years at 
#' \code{ 00:00:01 } o'clock:
#' 
#' 
#' \code{ get.yrs("2000-01-01") = 2000 }, and
#' \code{ get.yrs("2000-01-02") = 2000 + 1/365.242199 }. 
#' 
#' 
#' @seealso
#' \code{\link[Epi]{cal.yr}}, \code{\link{as.Date.yrs}}
#' 
#' @examples
#' 
#' test <- copy(sire)
#' test$dg_yrs <- get.yrs(test$dg_date)
#' summary(test$dg_yrs)
#' 
#' ## see: ?as.Date.yrs
#' dg_date2 <- as.Date(test$dg_yrs)
#' summary(as.numeric(dg_date2 - test$dg_date))
#' 
#' ## Epi's cal.yr versus get.yrs
#' Epi::cal.yr("2000-01-01") ## 1999.999
#' get.yrs("2000-01-01") ## 2000
#' 
get.yrs <- function(dates, format = "%Y-%m-%d", year.length = "approx") {
  year.length <- match.arg(year.length, c("actual", "approx"))
  y <- yrs <- NULL ## to instate as global variable to appease R CMD CHECK
  
  if (missing(dates) || length(dates) == 0L) 
    stop("'dates' not supplied or is a variable with length zero.")
  
  orle <- length(dates)
  nale <- sum(is.na(dates))
  
  dat <- data.table(dates=dates)
  if (is.character(dates)) {
    dat[, dates := as.IDate(dates, format = format)]
  } else {
    dat[, dates := as.IDate(dates)]
  }
  
  
  nale2 <- dat[, sum(is.na(dates))]
  
  
  if (year.length  == "actual") {
    ## fractional years using actual year lengths
    dat[, y    := year(dates)]
    ## calculate distance between first day of each year and given dates values
    dat[            , yrs := as.numeric(y + (yday(dates)-1L)/365L)]
    dat[is_leap_year(y), yrs := as.numeric(y + (yday(dates)-1L)/366L)]
    
  } else {
    ## fractional years using hypothetical year length    
    dat[, yrs := year(dates) + (yday(dates)-1L)/365.242199]
  }
  nale3 <- dat[, sum(is.na(yrs))]
  if (nale3 >  nale) {
    warning(nale3-nale, " values were coerced to NA by get.yrs, ",
            "of which ", nale2-nale, " were coerced to NA by as.Date(). ",
            "Please use any date class supported by as.Date(), or ",
            "supply the appropriate 'format' if 'dates' is a character ", 
            "vector. See ?as.Date for more information on formatting.")
  }
  
  setattr(dat$yrs, "year.length", year.length)
  setattr(dat$yrs, "class", c("yrs", "numeric"))
  
  return(dat$yrs)
}



#' @title Detect leap years
#' @author Joonas Miettinen
#' @description Given a vector or column of year values (numeric or integer), \code{\link{is_leap_year}} returns a vector of equal length
#' of logical indicators, i.e. a vector where corresponding leap years have value TRUE, and FALSE otherwise.
#' 
#' @param years a vector or column of year values (numeric or integer)
#' @examples
#' ## can be used to assign new columns easily, e.g. a dummy indicator column
#' df <- data.frame(yrs=c(1900,1904,2005,1995))
#' df$lyd <- as.integer(is_leap_year(df$yrs))
#' 
#' ## mostly it is useful as a condition or to indicate which rows have leap years
#' which(is_leap_year(df$yrs)) # 2
#' df[is_leap_year(df$yrs),] # 2nd row
#' 
#' @export is_leap_year
#' 
is_leap_year <- function(years) {
  if (!is.numeric(years)) {
    stop("years must be a numeric vector, preferably integer for speed. Use e.g. as.integer().")
  }
  
  years <- try2int(years)
  if (!is.integer(years)) stop("years could not be coerced to integer; don't use fractional years such as 2000.1234 but integers such as 2000")
  
  # divisible by four
  isLeap <- years %% 4L == 0L
  # not divisible by 100
  isLeap <- isLeap & years %% 100L != 0L
  # unless divisible by 400 also
  isLeap <- isLeap | years %% 400L == 0L
  isLeap
  
}
#' @title Test if object is a \code{Date} object
#' @description Tests if an object is a \code{Date} object and returns
#' a logical vector of length 1. \code{IDate} objects are also 
#' \code{Date} objects, but \code{date} objects from package \pkg{date}
#' are not. 
#' @author Joonas Miettinen
#' @param obj object to test on
#' @export is.Date
#' @seealso
#' \code{\link{get.yrs}}, \code{\link{is_leap_year}}, \code{\link{as.Date}}
#' @examples
#' ## the base "capital Date" format
#' da <- as.Date("2000-01-01")
#' is.Date(da) ## TRUE
#' date::is.date(da) ## FALSE
#' 
#' ## IDate format from data.table
#' da <- as.IDate("2000-01-01")
#' is.Date(da) ## TRUE
#' date::is.date(da) ## FALSE
#' 
#' ## from package "date"
#' da <- date::as.date("1jan2000")
#' is.Date(da) ## FALSE
#' date::is.date(da) ## TRUE
#'  
is.Date <- function(obj) {
  
  if (any(c("IDate","Date") %in% class(obj))) {
    return(TRUE)
  }
  
  return(FALSE)
}


#' @title Convert values to numeric robustly
#' @author Joonas Miettinen
#' 
#' @param num.values values to convert to numeric
#' @param force logical; if \code{TRUE}, returns a vector of values where values that cannot be interpreted as numeric are
#' set to \code{NA}; if \code{FALSE}, returns the original vector and gives a warning if any value cannot be interpreted as
#' numeric.
#' @param messages logical; if \code{TRUE}, returns a message of what was done with the \code{num.values}
#' @description Brute force solution for ensuring a variable is numeric by 
#' coercing a variable of any type first to factor and then to numeric
#' @export robust_values
#' @import data.table
#' @note
#' Returns \code{NULL} if given \code{num.values} is \code{NULL}. 
#' @examples
#' ## this works
#' values <- c("1", "3", "5")
#' values <- robust_values(values)
#' 
#' ## this works
#' values <- c("1", "3", "5", NA)
#' values <- robust_values(values)
#' 
#' ## this returns originals
#' values <- c("1", "3", "5", "a")
#' values <- robust_values(values)
#' 
#' ## this forces "a" to NA and works otherwise
#' values <- c("1", "3", "5", "a")
#' values <- robust_values(values, force=TRUE)
#' 

robust_values <- function(num.values, force = FALSE, messages = TRUE) {
  a <- NULL
  if (is.null(num.values)) {
    return(NULL)
  }
  dt <- data.table(num.values)
  nas <- dt[is.na(num.values), .N]
  
  suppressWarnings(
    dt[,a := fac2num(factor(num.values))]
  )
  dt[, a := try2int(a)]
  nas2 <- dt[is.na(a), .N]
  
  if (!force & nas2 > nas) {
    if (messages) warning("since force = FALSE and NAs were created, returning original values")
    return(dt$num.values)
  }
  if (force) {
    if (nas2 > nas) {
      if (messages) warning("some NAs were created")
    }
    return(dt$a)
  }
  
  
  return(dt$a)
  
  
}

#' @title Check if all names are present in given data
#' @author Joonas Miettinen
#' @param data dataset where the variable names should be found
#' @param var.names a character vector of variable names, e.g.
#' \code{c("var1", "var2")}
#' @param stops logical, stop returns exception
#' @param msg Custom message to return instead of default message.
#' Special: include \code{\%\%VARS\%\%} in message string and the missing 
#' variable names will be inserted there (quoted, separated by comma, e.g. 
#' \code{'var1'}, \code{'var2'} --- no leading or tracing white space). 
#' @description Given a character vector, checks if all names are present in \code{names(data)}.
#' Throws error if \code{stops=TRUE}, else returns \code{FALSE} if some variable name is not present.
#' @seealso
#' \code{\link{robust_values}}
#' @export all_names_present

all_names_present <- function(data, var.names, stops = TRUE, msg = NULL) {
  
  if (!is.null(var.names) && !is.character(var.names)) {
    stop("Argument 'var.names' must be NULL or a character vector of ",
         "variable names.")
  }
  if (length(var.names) && any(is.na(var.names))) {
    stop("There are ", sum(is.na(var.names)), " missing values in argument ",
         "'var.names'. Please only supply non-NA values.")
  }
  
  badNames <- setdiff(var.names, names(data))
  if (length(badNames) == 0L) return(TRUE)
  
  badNames <- paste0("'", badNames, "'", collapse = ", ")
  
  if (is.null(msg)) msg <- paste0("Cannot proceed - following given variable name(s) not present in dataset '",
                                  deparse(substitute(data)), "': ", badNames)
  if (!is.character(msg) || length(msg) > 1L) stop("Argument 'msg' must be a character string vector of length one.") else
    msg <- gsub(pattern = "%%VARS%%", replacement = badNames, x = msg)
  if (!is.logical(stops) || length(stops) > 1L) stop("Argument 'stops' must be either TRUE or FALSE.")
  
  if (stops) stop(msg)
  
  return(FALSE)
}


#' @title Return lower_bound value from char string (20,30]
#' @author Matti Rantanen
#' @description selects lowest values of each factor after cut() based
#' on that the value starts from index 2 and end in comma ",".
#' @param cut is a character vector of elements "(20,60]"
#' @export lower_bound

lower_bound <- function(cut) {
  cut <- as.character(cut)
  ind <- gregexpr(pattern=',',cut)
  ind <- as.numeric(ind) - 1
  t.sub <- as.numeric(substr(cut,2, ind))
  return(t.sub)
}


#' @title Change output values from cut(..., labels = NULL) output
#' @author Matti Rantanen
#' @param t is a character vector of elements, e.g. "(20,60]"
#' @param factor logical; TRUE returns informative character string, FALSE numeric (left value)
#' @description Selects lowest values of each factor after cut() based
#' on the assumption that the value starts from index 2 and end in comma ",".
#' @details type = 'factor': "[50,52)" -> "50-51" OR "[50,51)" -> "50"
#' 
#' type = 'numeric': lowest bound in numeric.
#' 
#' @export cut_bound
#' @examples
#' cut_bound("[1900, 1910)") ## "1900-1909"

cut_bound <- function(t, factor=TRUE) {
  if (!factor) {
    t <- as.character(t)
    ind <- gregexpr(pattern=',',t)
    ind <- as.numeric(ind) - 1
    t <- as.numeric(substr(t,2, ind))
    return(t)
  }
  if (factor) {
    t <- as.character(t)
    t <- gsub(',', '-' , substr(t, 2, nchar(t) - 1) )
    ind <-as.numeric( gregexpr(pattern='-',t) )
    if (any(as.numeric( substr(t,1,ind-1) ) +1 == as.numeric( substr(t,ind+1,nchar(t))) ) ) {
      t <- substr(t,1,ind-1)
      return(t)
    }
    t
    a <- substr(t, ind+1, nchar(t))
    t <- sub(a, as.character(as.numeric(a)-1), t)
    return(t)
  }
}




#' @title Set the class of an object (convencience function for
#'  \code{setattr(obj, "class", CLASS)}); can add instead of replace
#' @description Sets the class of an object in place to \code{cl}
#' by replacing or adding
#' @param obj and object for which to set class
#' @param cl class to set
#' @param add if \code{TRUE}, adds \code{cl} to the 
#' classes of the \code{obj}; otherwise replaces the class information
#' @param add.place \code{"first"} or \code{"last"}; adds \code{cl}
#' to the front or to the back of the \code{obj}'s class vector
#' @author Joonas Miettinen
setclass <- function(obj, cl, add=FALSE, add.place="first") {
  match.arg(add.place, c("first","last"))
  cl <- as.character(cl)
  
  if (add) {
    old_classes <- attr(obj, "class")
    
    if (add.place=="first") {
      setattr(obj, "class", c(cl, old_classes))
    } else {
      setattr(obj, "class", c(old_classes, cl))
    }
  } else {
    setattr(obj, "class", cl)
  }
}




#' @title Attempt coercion to integer
#' @author James Arnold
#' @description Attempts to convert a numeric object to integer, 
#' but won't if loss of information is imminent (if values after decimal
#' are not zero for even one value in \code{obj})
#' @param obj a numeric vector
#' @param tol tolerance; if each numeric value in \code{obj} deviate from
#' the corresponding integers at most the value of \code{tol}, they are considered
#' to be integers; e.g. by default \code{1 + .Machine$double.eps} is considered
#' to be an integer but \code{1 + .Machine$double.eps^0.49} is not.
#' @export try2int
#' @source \href{http://stackoverflow.com/questions/3476782/how-to-check-if-the-number-is-integer}{Stackoverflow thread}
try2int <- function(obj, tol = .Machine$double.eps^0.5) {
  if (!is.numeric(obj)) stop("obj needs to be integer or double (numeric)")
  if (is.integer(obj)) return(obj)
  
  test <- FALSE
  
  bad <- if (length(na.omit(obj)) == 0) TRUE else 
    min(obj, na.rm = TRUE) == -Inf || max(obj, na.rm = TRUE) == Inf
  if (bad) {
    return(obj)
  } else {
    test <- max(abs(obj) %% 1, na.rm = TRUE) < tol
  }
  
  if (is.na(test) || is.null(test)) test <- FALSE
  
  if (test) return(as.integer(obj))
  
  return(obj)
  
}


#' @title Shift a variable to create lag or lead values
#' @author Joonas Miettinen
#' @description 
#' \strong{DEPRECATED}: 
#' Intended to do what \code{\link[data.table]{shift}} from
#' \pkg{data.table} does better since \pkg{data.table} 1.9.6. 
#' Shifts the values of a variable forwards or 
#' backwards to create lag or lead values. Takes a copy of the whole data
#' and returns a new copy with the shifted variable.
#' @export shift.var
#' @param data a \code{data.frame} or \code{data.table}
#' @param id.vars a character string vector of variable names; \code{id.vars} are used to identify unique subjects,
#' for which shifting is done separately; e.g. with a panel data where \code{region} refers to different regions that
#' all have their own time series, using \code{id.vars = "region"} shifts the time series for each region separately
#' @param shift.var a character string vector of length one; specifies the variable according to which \code{value.vars}
#' are shifted; e.g. \code{id.vars = "year"} means shifting forward or backward in years (given one has a var name \code{"year"})
#' @param value.vars a character string vector; specifies the names of variables whose values that are shifted
#' @param shift.value an integer; specifies the direction and extent of shifting; e.g. \code{shift.value = -1L} shifts
#' one row backwards (a lag of one row) and \code{shift.value = 2L} creates a two-row lead
shift.var <- function(data, id.vars = NULL, shift.var = NULL, value.vars=NULL, shift.value=-1L) {
  .Deprecated(new = "shift", msg = "popEpi's shift.var is deprecated in 0.3.0 and will be removed in the next release; please use e.g. data.table's shift() function")
  
  merge_var <- makeTempVarName(data, pre = "merge_var")
  if (is.null(shift.var)||is.null(value.vars)) stop("shift.var and value.vars cannot be NULL")
  all_names_present(data, c(id.vars, shift.var, value.vars))
  if (shift.value == 0L) return(data)
  
  if (is.data.table(data)) old_key <- key(data)
  
  data <- data.table(data)
  setkeyv(data, c(id.vars, shift.var))
  data[, (merge_var) := as.integer(as.factor(get(shift.var)))]
  
  if (any(duplicated(data, by=c(id.vars,shift.var,value.vars)))) {
    stop("some levels of shift.var are duplicated in data, so shifting is not possible")
  }
  
  lagdata <- data[,c(id.vars, merge_var, value.vars), with=FALSE]
  lagdata[, (merge_var) := get(merge_var) - shift.value]
  
  if (shift.value<=0) {vn <- "lag"}
  if (shift.value >0) {vn <- "lead"}
  
  setnames(lagdata, value.vars, paste0(vn, abs(shift.value),"_", value.vars))
  
  setkeyv(data, c(id.vars, merge_var))
  setkeyv(lagdata, c(id.vars, merge_var))
  data <- lagdata[data]
  #   data <- merge(data, lagdata, all.x=TRUE, all.y=FALSE, by = c(id.vars, merge_var))
  
  data[, (merge_var) := NULL]
  setkeyv(data, old_key)
  return(data[])
}



#' @title Get rate and exact Poisson confidence intervals
#' @author epitools
#' @description Computes confidence intervals for Poisson rates
#' @param x observed
#' @param pt expected
#' @param conf.level alpha level
#' 
#' @export poisson.ci
#' 
#' 
#' @examples
#' 
#' poisson.ci(x = 4, pt = 5, conf.level = 0.95)
#'
poisson.ci <- function(x, pt = 1, conf.level = 0.95) {
  xc <- cbind(x, conf.level, pt)
  pt2 <- xc[, 3]
  results <- matrix(NA, nrow(xc), 6)
  f1 <- function(x, ans, alpha = alp) {
    ppois(x, ans) - alpha/2
  }
  f2 <- function(x, ans, alpha = alp) 1 - ppois(x, ans) + dpois(x, ans) - alpha/2
  for (i in 1:nrow(xc)) {
    alp <- 1 - xc[i, 2]
    interval <- c(0, xc[i, 1] * 5 + 4)
    uci <- uniroot(f1, interval = interval, x = xc[i, 1])$root/pt2[i]
    if (xc[i, 1] == 0) {
      lci <- 0
    }
    else {
      lci <- uniroot(f2, interval = interval, x = xc[i,1])$root/pt2[i]
    }
    results[i, ] <- c(xc[i, 1], pt2[i], xc[i, 1]/pt2[i], lci, uci, xc[i, 2])
  }
  coln <- c("x", "pt", "rate", "lower", "upper", "conf.level")
  colnames(results) <- coln
  data.frame(results)
}


#' @title Delete \code{data.table} columns if there
#' @author Joonas Miettinen
#' @description Deletes columns in a \code{data.table} conveniently.
#' May only delete columns that are found silently. Sometimes useful in e.g.
#' \code{on.exit} expressions.
#' @param DT a \code{data.table}
#' @param delete a character vector of column names to be deleted
#' @param keep a character vector of column names to keep; 
#' the rest will be removed; \code{keep} overrides \code{delete}
#' @param colorder logical; if \code{TRUE}, also does \code{setcolorder} using
#' \code{keep}
#' @param soft logical; if \code{TRUE}, does not cause an error if any variable
#' name in \code{keep} or \code{delete} is missing; \code{soft = FALSE} useful 
#' for programming sometimes
#' 
#' 
#' @export setcolsnull
setcolsnull <- function(DT=NULL, delete=NULL, keep=NULL, colorder=FALSE, soft=TRUE) {
  if (!is.data.table(DT)) stop("not a data.table")
  if (!soft) {
    all_names_present(DT, keep, msg = "Expected")
    all_names_present(DT, delete)
  }
  del_cols <- NULL
  del_cols <- intersect(delete, names(DT))
  if (!is.null(keep)) {
    del_cols <- setdiff(names(DT), keep)
  }
  if (length(del_cols) > 0) {
    set(DT, j = (del_cols), value = NULL)
  }
  if (colorder) {
    setcolorder(DT, intersect(keep, names(DT)))
  }
  return(invisible())
}






#' @title Coerce a \code{ratetable} Object to Class \code{data.frame}
#' @description
#' \code{ratatable} objects used in e.g. \pkg{survival} and \pkg{relsurv}
#' can be conveniently coerced to a long-format \code{data.frame}.
#' However, the names and levels of variables in the result
#' may not match names and levels of variables in your data.
#' @author Joonas Miettinen
#' @param x a \code{ratetable}
#' @param ... unused but added for compatibility with \code{as.data.frame}
#' @examples
#' library(relsurv)
#' data(slopop)
#' df <- as.data.frame(slopop)
#' head(df)
#' @seealso 
#' \code{\link[survival]{ratetable}}, 
#' \code{\link{as.data.table.ratetable}}
#'
#' @export
as.data.frame.ratetable <- function(x, ...) {
  dimids <- attr(x, "dimid")
  x <- as.data.frame.table(as.table(as.array(x)))
  names(x) <- c(dimids, "haz")
  x[]
}


#' @title Coerce a \code{ratetable} Object to Class \code{data.table}
#' @author Joonas Miettinen
#' 
#' @description
#' \code{ratatable} objects used in e.g. \pkg{survival} and \pkg{relsurv}
#' can be conveniently coerced to a long-format \code{data.frame}.
#' However, the names and levels of variables in the result
#' may not match names and levels of variables in your data.
#' @param x a \code{ratetable}
#' @param ... other arguments passed on to \code{as.data.table}

#' @seealso 
#' \code{\link[survival]{ratetable}}, 
#' \code{\link{as.data.frame.ratetable}}
#'
#' @examples
#' library(relsurv)
#' data(slopop)
#' dt <- as.data.table(slopop)
#' dt
#' @export
as.data.table.ratetable <- function(x, ...) {
  dimids <- attr(x, "dimid")
  x <- as.data.table(as.table(as.array(x)), ...)
  x[, names(x) := lapply(.SD, robust_values, messages = FALSE, force = FALSE)]
  setnames(x, c(dimids, "haz"))
  x[]
}


#' @title \strong{Experimental}: Coerce a long-format \code{data.frame} to a \code{ratetable} object
#' @author Joonas Miettinen
#' @description Coerces a long-format \code{data.frame} of population hazards
#' to an array, and in turn to a \code{\link[survival]{ratetable}},
#' which can be used in e.g. \pkg{survival}'s expected survival computations
#' and \pkg{relsurv}'s relative survival computations.
#' @param DF a \code{data.frame}
#' @param value.var name of values variable in quotes
#' @param by.vars names vector of variables by which to create (array) dimensions
#' @seealso 
#' \code{\link[survival]{ratetable}}, 
#' \code{\link{as.data.table.ratetable}}, 
#' \code{\link{as.data.frame.ratetable}}
#'
longDF2ratetable <- function(DF, value.var = "haz", by.vars = setdiff(names(DF), value.var)) {
  univals <- lapply(DF[, by.vars], unique)
  names(univals) <- NULL
  dimvec <- sapply(DF[,by.vars], function(x) {length(unique(x))},
                   simplify=TRUE)
  ar <- array(DF[, value.var], dim = dimvec)
  dimnames(ar) <- univals
  attr(ar, "class") <- "ratetable"
  attr(ar, "dimid") <- colnames(DF)
  ar
}

temp_var_names <- function(n = 1L, avoid = NULL, length = 10L) {
  ## INTENTION: make temporary variable names that don't exist in
  ## char vector "avoid", e.g. avoid = names(data).
  if (n < 1L || !is.integer(n)) {
    stop("n must an integer > 0")
  }
  if (length < 1L || !is.integer(length)) {
    stop("length must an integer > 0")
  }
  if (!is.null(avoid)) avoid <- as.character(avoid)
  
  pool <- c(0:9, letters, LETTERS)
  
  formTemp <- function(int) {
    v <- sample(x = pool, size = length, replace = TRUE)
    paste0(v, collapse = "")
  }
  
  l <- lapply(1:n, formTemp)
  dupll <- duplicated(l) | l %in% avoid
  tick <- 1L
  while (any(dupll) && tick <= 100L) {
    l[dupll] <- lapply(1:sum(dupll), formTemp)
    dupll <- duplicated(l) | l %in% avoid
    tick <- tick + 1L
  }
  if (tick >= 100L) {
    stop("ran randomization 100 times and could not create unique temporary",
         " names. Perhaps increase length?")
  }
  unlist(l)
}

#' @import stats
makeTempVarName <- function(data=NULL, names=NULL, 
                            pre=NULL, post=NULL, length = 10L) {
  DN <- NULL
  DN <- c(DN, names(data))
  DN <- c(DN, names)
  DN <- unique(DN)
  
  if (length(pre) != length(post) && length(post) > 0L && length(pre) > 0L) {
    stop("Lengths of arguments 'pre' and 'post' differ (", length(pre), " vs. ",
         length(post), "). (Tried to create temporary variables, so this is ",
         "most likely an internal error and the pkg maintainer should be ",
         "complained to.)")
  }
  useN <- max(length(pre), length(post), 1L)
  useL <- length
  tv <- temp_var_names(avoid = DN, n = useN, length = useL)
  tv <- paste0(pre, tv, post)
  tv
}


setDFpe <- function(x) {
  ## intended to only be used to set data.table to data.frame in place
  ## when option("popEpi.datatable") == FALSE
  if (!is.data.table(x)) stop("only accepts data.table as input")
  
  cl <- class(x)
  wh <- which(cl == "data.table")
  cl = c(cl[1:(wh-1)], cl[(wh+1):length(cl)])
  setattr(x, "class", cl)
  
  setattr(x, "sorted", NULL)
  setattr(x, ".internal.selfref", NULL)
}



evalLogicalSubset <- function(data, substiset, n = 2, enclos = parent.frame(n)) {
  ## NOTE: subset MUST be substitute()'d before using this function!
  ## we allow substiset to be a logical condition only
  ## ALWAYS returns a logical vector of length nrow(data)
  
  substiset <- eval(substiset, envir = data, enclos = enclos)
  if (!is.null(substiset)) {
    if (!is.logical(substiset)) stop("Expression to subset by must be a logical condition, e.g. var1 == 0, var1 %in% 1:2, var1 > 0, etc.")
    substiset <- substiset & !is.na(substiset)
    if (sum(substiset) == 0) stop("zero rows in data after subset")
  } else {
    substiset <- rep(TRUE, nrow(data))
  }
  substiset
}


subsetDTorDF <- function(data, subset=NULL, select=NULL) {
  ## INTENTION: subsetting either a data.table or a data.frame
  ## and returning only selected variables for lazy people.
  if (!is.data.frame(data)) stop("data must be a data.table/data.frame")
  if (!is.logical(subset) && !is.null(subset)) stop("subset must be a logical vector or NULL")
  
  if (is.null(select)) {
    select <- names(data)
  } else {
    all_names_present(data, select)
  }
  
  e <- "data["
  if (!is.null(subset) && !all(subset)) e <- paste0(e, "subset") 
  if (!is.null(select) && (length(select) < names(data) || any(select != names(data)))) {
    e <- paste0(e, ", eval(select)")
    if (is.data.table(data)) e <- paste0(e, ", with = FALSE")
  }
  e <- paste0(e, "]")
  
  e <- parse(text = e)
  
  eval(e)
  
}

subsetRolling <- function(data, subset = NULL, select = NULL) {
  ## INTENTION: subsets a data.table column by column and by deleting columns
  ## in the old data.table.
  if (!is.data.table(data)) stop("data must be a data.table")
  if (!is.logical(subset)) stop("subset must be a logical vector")
  
  if (is.null(select)) {
    select <- names(data)
  } else {
    all_names_present(data, select)
  }
  
  if (length(select) == 0L) stop("select is of length zero, which would delete all columns in data")
  
  setcolsnull(data, keep = select)
  
  dt <- data[subset, select[1L], with = FALSE]
  
  setcolsnull(data, delete = select[1L])
  select <- select[-1L]
  
  for (v in select) {
    set(dt, j = v, value = data[[v]][subset])
    set(data, j = v, value = NULL)
  }
  
  rm(list = deparse(substitute(data)), envir = parent.frame(1L))
  
  dt
}



setDT2DF <- function(x) {
  if (!is.data.table(x)) stop("only accepts data.table as input")
  
  cl <- class(x)
  cl <- setdiff(cl, "data.table")
  setattr(x, "class", cl)  
  setattr(x, "sorted", NULL)
  setattr(x, ".internal.selfref", NULL)
  invisible(x)
}

setDF2DT <- function(x) {
  if (!is.data.frame(x) || is.data.table(x)) stop("only accepts data.frame as input")
  
  cl <- class(x)
  whDF <- which(cl == "data.frame")
  cl <- c(cl[1:(whDF-1)], "data.table", "data.frame", cl[whDF:length(cl)])
  
  setattr(x, "class", cl)
  alloc.col(x)
  
  invisible(x)
}




p.round <- function(p, dec=3) {
  th <- eval( parse(text=paste0('1E-', dec ) ))
  if( is.null(p)) return( '= NA') 
  if( is.na(p))   return( '= NA') 
  if( p < th ){
    p <- paste0('< ', th  )
  } else {
    p <- paste0('= ', round(p, dec) )
  }
  p 
}

popArg2ModelNames <- function(arg, type) {
  ## INTENTION: given a quoted/substituted expression,
  ## digs out the expression(s) creating a/multiple column(s)
  ## and returns the deparsed expression(s) to be used as names
  ## of columns the same way that models such as lm() display
  ## the names of expressions used within formula
  
  ## some exceptions
  if (is.data.frame(arg)) return(names(arg))
  if (is.character(arg)) return(arg)
  
  type <- match.arg(type[1L], c("NULL", "character", "list", "expression", "formula"))
  
  lang <- NULL
  lang <- try(is.language(arg) || inherits(arg, "formula"), silent = TRUE)
  
  
  if (inherits(lang, "try-error") || !lang) stop("arg must be a quoted or substituted expression or a formula. Error message: ", lang, ". type of arg: ", typeof(arg), ". Class: ", class(arg), ". Mode: ", mode(arg), ".")
  
  d <- oneWhitespace(paste0(deparse(arg)))
  
  if (type == "expression") return(d) else 
    if (type == "NULL") return(NULL) else 
      if (type == "character") return(eval(arg)) else 
        if (type == "list") {
          d <- substr(d, 6, nchar(d)-1L) ## removes "list(" and ")"
          d <- strsplit(d, ", ")
          return(unlist(d))
        } else if (type == "formula") {
          arg <- eval(arg)
          d <- names(RHS2list(arg))
          if (length(d) == 0L) return(NULL) ## e.g. y ~ 1
          return(d)
        }
  stop("could not determine deparsed-expression-names")
}

uses_dollar <- function(q, data.names) {
  ## INTENTION: determine whether q is an expressions that is evaluated
  ## outside a data.frame, i.e. one that uses the dollar operator.
  ## e.g. TF$V1 should not be evaluated in a data.frame even if it has
  ## the variables TF and V1 since it wont work and was not intended.
  if (!is.language(q) || inherits(q, "formula")) {
    return(FALSE)
  }
  
  d <- deparse(q)
  ## sometimes d is of length > 1 for some reason...
  d <- paste0(d, collapse = "")
  d <- oneWhitespace(d)
  
  if (substr(d, 1, 4) == "list") {
    ## lists are not allowed to work in this manner for now.
    return(FALSE)
  }
  
  if (!grepl(x = d, pattern = "\\$")) {
    ## does not use dollar operator.
    return(FALSE)
  }
  
  ## detect if word$word is used in d
  t <- regexec(pattern = "\\w+\\$\\w+", text = d)
  if (t != -1) {
    ## ok, used word$word
    ## is there are variable with that name in data.names?
    m <- unlist(regmatches(d, t))
    if (m %in% data.names) {
      return(FALSE)
    }
    ## if not, it should be evaluated outside the data.
    return(TRUE)
  } 
  
  return(FALSE)
}


evalPopArg <- function(data, arg, n = 1L, DT = TRUE, enclos = NULL, recursive = TRUE, types = c("NULL","character", "list", "expression"), naming = c("DT", "model")) {
  ## arg: an unevaluated AND substitute()'d argument within a function, which may be
  ## * an expression
  ## * a list of expressions
  ## * a character vector of variable names (in a given data set)
  ## n: steps upstream as in parent.frame(n); 0L refers to calling environment
  ## of evalPopArg, 1L to calling environment of e.g. sir which uses evalPopArg, etc.
  ## hence n = 1L should be almost always the right way to go.
  ## ALTERNATIVELY supply an environment by hand via enclos.
  ## enclos will override n.
  ## recursive: if TRUE, evals arg as many times as it is of type language.
  ## output:
  ## * vector as a result of an expression
  ## * list as a result of a list
  ## * character vector of names
  ## OR with DT = TRUE, a data.table based on aforementioned results.
  ## intention: output to be used in by argument of data.table.
  ## a data.table output is directly usable in by.
  ## if column names cannot be easily found, BV1, BV2, ... are imputed
  ## for missing names (unrobustly: such names may already exist, resulting in duplicates)
  
  ## naming: DT style uses first element of all.names() where 
  ## a name has to be created; model style keeps the whole deparsed
  ## expression. Only applied when DT = TRUE
  naming <- match.arg(naming[1L], c("DT", "model"))
  
  ## types: allowed popArg types of arguments.
  types <- match.arg(types, c("NULL","character", "list", "expression", "formula"), several.ok = TRUE)
  
  if (!is.null(enclos) && !is.environment(enclos)) {
    stop("enclos must be NULL or an environment")
  }
  if (!is.environment(enclos)) enclos <- parent.frame(n + 1L)
  
  ## used data may change if expression uses dollar operator, hence
  ## arg should not be evaluated within data but only its surroundings.
  use_data <- data
  use_enc <- enclos
  dataNames <- names(data)
  
  if (uses_dollar(arg, data.names = dataNames)) {
    use_data <- enclos
    use_enc <- baseenv()
  }
  e <- eval(arg, envir = use_data, enclos = use_enc)
  if (is.language(e) && !inherits(e, "formula")) {
    if (!recursive) stop("arg is of type language after evaluating, and recursive = FALSE")
    
    tick <- 1L
    while (is.language(e) && !inherits(e, "formula") && tick < 100L) {
      arg <- e
      use_data <- data
      use_enc <- enclos
      if (uses_dollar(arg, data.names = dataNames)) {
        use_data <- enclos
        use_enc <- baseenv()
      }
      e <- eval(arg, envir = use_data, enclos = use_enc)
      tick <- tick + 1L
    }
    if (tick == 100L) stop("arg was of type language even after 100 evaluations. Something went wrong here...")
    
    
    
  } 
  argType <- "NULL"
  if (is.list(e)) argType <- "list" else 
    if (is.character(e)) argType <- "character" else 
      if (is.vector(e) || is.factor(e)) argType <- "expression" else 
        if (inherits(e, "formula")) argType <- "formula"
  
  if (!argType %in% types) stop("popArg type of evaluated arg not one of the allowed types (set via argument types). Detected type: '", argType, "'. Allowed types: ", paste0("'", types, "'", collapse = ", "))
  
  if (argType == "NULL") return(NULL)
  
  av <- all.vars(arg)
  if (argType == "character") av <- e
  
  ## byNames: names of columns resulting from aggre argument, by which
  ## pyrs and such are aggregated. same functionality
  ## as in results seen in e.g.DT[, .N, by = list(factor(x), y, z = w)] ## factor, y, z
  ## note: first object in ags with list or expression aggre is "list"
  byNames <- NULL
  
  if (is.character(e)) byNames <- e
  else if (argType == "list" && substr(paste0(deparse(arg)), 1, 5) == "list(") byNames <- sapply(arg[-1], function(x) all.names(x)[1]) 
  else if (argType == "expression") byNames <- all.names(arg)[1]
  
  badNames <- c("$", ":")
  
  byNames[byNames %in% badNames] <- paste0("BV", 1:length(byNames))[byNames %in% badNames]
  
  
  if (argType == "formula") {
    arg <- e
    use_data <- data
    use_enc <- enclos
    e <- RHS2DT(formula = e, data = use_data, enclos = use_enc)
    if (ncol(e) == 0L || nrow(e) == 0L) e <- data.table() ## e.g. y ~ 1
    
  } else if (is.character(e)) {
    all_names_present(data, e)
    if (DT) {
      ## note: e contains variable names in character strings,
      ## ergo fully named list & DT created
      l <- lapply(e, function(x) data[[x]])
      setattr(l, "names", e)
      setDT(l)
      e <- l; rm(l)
    }
  } else if (is.list(e)) {
    ## note: fully unnamed list has NULL names()
    ## partially named list has some "" names
    ne <- names(e)
    
    if (DT && any(sapply(e, is.null))) stop("at least one object in list arg is NULL; cannot form data.table with such list")
    
    if (is.null(ne)) ne <- rep("", length(e))
    
    
    wh_bad <- which(ne == "")
    if (length(wh_bad) > 0) {
      if (is.null(byNames)) {
        byNames <- paste0("BV", 1:length(e))
      }
      
      ne[wh_bad] <- byNames[wh_bad]
      setattr(e, "names", ne)
    }
    
    if (DT) {
      ## NOTE: used to be setDT, but length of different elements
      ## in list may differ, which as.data.table handles correctly
      e <- as.data.table(e)
    }
  } else if ((is.vector(e) || is.factor(e))) {
    ## is e.g. a numeric vector or a factor
    if (DT) {
      e <- data.table(V1 = e)
      setnames(e, 1, byNames)
    }
  }
  
  ## NOTE: e may be of type language at this point if arg was double-quoted
  ## and recursive = FALSE
  
  if (DT) {
    setDT(e)
    setattr(e, "all.vars", av)
    setattr(e, "quoted.arg", arg)
    setattr(e, "arg.type", argType)
    if (naming == "model" && ncol(e) > 0L) setnames(e, 1:ncol(e), popArg2ModelNames(arg, type = argType))
  }
  e
}


popArgType <- function(arg, data = NULL, n = 1L, enclos = NULL, recursive = TRUE) {
  ## input: a substitute()'d expression / argument
  ## NOTE: recursive useful when arg might be quoted twice and want the eventual
  ## result; need to supply data for it though
  ## output: type of thingie that was substitute()'d
  ##  * list (of expressions)
  ##  * character string vector
  ##  * an expression (includes symbol)
  av <- all.vars(arg, unique = TRUE) ## all variables
  av <- setdiff(av, c("$", "T", "F"))
  an <- all.names(arg, unique = TRUE) ## all variables and functions
  af <- setdiff(an, av) ## all functions used
  
  a <- deparse(arg)
  a <- paste0(a, collapse = "") ## lists may somehow produce length > 1 here
  if (substr(a, 1, 5) == "list(") return("list")
  if (a == "NULL") return("NULL")
  ## detection of character arguments is not easy and should not be considered
  ## fool proof since user may pass e.g. a vector of character strings as a 
  ## symbol, which can only really be interpreted as an expression
  if (sum(grep('\\"', a)) && length(setdiff(af, "c")) == 0) return("character")
  
  if (is.data.frame(data)) {
    if (is.symbol(arg) && a %in% names(data)) return("expression")
    if (length(av) == 1L && av %in% names(data)) return("expression")
    e <- eval(arg, envir = data[1:min(nrow(data), 20L), ], 
              enclos = if (is.environment(enclos)) enclos else parent.frame(n + 1L))
    if (inherits(e, "formula")) return("formula")
    if (is.null(e)) return("NULL")
    if (is.list(e)) return("list")
    if (is.character(e) && all(e %in% names(data))) return("character")
    if (is.vector(e) || is.factor(e)) return("expression")
    
    if (recursive && is.language(e)) return(popArgType(e, data = data, n = n + 1L, enclos = enclos))
  }
  
  "expression"
  
}

cutLow <- function(x, breaks, tol =  .Machine$double.eps^0.5) {
  ## a cut function that returns the lower bounds of the cut intervals (as numeric) as levels
  
  breaks <- sort(breaks)
  x <- cut(x + tol, right = FALSE, breaks = breaks, labels = FALSE)
  x <- breaks[-length(breaks)][x]
  x
}




cutLowMerge <- function(x, y, by.x = by, by.y = by, by = NULL, all.x = all, all.y = all, all = FALSE, mid.scales = c("per", "age"), old.nums = TRUE) {
  ## INTENTION: merges y to x by by.x & by.y after cutLow()'ing appropriate
  ## variables in x so that y's values match with x's values
  ## requirements;
  ## * numeric variables in y correspond to lower limits of some intervals OR
  ##   are group variables (e.g. sex = c(0,1))
  ## inputs: two datas as in merge, preferably both data.table, and other args
  ## to merge()
  ## output: a data.table where y has been merged to x after cutLow()
  ## example: merging popmort to a split Lexis object, where popmort's variables
  ## correspond to at least some Lexis time scales
  ## old.nums: return old numeric variable values used in cutLow()'ing?
  ## mid.scales: use mid-point of interval when merging by these Lexis time scales
  ## computed by adding + 0.5*lex.dur, which must exist
  
  if (!is.data.table(x)) {
    stop("x must be a data.table")
  }
  
  if ((is.null(by.x) && !is.null(by.y)) || (!is.null(by.x) && is.null(by.y))) {
    stop("one but not both of by.x / by.y is NULL")
  }
  if (!is.null(by)) by.x <- by.y <- by 
  
  if (length(by.x) != length(by.y)) stop("lengths differ for by.y & by.x")
  all_names_present(x, by.x)
  all_names_present(y, by.y)
  names(by.x) <- by.y
  names(by.y) <- by.x
  
  if (length(mid.scales)>0) all_names_present(x, c("lex.dur", mid.scales))
  
  whScale <- by.x %in% mid.scales
  xScales <- by.x[whScale]
  yScales <- by.y[whScale]
  
  if (length(yScales) > 0) {
    
    oldVals <- copy(with(x, mget(xScales)))
    on.exit(set(x, j = xScales, value = oldVals))
    setattr(oldVals, "names", yScales)
    
    for (yVar in yScales) {
      xVar <- xScales[yVar]
      xBr <- sort(unique(y[[yVar]]))
      xBr <- unique(c(xBr, Inf))
      set(x, j = xVar, value = cutLow(x[[xVar]] + x$lex.dur*0.5, breaks = xBr))
    }
    
  }
  
  ## ensure x retains order (no copy taken of it)
  xKey <- key(x)
  if (length(xKey) == 0) {
    xKey <- makeTempVarName(x, pre = "sort_")
    on.exit(if ("x" %in% ls()) setcolsnull(x, delete = xKey, soft = TRUE), add = TRUE)
    on.exit(if ("z" %in% ls()) setcolsnull(z, delete = xKey, soft = TRUE), add = TRUE)
    x[, (xKey) := 1:.N]
  }
  
  if (any(duplicated(y, by = by.y))) {
    stop("y is duplicated by the inferred/supplied by.y variables (",
         paste0("'", by.y, "'", collapse = ", "), "). ",
         "First ensure this is not so before proceeding.")
  }
  
  ## avoid e.g. using merge.Lexis when x inherits Lexis
  xClass <- class(x)
  on.exit({
    setattr(x, "class", xClass)
    }, add = TRUE)
  setattr(x, "class", c("data.table", "data.frame"))
  
  ## return old numeric values of variables that were cutLow()'d
  ## by keeping them 
  if (old.nums && length(xScales)) {
    tmpXScales <- makeTempVarName(names = c(names(x), names(y)), pre = xScales)
    set(x, j = tmpXScales, value = oldVals)
    on.exit({
      xOrder <- setdiff(names(x), tmpXScales)
      setcolsnull(x, delete = xScales, soft = TRUE)
      setnames(x, tmpXScales, xScales)
      setcolorder(x, xOrder)
      
    }, add = TRUE)
  }
  
  ## merge
  z <- merge(x, y, by.x = by.x, by.y = by.y, 
             all.x = all.x, all.y = all.y, all = all, 
             sort = FALSE)
  
  setDT(z)
  if (old.nums && length(xScales)) {
    ## avoid warning due to coercing double to integer
    set(z, j = xScales, value = NULL)
    setnames(z, tmpXScales, xScales)
  }
  
  zOrder <- intersect(names(x), names(z))
  zOrder <- c(zOrder, setdiff(names(z), names(x)))
  setcolorder(z, zOrder)
  if (length(xKey) > 0) setkeyv(z, xKey)
  z[]
  
}


getOrigin <- function(x) {
  ## input: Date, IDate, or date variable
  ## output: the origin date in Date format,
  ## the origin date being the date where the underlying index is zero.
  if (inherits(x, "Date") || inherits(x, "IDate")) {
    as.Date("1970-01-01")
  } else if (inherits(x, "date")) {
    as.Date("1960-01-01")
  } else if (inherits(x, "dates")) {
    as.Date(paste0(attr(x, "origin"), collapse = "-"), format = "%d-%m-%Y")
  } else {
    stop("class '", class(x), "' not supported; usage of Date recommended - see ?as.Date")
  }
  
}


setcols <- function(x, j, value) {
  ## intention: add new columns to DT via modifying in place, and to DF
  ## via DF$var <- value; both conserve memory (don't take copy of whole data)
  
  if (!is.data.frame(x)) stop("x must be a data.frame")
  if (!is.list(value)) stop("value must be a list of values (columns to add)")
  if (missing(j)) j <- names(value)
  
  if (!is.data.table(x)) {
    x[j] <- value
  } else {
    set(x, j = j, value = value)
  }
  x
}



`%.:%` <- function(x, y) {
  ## INTENTION: hacking formula calls using `:`
  ## which is apparently normally evaluated in C... (in model.matrix.default)
  ## USAGE: e.g. c(1,2) %.:% c(3,4) = c(3, 8)
  ## (instead of getting warning)
  if (length(x) > 1L && length(y) > 1L && is.numeric(x) && is.numeric(y)) {
    return(x*y)
  } else if (length(x) == 1L && length(y) == 1L && is.numeric(x) && is.numeric(y)) {
    return(x:y)
  }
  as.factor(x):as.factor(y)
  
}



RHS2list <- function(formula, handle.adjust=TRUE) {
  ## INTENTION: turns the right-hand side of a formula
  ## into a list of substituted expressions;
  ## each element in list is an expressions separated
  ## by a '+' in the formula. needs to be eval()'d,
  ## preferably using the appropriate data set.
  if (!inherits(formula, "formula")) stop("not a formula")
  
  ## no response
  formula <- formula[c(1, length(formula))]
  
  te <- terms(formula)
  tl <- attr(te, "term.labels")
  
  ## handle adjusting variables (e.g. adjust(V1, V2) -> c("V1", "V2"))
  adj <- tl[substr(tl, 1, 7) == "adjust("]
  if (length(adj) == 0L) adj <- NULL
  if (handle.adjust && !is.null(adj)) {
    tl <- setdiff(tl, adj)
    adjNames <- substr(adj, 8, nchar(adj)-1L)
    adj <- lapply(adj, function(x) parse(text = x))
    adj <- unlist(lapply(adj, eval), recursive = FALSE)
    adj <- lapply(adj, deparse)
    adj <- unlist(adj)
    
    tl <- c(tl, adj)
  }
  
  ## to avoid e.g. c(1,2):c(3,4) NOT evaluating as c(3, 8)
  l <- lapply(tl, function(x) gsub(pattern = ":", x = x, replacement = "%.:%"))
  l <- lapply(l, function(x) parse(text = x)[[1L]])
  
  names(l) <- tl
  
  setattr(l, "adjust", adj)
  
  l
}

RHS2DT <- function(formula, data = data.frame(), enclos = parent.frame(1L)) {
  l <- RHS2list(formula)
  if (length(l) == 0L) return(data.table())
  adj <- attr(l, "adjust")
  
  dana <- names(data)
  dana <- gsub(x=dana, pattern=" %.:% ", replacement = ":")
  dana <- gsub(x=dana, pattern="%.:%", replacement = ":")
  
  ld <- lapply(l, deparse)
  ld <- lapply(ld, function(ch) gsub(x=ch, pattern=" %.:% ", replacement = ":"))
  ld <- lapply(ld, function(ch) gsub(x=ch, pattern="%.:%", replacement = ":"))
  ld <- lapply(ld, function(ch) if (ch %in% dana) which(dana %in% ch) else ch)
  ld <- lapply(ld, function(el) if (is.integer(el)) data[[names(data)[el]]] else NULL)
  ld[which(unlist(lapply(ld, is.null)))] <- NULL
  l[names(ld)] <- ld
  
  l <- lapply(l, function(elem) eval(expr = elem, envir = data, enclos = enclos))
  
  l <- as.data.table(l)
  setattr(l, "adjust", adj)
  l
}

Surv2DT <- function(Surv) {
  sa <- attributes(Surv)
  dt <- copy(Surv)
  setattr(dt, "class", "array")
  dt <- data.table(dt)
  
  type <- attr(Surv, "type")
  statNA <- sum(is.na(dt$status))
  if (statNA) 
    stop("Some status indicators (", statNA  ," values in total) were NA as a result of using Surv(). Usual suspects: original status variable has NA values, or you have numeric status variable with more than two levels and you did not assign e.g. type = 'mstate' (e.g. Surv(time = c(1,1,1), event = c(0,1,2), type = 'mstate') works).")
  
  
  setattr(dt, "type", type)
  testClass <- sa$inputAttributes$time2$class
  if (!is.null(testClass) && testClass == "factor") dt[, status := factor(status, labels = sa$inputAttributes$time2$levels)]
  testClass <- sa$inputAttributes$event$class
  if (!is.null(testClass) && testClass == "factor") dt[, status := factor(status, labels = sa$inputAttributes$event$levels)]
  
  dt[]
}

promptYN <- function(q) {
  
  rl <- readline(prompt = paste0(q, " (Y/N) ::: "))
  y <- c("y", "Y")
  n <- c( "n", "N")
  if (!rl %in% c(y,n)) {
    cat("Answer must be one of the following (without ticks):", paste0("'",c(y, n),"'", collapse = ", "))
    promptYN(q = q)
  }
  
  if (rl %in% y) TRUE else FALSE
  
}

#' @title Adjust Estimates by Categorical Variables
#' @description This function is only intended to be used within a formula
#' when supplied to e.g. \code{\link{survtab_ag}} and should not be
#' used elsewhere. 
#' @param ... variables to adjust by, e.g. \code{adjust(factor(v1), v2, v3)}
#' @return Returns a list of promises of the variables supplied which can be
#' evaluated.
#' @examples 
#' 
#' y ~ x + adjust(z)
#' @export
adjust <- function(...) {
  
  call <- sys.call(1L)
  call <- as.list(call)[1L]
  
  if (deparse(call) %in% c("adjust", "list(adjust)")) stop("Function adjust() only intended to be used within the formulas of certain functions of package popEpi. See e.g. ?survtab_ag for usage.")
  
  mc <- as.list(match.call())[-1L]
  mc
}

evalPopFormula <- function(formula, data = data.frame(), enclos = parent.frame(2L), subset = NULL, Surv.response = TRUE) {
  
  ## INTENTION: given a formula object, returns a DT where each column
  ## is an evaluated expression from formula (separated by  + )
  
  fe <- environment(formula)
  
  either <- FALSE
  if (is.character(Surv.response)) {
    Surv.response <- match.arg(Surv.response, "either")
    Surv.response <- TRUE
    either <- TRUE
  } else if (!is.logical(Surv.response)) {
    stop("Surv.response must be either logical or 'either'")
  }
  
  ## subset if needed ----------------------------------------------------------
  if (!is.null(subset) && !is.logical(subset)) {
    stop("subset must be NULL or a logical vector and not an expression at ",
         "this point. If you see this, complain to the package maintainer.")
  }
  if (!is.null(subset)) {
    keepVars <- c(all.vars(formula), "lex.Xst")
    data <- subsetDTorDF(data, subset = subset, select = keepVars)
  }
  
  
  ## formula -------------------------------------------------------------------
  if (!inherits(formula, "formula")) {
    stop("formula is not of class 'formula'; supply it as e.g. y ~ x")
  }
  if (length(formula) < 3L) {
    stop("formula appears to be one-sided, which is not supported; ",
         "supply it as e.g. y ~ x")
  }
  
  ## response
  y <- eval(formula[[2L]], envir = data, enclos = enclos)
  if (inherits(y, "Surv") && !either && !Surv.response) {
    stop("Response is a result of using Surv(), which is not allowed in ",
         "this context.")
  }
  
  if (!inherits(y, "Surv") && !either && Surv.response) {
    stop("The response of the formula must be a Surv object; ",
         "see ?Surv (in package survival).")
  }
  
  if (inherits(y, "Surv")) {
    y <- Surv2DT(y)
    setcolsnull(y, keep = c("time", "start", "status"), colorder = TRUE)
    if (!any(c("time", "start") %in% names(y))) {
      stop("You must supply function Surv a value to the 'time' ",
           "argument. See ?Surv")
    }
    setnames(y, names(y), c("time", "status")[1:ncol(y)])
  } else {
    y <- data.table(y)
    setnames(y, 1, deparse(formula[[2L]]))
    if (either && inherits(data, "Lexis")) {
      ## we assume the unmodified lex.Xst to be a useful status variable.
      if (!"lex.Xst" %in% names(data)) {
        stop("Supplied a formula without using Surv(), and data was a Lexis ",
             "object, so assumed you intended to use 'lex.Xst' in data as the ",
             "status variable in this context, but that column was missing ",
             "from data.")
      }
      setnames(y, 1, "time")
      y[, "status"] <- data$lex.Xst
      
    }
  }
  
  
  ## RHS
  l <- RHS2DT(formula, data = data, enclos = enclos)
  adj <- attr(l, "adjust")
  
  ## combine
  l <- if (length(l) > 0L) cbind(y, l) else y
  
  setattr(l, "adjust.names", adj)
  setattr(l, "print.names", setdiff(names(l), c(adj, names(y))))
  setattr(l, "Surv.names", names(y))
  setattr(l, "formula", formula)
  
  l
}

oneWhitespace <- function(x) {
  if (!is.character(x)) stop("x not a character")
  x <- paste0(x, collapse = " ")
  while(sum(grep(pattern = "  ", x = x))) {
    x <- gsub(pattern = "  ", replacement = " ", x = x)
  }
  x
}



evalRecursive <- function(arg, env, enc, max.n = 100L) {
  ## INTENTION: digs out actual evaluatable value and expression
  
  if (missing(env)) env <- environment()
  if (missing(enc)) enc <- parent.frame(1L)
  
  if (is.data.frame(env)) {
    na <- names(env)
    env <- env[1:(min(10L, nrow(env))), ]
    
    env <- data.frame(env)
    setattr(env, "names", na)
    
  }
  argSub <- arg
  
  tick <- 1L
  while (!inherits(arg, "try-error") && is.language(arg) && 
         !inherits(arg, "formula") && tick < max.n) {
    
    argSub <- arg
    arg <- try(eval(argSub, envir = env, enclos = enc), silent = TRUE)
    
    tick <- tick + 1L
  }
  
  if (tick == max.n) {
    stop("evaluated expression ", max.n, 
         " times and still could not find underlying expression")
  }
  if (!is.language(argSub)) argSub <- substitute(arg)
  list(arg = arg, argSub = argSub, all.vars = all.vars(argSub))
}


usePopFormula <- function(form = NULL, adjust = NULL, data = data.frame(), 
                          enclos, Surv.response = TRUE) {
  ## INTENTION: evaluates form and combines with adjust appropriately
  ## returns a list of the elements dug out from the formula and adjust 
  ## arguments.
  # formSub <- substitute(form)
  al <- evalRecursive(arg = form, env = data, enc = enclos)
  
  if (!inherits(al$arg, "formula")) stop("'form' is not a formula object")
  
  dt <- evalPopFormula(formula = al$arg, data = data, enclos = enclos, 
                       Surv.response = Surv.response)
  adNames <- attr(dt, "adjust.names")
  prNames <- attr(dt, "print.names")
  suNames <- attr(dt, "Surv.names")
  
  adjust <- evalPopArg(data, adjust, DT = TRUE, recursive = TRUE, 
                       enclos = new.env(), naming = "model",
                       types = c("NULL", "character", "list", "expression"))
  
  if (is.data.frame(adjust) && (nrow(adjust) == 0L || ncol(adjust) == 0L)) {
    stop("adjust evaluated to an empty data.frame")
  }
  if (!is.null(adjust) && ncol(adjust) > 0L && length(adNames) > 0L) {
    stop("Cannot both use argument 'adjust' AND use an adjust() term within ",
         "the formula argument. Please only use one.")
  }
  if (is.null(adjust) && length(adNames) > 0L) {
    adjust <- dt[, .SD, .SDcols = c(adNames)]
  }
  
  
  print <- NULL
  if (length(prNames > 0L)) print <- dt[, .SD, .SDcols = eval(prNames)]
  
  list(y = dt[, .SD, .SDcols = c(suNames)], 
       print = print, 
       adjust = adjust, formula = al$arg)
}


aliased_cols <- function(data, cols) {
  
  if (missing(cols)) cols <- names(data)
  all_names_present(data, cols)
  
  if (length(cols) < 2L) return(invisible())
  
  x <- with(data, mget(cols))
  x <- lapply(x, duplicated)
  
  sub_cols <- cols
  tl <- list()
  ## loop: each step reduce vector of names by one
  ## to avoid testing the same variables twice (in both directions)
  tick <- 0L
  aliased <- FALSE
  while (!aliased && length(sub_cols) > 1L && tick <= length(cols)) {
    
    currVar <- sub_cols[1L]
    sub_cols <- setdiff(sub_cols, currVar)
    tl[[currVar]] <- unlist(lapply(x[sub_cols], function(j) identical(x[[currVar]], j)))
    aliased <- sum(tl[[currVar]])
    
    tick <- tick + 1L
  }
  
  if (tick == length(cols)) warning("while loop went over the number of columns argument cols")
  
  ## result: list of logical vectors indicating if a column is aliased
  ## with other columns
  tl[vapply(tl, function(j) sum(j) == 0L, logical(1))] <- NULL
  
  if (length(tl) == 0L) return(invisible())
  
  ## take first vector for reporting
  var <- names(tl)[1L]
  aliases <- names(tl[[1L]])[tl[[1]]]
  aliases <- paste0("'", aliases, "'", collapse = ", ")
  stop("Variable '", var, "' is aliased with following variable(s): ", aliases, ".")
  
  invisible()
}



