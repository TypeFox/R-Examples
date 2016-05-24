#  ------------------------------------------------------------------------

#' From CliFlo to \pkg{clifro}: Enhancing The National Climate Database With \R
#'
#' Import data from New Zealand's National Climate Database via CliFlo into \R
#' for exploring, analysis, plotting, exporting to KML, CSV, or other software.
#'
#' The \pkg{clifro} package is intended to simplify the process of data
#' extraction, formatting and visualisation from the
#' \href{http://cliflo.niwa.co.nz/}{CliFlo web portal}. It
#' requires the user to build a query consisting of 3 main components; the user,
#' the datatype(s) and the station(s). These are
#' then combined using the \code{\link{cf_query}} function that sends the query
#' to the CliFlo database and returns the results that can easily be plotted
#' using generic plotting functions.
#'
#' This package requires the user to already have a current subscription to the
#' National Climate Database unless a public user is sought, where data is
#' limited to Reefton Ews. Subscription is free and can obtained from
#' \url{http://cliflo.niwa.co.nz/pls/niwp/wsubform.intro}.
#'
#' @seealso \code{\link{cf_user}}, \code{\link{cf_datatype}}, and
#'   \code{\link{cf_station}} for choosing the clifro user, datatypes and
#'   stations, respectively.
#' @name clifro
#' @aliases clifro-package
#' @docType package
#' @keywords package
#' @examples
#' \dontrun{
#' # Create a public user ----------------------------------------------------
#'
#' public.user = cf_user() # Defaults to "public"
#' public.user
#'
#' # Select datatypes --------------------------------------------------------
#'
#' # 9am Surface wind (m/s)
#' wind.dt = cf_datatype(2, 1, 4, 1)
#'
#' # Daily Rain
#' rain.dt = cf_datatype(3, 1, 1)
#'
#' # Daily temperature extremes
#' temp.dt = cf_datatype(4, 2, 2)
#'
#' # Combine them together
#' all.dts = wind.dt + rain.dt + temp.dt
#' all.dts
#'
#' # Select the Reefton Ews station ------------------------------------------
#'
#' reefton.st = cf_station()
#' reefton.st
#'
#' # Submit the query --------------------------------------------------------
#'
#' # Retrieve all data from ~ six months ago at 9am
#' reefton.data = cf_query(public.user, all.dts, reefton.st,
#'                         paste(as.Date(Sys.time()) - 182, "9"))
#' reefton.data
#'
#'
#' # Plot the data -----------------------------------------------------------
#'
#' # Plot the 9am surface wind data (first dataframe in the list) ---
#' reefton.data[1]
#'
#' # all identical - although passed to different methods
#' plot(reefton.data)    #plot,cfDataList,missing-method
#' plot(reefton.data, 1) #plot,cfDataList,numeric-method
#' plot(reefton.data[1]) #plot,cfData,missing-method --> plot,cfWind,missing-method
#'
#' speed_plot(reefton.data)
#' direction_plot(reefton.data)
#'
#' # Plot the daily rain data (second dataframe in the list) ---
#' reefton.data[2]
#'
#' # With runoff and soil deficit
#' plot(reefton.data, 2)
#'
#' # Just plot amount of rain (mm)
#' plot(reefton.data, 2, include_runoff = FALSE)
#'
#' # Plot the hourly temperature data (third dataframe in the list) ---
#' plot(reefton.data, 3)
#'
#' # Pass an argument to ggplot2::theme
#' library(ggplot2) # for element_text()
#' plot(reefton.data, 3, text = element_text(size = 18))
#' }
NULL

# Validation (internals) --------------------------------------------------

#' Validation Functions For The \code{cfUser} Class
#'
#' These internal functions are used by the \code{\link{cf_user}} constructor
#' function to ensure the user has a valid subscription to CliFlo.
#'
#' \code{cf_login} initiates a curl handle storing the cookies in the current
#' \R session's temporary directory. It then POSTs the user credentials to the
#' CliFlo login page and stores the resultant \code{h1} heading to check for the
#' string 'Info'. The cookies are kept for future (immediate) use.
#'
#' \code{cf_logout} points the curl handle to the existing cookie session
#' initiated with \code{cf_login}. It reads the header information from the
#' cliflo logout page to ensure no HTTP error and logs the user out on
#' cliflo and deletes the cookies. This should be (is) called immediately after
#' \code{cf_login} in any function requiring a login, using
#' \code{\link{on.exit}} to ensure the user isn't still logged in on the server,
#' after the function call, for any reason.
#'
#' \code{valid_cfuser} is the validation function for the \code{cfUser} class
#' and uses  \code{cf_login} to ensure the credentials are authenticated on the
#' CliFlo server and then (\code{cf_})logs out immediately afterwards. It also
#' ensures the user provides exactly one username and password - except for
#' 'public' users.
#'
#' @param object S4 object which inherits the \code{cfUser} class
#'
#'@param msg Display a 'successful logout' message, defaults to
#' \code{TRUE}.
#'
#' @importFrom RCurl getCurlHandle postForm getURL
#' @importFrom XML htmlParse xmlValue
#' @importFrom selectr querySelector
#' @keywords internal
#' @aliases cf_logout cf_login
#' @name valid_cfuser
#' @rdname valid_cfuser
#' @examples
#' \dontrun{
#' cf_user("public")                    # Returns a valid object
#' cf_user("bad_name", "bad_password")    # Bad Login
#' }

cf_login = function(object){
  cookies = file.path(tempdir(), object@username)
  curl = getCurlHandle(followlocation = TRUE,
                       cookiejar = cookies,
                       cookiefile = cookies,
                       useragent = paste("clifro", R.Version()$version.string),
                       timeout = 100)
  if (object@username == "public"){
    login_html = htmlParse(getURL(
      "http://cliflo.niwa.co.nz/pls/niwp/wgenf.genform1",
      curl = curl
    ))
    result = "Info"
  }
  else{
    login_html = htmlParse(postForm(
      "http://cliflo.niwa.co.nz/pls/niwp/wa.logindb",
      cusername = object@username,
      cpwd = rot(object@password, 3),
      ispopup = "false",
      submit = "login",
      curl = curl))
    result = xmlValue(querySelector(login_html, "h1"))
  }
  rm(curl)
  gc()
  return(grepl("Info", result))
}

#' @rdname valid_cfuser
#' @importFrom RCurl getCurlHandle getURLContent getURL
cf_logout = function(object, msg = TRUE){
  cookies = file.path(tempdir(), object@username)
  curl = getCurlHandle(followlocation = TRUE,
                       timeout = 100,
                       useragent =
                         paste("clifro", R.Version()$version.string),
                       cookiefile = cookies,
                       cookiejar = cookies)

  header = getURLContent("http://cliflo.niwa.co.nz/pls/niwp/wa.logout",
                         curl = curl, header = TRUE)
  if (!grepl("OK", header$header[11]))
    stop("HTTP error")

  getURL("http://cliflo.niwa.co.nz/pls/niwp/wa.logout", curl = curl)

  file.remove(cookies)
  if (msg)
    message("Logout successful")
}

#' @rdname valid_cfuser
valid_cfuser = function(object){
  length_username = length(object@username)
  length_password = length(object@password)
  errors = character()

  if (length_username != 1){
    msg = "Exactly one username must be specified"
    errors = c(errors, msg)
  }

  if (tolower(object@username) != "public" && length_password != 1){
    msg = "Exactly one password must be specified"
    errors = c(errors, msg)
  }

  if (tolower(object@username) != "public"){
    login_OK = cf_login(object)
    if (login_OK)
      on.exit(cf_logout(object, msg = FALSE))
  } else
    login_OK = TRUE

  if (!login_OK){
    msg = "Bad Login"
    errors = c(errors, msg)
  }

  if (length(errors) == 0)
    TRUE
  else
    errors
}

# cfUser class ------------------------------------------------------------

#' @rdname cfUser-class
#' @name cfUser-class
#' @aliases cf_user
#' @importFrom methods setClass
setClass("cfUser",
         representation = representation(username = "character",
                                         password = "character"),
         validity = valid_cfuser)

#' The Clifro User Object
#'
#' Create a \code{cfUser} object to allow the user to log into CliFlo from \R
#' and  build their query.
#'
#' An object inheriting from the \code{cfUser} class is created by the constructor
#' function \code{cf_user}. The user must have an active subscription to cliflo
#' in order to create a valid object, unless a 'public' user is sought.
#' Visit \url{http://cliflo.niwa.co.nz/} for more information and to subscribe
#' to cliflo.
#'
#' @param username a character string to be used as the cliflo username
#' @param password a character string to be used as the cliflo password
#'
#' @note For the 'public' user (see examples) only the Reefton Ews station data
#' is available.
#'
#' @importFrom methods new
#' @rdname cfUser-class
#' @name cfUser-class
#' @aliases cfUser
#' @aliases cfUser-class
#' @return \code{cfUser} object
#' @export
#' @seealso \code{\link{valid_cfuser}} for details on the validation of
#' \code{cfUser} and \code{\link{summary,cfUser-method}} to summarise user
#' information.
#' @examples
#' \dontrun{
#' public.cfuser = cf_user(username = "public")
#' public.cfuser
#' }
cf_user = function(username = "public", password = character()){
  new("cfUser", username = username, password = password)
}

# Initialize the cfUser with a cryptic password
#' @importFrom methods setMethod validObject
setMethod("initialize", "cfUser", function(.Object, username, password){
  .Object@username = username
  if (length(password) == 1){
    .Object@password = rot(password, 60)
  }
  else
    .Object@password = password

  validObject(.Object)
  return(.Object)
})

# Methods -----------------------------------------------------------------

#'@importFrom methods setGeneric
if (!isGeneric("summary"))
  setGeneric("summary", function(object, ...) standardGeneric("summary"))

#' Summarise User Information
#'
#' Show the subscription status for the \pkg{clifro} user
#'
#' @param object an object of class \code{cfUser}.
#'
#' @importFrom RCurl getCurlHandle getForm
#' @importFrom selectr querySelectorAll querySelector
#' @importFrom XML htmlParse xmlValue
#' @importFrom lubridate dmy round_date now with_tz
#' @aliases summary,cfUser-method
#' @export
setMethod("summary", signature(object = "cfUser"),
          function(object){
  if (object@username == "public")
    return(object)
  cf_login(object)
  on.exit(cf_logout(object, msg = FALSE))
  cookies = file.path(tempdir(), object@username)
  curl = getCurlHandle(followlocation = TRUE,
                       timeout = 100,
                       useragent =
                         paste("clifro", R.Version()$version.string),
                       cookiefile = cookies,
                       cookiejar = cookies)
  user_info_xml =
    getForm("http://cliflo.niwa.co.nz/pls/niwp/wa.subscr_info",
            sub = "t",
            curl = curl)
  user_info_html = querySelectorAll(htmlParse(user_info_xml),
                                    "body.popup > div")
  info = gsub("  |   |    |     ", " ", sapply(user_info_html, xmlValue))
  rows = sapply(querySelectorAll(user_info_html[[3]], "b"), xmlValue)
  rows = gsub(",", "", rows)
  subscription_level = xmlValue(querySelector(user_info_html[[5]], "b"))
  expiry = strsplit(info[2], ": ")[[1]][2]
  rows_used = as.numeric(rows[1])
  total_rows = as.numeric(rows[3])
  time_diff = dmy(expiry, tz = "Pacific/Auckland") -
    with_tz(round_date(now(), "month"), "Pacific/Auckland")
  cat(paste0("Username is: ", object@username, "\n",
               "Subscription status:\n\n",
               "Your subscription expires on: ", expiry, " (", format(time_diff),
               ")\n", "You have used ",
               format(rows_used, big.mark = ",", scientific = FALSE),
               " rows (", round(rows_used / total_rows * 100, 1), "%)\n",
               "from a subscription total of ",
               format(total_rows, big.mark = ",", scientific = FALSE), " rows.\n",
               "Remaining rows: ",
               format(total_rows - rows_used, big.mark = ",",
                      scientific = FALSE), ".\n",
               "Your subscription level is: ", subscription_level, "\n"))
})


# Show
#' @importFrom methods setMethod
setMethod("show", "cfUser", function(object){
  status = "Authenticated clifro User\n"
  if (tolower(object@username) == "public")
    message("public user - only data from Reefton Ews (3925) available")
  else
    cat(paste0(status, "Username is: ", object@username, "\n"))
})

## Internal function to hide password
rot <- function(ch, k) {
  p0 <- function(...) paste(c(...), collapse = "")
  A <- c(letters, LETTERS, " '", paste(0:9))
  I <- seq_len(k)
  chartr(p0(A), p0(c(A[-I], A[I])), ch)
}
