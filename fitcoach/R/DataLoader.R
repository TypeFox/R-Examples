#' R6 class for Loading Fitbit data 
#' 
#' DataLoader is an R6 Class that connects to the Fitbit API with the  credentials, requests the data, and writes the response to JSON files, 
#'
#' @docType class
#' @format A \code{\link{R6Class}} generator object
#' @keywords data
#' 
#' @importFrom R6 R6Class
#' @importFrom httr content
#' @export DataLoader
#' 
#' @section Methods:
#' \describe{
#'   \item{\code{connect(appname, key, secret, cache.file)}}{This method connects to the Fitbit API and to your application. 
#'   \cr \code{appname}: Name of the Fitbit App
#'   \cr \code{key}: Fitbit API Client key
#'   \cr \code{secret}: Fibit API Client secret
#'   \cr \code{cache.file}: Path to a cached token file, instead of providing credentials in the function call}
#'   \item{\code{request(type = "day", activities = "", start.date = Sys.Date(), end.date = "", path = "./json/"))}}{This method builds the request URLs, sends the requests and writes response to JSON files, in the specified folder.
#'   \cr \code{type}: Type of time series. Must be 'day' or 'intraday'.
#'   \cr \code{activities}: A list of the Fitibit activities to be retrieved.
#'   \cr \code{start.date}: Start date in format YYYY-mm-dd.
#'   \cr \code{end.date}: End date in format YYYY-mm-dd.
#'   \cr \code{path}: Folder where the JSON files will be written.}
#' }
#' 
#' @examples \dontrun{
#' testObject <- DataLoader$new()
#' 
#' testObject$connect(appname = "abcd",
#'                    key = "123ABC",
#'                    secret = "3089e3h1ac9dde0aa67b54ajc8691j44")
#' 
#' testObject$request(
#'     type = 'day', 
#'     activities = list("calories", "steps", "distance", "minutesVeryActive"), 
#'     start.date = "2016-01-01", 
#'     end.date = "2016-02-01", 
#'     path = "~/fitbit-daily/")
#' }


##
## Begin lassence@ code
##

DataLoader <- R6::R6Class (
    "DataLoader",
    
    public = list (
        
        ### Public variables

        # API Token
        api.token = NA,
        # Request response
        response = NA,
        
        ### METHOD initialize
        ### Standard R6 Initialize function

        initialize = function () {
            message("Object DataLoader initialized")
        },
        
        ### METHOD connect
        ### Connects to the API with credentials

        connect = function (appname, key, secret, cache.file) {
            
            # If cache file provided, use it
            if(!missing(cache.file)) {
                self$api.token <- readRDS(cache.file)[[1]]

            # Else, check if cache file exists    
            } else {
                if (file.exists('.httr-oauth')) {
                    if (difftime(Sys.time(), file.info('.httr-oauth')$mtime, units = "mins") < 60) {
                        self$api.token <- readRDS('.httr-oauth')[[1]]
                    } else {
                        # Known bug: autorefresh does not work in basic mode
                        # https://github.com/hadley/httr/pull/320
                        file.remove('.httr-oauth')
                        self$api.token <- connectToAPI(appname, key, secret)
                    }
                } else {
                    self$api.token <- connectToAPI(appname, key, secret)
                }
            }
        },
        
        ### METHOD request
        ### Build URL, send request and write response to JSON file

        request = function (type = "day",
                           activities = "",
                           start.date = Sys.Date(),
                           end.date = "",
                           path = "./json/") {
         
            # Check 'type' argument
            if (!(type %in% c("day", "intraday")))
                stop("Invalid 'type'. Must be 'day' or 'intraday'")
            
            # Check 'start.date' argument
            if (!(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", start.date)))
                stop("Invalid 'start.date'. Must be in the following format: 'YYYY-MM-dd'")
            
            # Call request function for each activity
            for (acty in activities) {

                self$response <- makeAPIRequest(
                    type = type,
                    activity = acty,
                    start.date = start.date,
                    end.date = end.date,
                    api.token = self$api.token
                )
                
                writeToJSON(content = httr::content(self$response, as = "text"),
                            path = path,
                            type = type, 
                            activity = acty,
                            start.date = start.date)
                
            }
        }
    )
)

##
## End lassence@ code
##

