#######################################################################################################
#               Written by Andrew R Sommerlot <andrewrs@vt.edu>, March 2016                           #
#######################################################################################################
#' Gets met data from the wundergorund 10 day forecast and outputs SWAT IO format meteorological input files
#' @param centroids - data.frame object or location of csv file with two columns, the first being lattitude and second being longidute in decimal degrees. These should be the centroids of the swat model subbasins and are the locations of the ouput met data.
#' @param outDir - Directory where ouput files will be saved.
#' @param apiKey - String input of your freely available api key from wunderground.com. See https://www.wunderground.com/weather/api/d/pricing.html for more details. Select the ANVIL version of the developer plan, record the resulting key and use it in this function.
#' @return returns wudnerground 10 day forecast met data in swat IO format
#' @examples
#' \dontrun{
#' centroids = data.frame(lat = 38, lon = 79)
#' outDir = "test"
#' getSWATwunderForecast(centroids=centroids, outDir=outDir)
#' }
#' @importFrom jsonlite fromJSON
#' @importFrom utils read.csv write.table URLencode
#' @export

getSWATwunderForecast <- function(centroids, outDir = getwd(), apiKey){
  retDir = getwd()
  setwd(outDir)

  if(class(centroids) == 'character'){
    cent <- read.csv(centroids)
  } else if(class(centroids) == 'data.frame'){
    cent <- centroids
  }

  ln <- nrow(cent)

  forecast <- list()
  for(i in 1:ln){
    url <- paste('http://api.wunderground.com/api/', apiKey, '/conditions/forecast10day/q/',cent[i,1], ',', cent[i,2],'.json', sep = '')
    lookup <- URLencode(url)
    readin <- readLines(lookup, warn = "F")
    data <- fromJSON(readin)
    forecast[[i]] <- data[[3]][[2]][[1]]
  }

  # initialize lists
  temp = list()
  rhum = list()
  pcp = list()
  wnd = list()

  ## format the date
  date.days <- formatC(forecast[[1]]$date$yday, width = 3, flag = 0)
  date.year <- formatC(forecast[[1]]$date$year, width = 4)
  fulldate <- paste(date.year, date.days, sep = '')

  for(i in 1:length(forecast)){


    # get the info
    tmax <- formatC((as.numeric(forecast[[i]]$high$fahrenheit) - 32)*(5/9), digits = 1, width = 5, flag = '0', format = 'f')
    tmin <- formatC((as.numeric(forecast[[i]]$low$fahrenheit) - 32)*(5/9), digits = 1, width = 5, flag = '0', format = 'f')
    temp[[i]] <- paste(tmax, tmin, sep = '')

    # relative humididty
    rhum[[i]] <- formatC(forecast[[i]]$avehumidity, digits = 3, width = 8, flag = '0', format ='f')

    ##for pcp # !!! might need to add snow dont know yet
    pcp[[i]] <- formatC(forecast[[i]]$qpf_allday[,1]*25.4, digits = 1, width = 5, flag = '0', format ='f')

    ##for wnd
    wnd[[i]] <- formatC(forecast[[i]]$avewind$kph*(1000/3600), digits = 3, width = 8, flag = '0', format = 'f')

   ################################################################################
  }

  tmpprint = data.frame(fulldate, as.data.frame(temp))
  pcpprint = data.frame(fulldate, as.data.frame(pcp))
  rhumprint = data.frame(fulldate, as.data.frame(rhum))
  wndprint = data.frame(fulldate, as.data.frame(wnd))

  Lati <- cent[,1]
  Long <- cent[,2]

  headtmp <- c('Temp	Source: wunderground 10 day forecast api', paste('Lati', paste(Lati, collapse = '   '), sep = '   '), paste('Long', paste(Long, collapse = '   '), sep = '   '), 'Elev')
  headpcp <- c('Precip Source: wunderground 10 day forecast api', paste('Lati', paste(Lati, collapse = '   '), sep = '   '), paste('Long', paste(Long, collapse = '   '), sep = '   '), 'Elev')
  headhmd <- 'Relative Humidity % 	Source: wunderground 10 day forecast api'
  headwnd <- 'Wind Speed m/s 	Source: wunderground 10 day forecast api'

  ############################################################################################################
  #write tmp and pcp header files and append the data
  writeLines(headtmp, 'tmp1.tmp')
  writeLines(headpcp, 'pcp1.pcp')
  writeLines(headhmd, 'hmd.hmd')
  writeLines(headwnd, 'wnd.wnd')
  write.table(tmpprint, 'tmp1.tmp', append = TRUE, sep = '', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(pcpprint, 'pcp1.pcp', append = TRUE, sep = '', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(rhumprint, 'hmd.hmd', append = TRUE, sep = '', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(wndprint, 'wnd.wnd', append = TRUE, sep = '', row.names = FALSE, col.names = FALSE, quote = FALSE)

  setwd(retDir)
}

