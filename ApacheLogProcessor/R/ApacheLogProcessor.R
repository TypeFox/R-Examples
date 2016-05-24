
#' read.apache.log.combined
#'
#' Reads the Apache Log Combined Format and return a data frame with the log data.
#'
#' The functions recives a full path to the log file and process the default log combined format of Apache. 
#' LogFormat "\%h \%l \%u \%t \\"\%r\\" \%>s \%b \\"\%\{Referer\}i\\" \\"\%\{User-Agent\}i\\"" combined
#'
#' @param file string. Full path to the log file.
#' @param url_includes regex. If passed only the urls that matches with the regular expression passed will be returned.
#' @param url_excludes regex. If passed only the urls that don't matches with the regular expression passed will be returned.
#' @param url_query_string boolean. If the query string will be included in url column.
#' @param url_http_version booelan. If the http version will be included in url column.
#' @param url_http_method boolean. If the http method will be included in url column.
#' @param columns list. List of columns names that will be included in data frame output. All columns is the default value. c("ip", "datetime", "url", "httpcode", "size" , "referer", "useragent")
#' @param num_cores number. Number of cores for parallel execution, if not passed 1 core is assumed. Used only to convert datetime form string to datetime type.
#' @param prevent_quotes_inside_url boolean. Search and remove the quotes inside the url field. 
#' @return a data frame with the apache log file information.
#' @author Diogo Silveira Mendonca
#' @seealso \url{http://httpd.apache.org/docs/1.3/logs.html}
#' @examples
#' path = system.file("examples", "apache_example.txt", package = "ApacheLogProcessor")
#' 
#' #Read the full log with all lines and columns and return a data frame
#' df1 = read.apache.log.combined(path)
#' 
#' #Read only the lines that url matches with the pattern passed
#' df2 = read.apache.log.combined(path, url_includes="infinance")
#' 
#' #Read only the lines that url matches with the pattern passed, but do not matche the exclude pattern
#' df3 = read.apache.log.combined(path, url_includes="infinance", url_excludes="infinanceclient")
#' 
#' #Removes the method and query string from urls
#' df4 = read.apache.log.combined(path, url_http_method=FALSE, url_query_string=FALSE)
#' 
#' #Return only the ip, url and datetime columns
#' df5 = read.apache.log.combined(path, columns=c("ip", "url", "datetime"))
#' 
#' #Process using 2 cores in parallel for speed up. 
#' df6 = read.apache.log.combined(path, num_cores=2)
#' 
#' @import foreach
#' @import parallel
#' @import doParallel
#' @export
read.apache.log.combined <- function(file, url_includes = "", url_excludes = "", 
                                     url_query_string = TRUE, url_http_version = TRUE, url_http_method = TRUE,
                                     columns = c("ip", "datetime", "url", "httpcode", "size" , 
                                                 "referer", "useragent"), num_cores = 1,
                                     prevent_quotes_inside_url = FALSE){
  
  #=== REMOVE QUOTES INSIDE QUOTES IN URL FIELD ===================================
  if(prevent_quotes_inside_url == TRUE){
    text <- readLines(file)
    text <- gsub("\\\\\"", "'", text)
    tConnection <- textConnection(text)
  }else{
    tConnection <- file
  }
  #=== LOAD THE APACHE ACCESS LOG FILE AS CSV =====================================
  
  logDf = read.csv(tConnection, header = FALSE, sep = " ", quote = "\"",
                   dec = ".", fill = FALSE, stringsAsFactors = FALSE)
  
  if (prevent_quotes_inside_url == TRUE){
    close(tConnection)
  }
  
  #=== SET UP THE COLUMNS =========================================================
  
  #remove the columns that will not be used
  logDf$V2 <- NULL; 
  logDf$V3 <- NULL; 
  
  #set the column names
  colnames(logDf) <- c("ip", "datetime", "timezone", "url", "httpcode", "size" , "referer", "useragent")
  
  #include only the columns required
  c_include = c()
  for (col in colnames(logDf)){
    if (col %in% columns){
      c_include <- c(c_include, col)
      if(col == "datetime"){
        c_include <- c(c_include, "timezone")
      }
    }
  }
  logDf <- logDf[,c_include]
  
  #=== APPLY RULES FROM LINES ====================================================
  
  #filter the lines to be included
  line_numbers <- grep(url_includes, t(logDf["url"]))
  
  #filter the lines to be excluded
  if(url_excludes != ""){
    line_numbers <- 
      line_numbers[ !line_numbers %in% 
                              grep(url_excludes, t(logDf["url"]))]
  }
  
  #Get only the necessary lines
  logDf <- logDf[line_numbers,]
  
  #=== CLEAR THE DATETIME AND TIMEZONE COLUMNS ====================================
  if ("datetime" %in% c_include){
    
    #Create a vector of dates
    dates = seq( as.POSIXlt(Sys.Date()), by=1, len=nrow(logDf))
    lct <- Sys.getlocale("LC_TIME") 
    Sys.setlocale("LC_TIME", "C")
    
    #CREATE CLUSTERS FOR PARALLEL EXECUTION
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    
    #parse the dates form data frame to dates vector
    dates <- foreach (i = 1:nrow(logDf), .combine=rbind) %dopar%{
      Sys.setlocale("LC_TIME", "C")
      datetimeWithTimezone = paste(logDf$datetime[i], logDf$timezone[i])
      dates[i] <- strptime(datetimeWithTimezone, format="[%d/%b/%Y:%H:%M:%S %z]")
      return(dates[i])
    }
    
    Sys.setlocale("LC_TIME", lct)
    dates <- as.POSIXct(dates, origin="1970-01-01")
    
    #Shutdown the cluster
    stopCluster(cl)
    
    #Create a new data frame with coverted dates
    datesFrame <- data.frame(dates)
    colnames(datesFrame) <- c("datetime")
    
    #Remove old date and timezone columns
    logDf$datetime <- NULL; 
    logDf$timezone <- NULL;
    
    #Inserts the converted dates in logDf data frame
    logDf <- cbind(logDf, datesFrame)
  }
  
  
  #==== CLEAR THE URL FIELD ===================================================
  if("url" %in% c_include){
      
      for(i in 1:nrow(logDf)){
        
        #removes the query string
        if (!url_query_string){
          logDf$url[i] <- sub("\\?.* ", " ", logDf$url[i])
        }
        #removes the http version
        if(!url_http_version){
          logDf$url[i] <- sub(" HTTP/.*", "", logDf$url[i])
        }
        #removes the url method
        if (!url_http_method){
          logDf$url[i] <- sub("[A-Z]* ", "", logDf$url[i])
        }
      }
      
  }
  
  #=== CONVERT THE SIZE COLUMN FROM TEXT TO NUMERIC ===========================
  
  if("size" %in% c_include){
    sizes <- logDf$size
    sizes <- as.numeric(sizes)
    logDf$size <- NULL
    sizesFrame <- data.frame(sizes)
    colnames(sizesFrame) <- c("size")
    logDf <- cbind(logDf, sizesFrame)
  }
  
  #=== RETURN THE DATA FRAME ==================================================
  logDf
  
}

#' Apache log combined file example.
#'
#' A set of 12 log lines in Apache Log Combined Format
#'
#' @format LogFormat "\%h \%l \%u \%t \\"\%r\\" \%>s \%b \\"\%\{Referer\}i\\" \\"\%\{User-Agent\}i\\"" combined
#' @source \url{http://www.infinance.com.br/}
#' @name apache_example
NULL
