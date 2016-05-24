#' Search for and collect BISON data.
#' 
#' @importFrom plyr compact llply ldply rbind.fill
#' @param input Output from bison function.
#' @param datatype One of counties, states, data_df, data_list, or NULL.
#' @return A data.frame or list.
#' @description datatype = data_df returns all data in a data.frame, while data_list
#'    returns data in a list. If you are getting data from a call to 
#'    \code{\link{bison_solr}} you must set datatype="data_df" - other options will 
#'    throw an error. 
#' @examples \donttest{
#' # output data
#' out <- bison(species="Bison bison", type="scientific_name", count=10)
#' class(out) # check right class
#' bison_data(out) # summary of output
#' bison_data(input=out, datatype="data_df") # as a data.frame
#' bison_data(input=out, datatype="data_list") # as a list
#' bison_data(input=out, datatype="counties") # summary by counties
#' bison_data(input=out, datatype="states") # summary of states
#' 
#' # SOLR occurrences endpoint
#' out <- bison_solr(scientific_name='"Ursus americanus"')
#' bison_data(input=out)
#' }
#' @export
bison_data <- function(input = NULL, datatype=NULL) UseMethod("bison_data")

#' @S3method bison_data bison
bison_data.bison <- function(input = NULL, datatype=NULL)
{
  if(!is.bison(input))
    stop("Input is not of class bison")
  
  if(is.null(datatype)){
    data.frame(c(input[1], input$occurrences$legend))
  } else
    if(datatype=="counties"){
      if(class(input$counties$data[[1]])=="character"){
        df <- ldply(input$counties$data)
      } else
      {
        df <- ldply(input$counties$data, function(x) data.frame(x))
      }
      names(df)[c(1,3)] <- c("record_id","county_name")
      return(df)
    } else
      if(datatype=="states"){
        df <- ldply(input$states$data, function(x) data.frame(x))
        names(df)[c(1,3)] <- c("record_id","county_fips")
        return(df)
      } else
        if(datatype=="data_df"){
          withlatlong <- input$data[sapply(input$data, length, USE.NAMES=FALSE) == 8]
          data_out <- ldply(withlatlong, function(x){
            x[sapply(x, is.null)] <- NA
            data.frame(x[c("id","name","longitude","latitude","provider")])
          })
          data_out$longitude <- as.numeric(as.character(data_out$longitude))
          data_out$latitude <- as.numeric(as.character(data_out$latitude))
          return(data_out)
        } else
          if(datatype=="data_list"){
            withlatlong <- input$data[sapply(input$data, length, USE.NAMES=FALSE) == 8]
            data_out <- llply(withlatlong, function(x){
              x[sapply(x, is.null)] <- NA
              temp <- x[c("id","name","longitude","latitude","provider")]
              temp$longitude <- as.numeric(as.character(temp$longitude))
              temp$latitude <- as.numeric(as.character(temp$latitude))              
              temp
            })
            return(data_out)
          }
}

#' @S3method bison_data bison_solr
bison_data.bison_solr <- function(input = NULL, datatype="data_df")
{
  num_found=rec=high=fac=NULL
  
  if(!is.bison_solr(input))
    stop("Input is not of class bison_solr")
  if(!datatype=="data_df")
    stop("datatype must equal data_df")
  if(!length(input$num_found) == 0)
    num_found <- input$num_found
  if(!length(input$records) == 0)
    rec <- do.call(rbind.fill, lapply(input$records, data.frame, stringsAsFactors=FALSE))
  if(!length(input$highlight) == 0)
    high <- do.call(rbind.fill, lapply(input$highlight, data.frame, stringsAsFactors=FALSE))
  if(!length(input$facets) == 0)
    fac <- lapply(input$facets$facet_fields, function(x){
      data.frame(do.call(rbind, lapply(seq(1, length(x), by=2), function(y){
        x[c(y, y+1)]
      })), stringsAsFactors=FALSE)
    })
  list(num_found = num_found, records = rec, highlight = high, facets = fac)
}