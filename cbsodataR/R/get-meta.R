#' Get metadata of a cbs table
#' 
#' Retrieve the meta data of a CBS open data table.  Caching (\code{cache=TRUE}) may improve
#' the performance considerably. The meta data of a CBS table is 
#' @param id internal id of CBS table, can be retrieved with \code{\link{get_table_list}}
#' @param cache should the result be cached? C
#' @param base_url optionally specify a different server. Useful for
#' third party data services implementing the same protocal.
#' @return cbs_table object containing serveral \code{data.frames} with meta data.  
#' @importFrom whisker whisker.render
#' @importFrom jsonlite fromJSON
#' 
#' @export
get_meta <- function(id, cache=FALSE, base_url = CBSOPENDATA){
  url <- whisker.render("{{BASEURL}}/{{API}}/{{id}}"
                       , list( BASEURL = base_url
                             , API = API
                             , id = id
                             )
                       )
  meta_top <- resolve_resource( url
                              , "Retrieving meta data for table '", id, "'"
                              , cache = cache
                              )
  
  #exclude data sets
  meta_idx = !(meta_top$name %in% c("TypedDataSet", "UntypedDataSet"))
  meta_top <- meta_top[meta_idx,]
  
  meta <- lapply(meta_top$url, resolve_resource, "Retrieving ", url, cache=cache)
  names(meta) <- meta_top$name
  structure( meta
           , class="cbs_table"
           )
}

#' Load meta data from a downloaded table
#' 
#' @param dir Directory where data was downloaded
#' @return cbs_table object with meta data
#' @export
get_meta_from_dir <- function(dir){
  wd <- setwd(dir)
  on.exit(setwd(wd))
  
  meta_files <- list.files(".", "*.csv")
  meta_files <- meta_files[meta_files != "data.csv"]
      
  meta <- lapply(meta_files, read.csv, colClasses="character")
  names(meta) <- sub("\\.csv$","", meta_files)
  meta$directory <- dir
  structure( meta
             , class="cbs_table"
  )
}
### testing


#dump_meta("81819NED")
#get_meta_from_dir("D:/data/StatLine/00370/")

### testing
#meta <- get_meta("81819NED")
