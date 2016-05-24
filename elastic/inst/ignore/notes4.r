library("dplyr")
library("jsonlite")
library("stringr")

foo <- function(x, file){
  other <- as.list(x[ !names(x) %in% c('latitude','longitude','iter')])
  other$location <- list(type="point", coordinates=sprintf(gsub("\\s", "", '[%s,%s]'), x[['longitude']], x[['latitude']]))
  out <- gsub("\\]\"", "\\]", gsub("\"\\[", "\\[", jsonlite::toJSON(other, auto_unbox = TRUE)))
  
  index <- list(index = list(`_index` = "geonames", `_type` = "record", `_id` = str_trim(x[["iter"]])))
  indexout <- gsub("\"}}", "}}", gsub("_id\":\"", "_id\":", jsonlite::toJSON(index, auto_unbox = TRUE)))
  
  cat(indexout, file=file, sep = "\n", append = TRUE)
  cat(out, file=file, sep = "\n", append = TRUE)
}

src <- src_sqlite("/Users/sacmac/Downloads/geonames.db")

options(scipen=999)
starts <- seq(4000000, 6648291, 100000)
indx <- 45:(44+length(starts))
for(i in seq_along(starts)){
  cat(i, sep = "\n")
  query <- sprintf("select * from geotab limit 100000 offset %s", starts[i])
  tab <- tbl(src, sql(query))
  df <- tab %>% 
    data.frame %>%
    tbl_df()
  df$iter <- seq(starts[i], starts[i]+(NROW(df)-1), 1)
  apply(df, 1, foo, file=sprintf("~/github/ropensci/elastic/inst/geonames/geonames%s.json", indx[i]))
}
