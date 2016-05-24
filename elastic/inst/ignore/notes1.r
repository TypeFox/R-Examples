######
# parse_input <- function(...){
#   x <- as.character(dots(...))
#   neg <- gsub('-', '', x[grepl("-", x)])
#   pos <- x[!grepl("-", x, )]
#   list(neg=neg, pos=pos)
# }
# 
# dots <- function(...){
#   eval(substitute(alist(...)))
# }

#####
# select <- function (.data, ...) {
#   select_(.data, .dots = lazyeval::lazy_dots(...))
# }
# 
# select_.data.frame <- function (.data, ..., .dots) {
#     {
#       dots <- lazyeval::all_dots(.dots, ...)
#       vars <- select_vars_(names(.data), dots)
#       select_impl(.data, vars)
#     }
# 
# lazy_dots <- function (..., .follow_symbols = FALSE) {
#   if (nargs() == 0) 
#     return(structure(list(), class = "lazy_dots"))
#   .Call(make_lazy_dots, environment(), .follow_symbols)
# }


body <- '{
 "mappings": {
   "record": {
     "properties": {
         "location" : {"type" : "geo_shape"}
      }
   }
 }
}'
index_create(index='geonames', body=body)

body <- '{
 "query":{
   "geo_shape" : {
     "location" : {
         "shape" : {
           "type": "envelope",
            "coordinates": [[-30, 50],[30, 0]]
         }
       }
     }
   }
}'
out <- Search('geonames', body = body)
out$hits$total

index("geonames") %>% 
   geoshape(field = "location", type = "circle", radius = "5000km", 
            coordinates = c(-10, 45)) %>% 
   n()
