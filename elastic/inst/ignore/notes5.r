library("elastic")

body <- '{
 "mappings": {
   "record": {
     "properties": {
         "location" : {"type" : "geo_point"}
      }
   }
 }
}'
index_create(index='gbifgeopoint', body=body)
path <- system.file("examples", "gbif_geopoint.json", package = "elastic")
docs_bulk(path)

### Points within a bounding box
body <- '{
 "query":{
    "filtered":{
      "query": {
        "bool": {
          "must_not": {
            "term": {"publishingCountry": "BE"}
          }
        }
      },
      "filter":{
        "geo_bounding_box" : {
          "location" : {
            "top_left" : {
              "lat" : 60,
              "lon" : 40
            },
            "bottom_right" : {
              "lat" : 10,
              "lon" : -20
            }
          }
       }
     }
   }
 }
}'
out <- Search('gbifgeopoint', body = body)
out$hits$total
do.call(rbind, lapply(out$hits$hits, function(x) x$`_source`$location))



library("elastic")
body <- '{
  "query":{
    "filtered":{
      "filter":{
        "geo_distance":"200km",
            "place.geo_point": {
              "lon": 4,
              "lat": 50
            }
        }
    }
  }
}'
connect(es_base = "http://ubicity.ait.ac.at:8080/es", NULL)
out <- Search(index="twitter_all_geo",body=body)
out$hits$total



connect("http://192.168.59.103:9200")
connection() 
ping()
cat_indices()


body <- '{
"query":{
"filtered":{
"query":{
"match" : {"msg.lang":"EN"}
},
"filter":{
"geo_distance": {
"distance":"100km",
"place.geo_point": {
"lon": 30.52193,
"lat": 39.76548
}
}
}
}
}
}'

connect()
connection()
out <- Search(index="twitter_all_geo-2014-11-01", type="ctweet", body=body)
out <- Search(index="twitter_all_geo", type="ctweet", body=body)
out$hits$total

connect(es_base = "https://fishbase.ropensci.org",
        es_port = 9200,
        es_user = "elasticsearch", 
        es_pwd = "c7587757ba32511b365e43fdb6ea6905a6db8af2f4854b313daa081a5d0e6fe5")
connection()


cat_()
cat_count()
cat_indices()
nodes_info()
index_exists("logstash-2015.04.02")
Search_uri("logstash-2015.04.02")
Search_uri("logstash-2015.04.02", asdf = TRUE)
entries <- pluck(Search_uri("logstash-2015.04.02")$hits$hits, c("_source", "message"))
lapply(entries, jsonlite::fromJSON)

library("httr")
alias_get(config=c(verbose(), ssl.verifypeer = FALSE))
alias_get(ssl.verifypeer = FALSE)

foo(verbose())

foo <- function(...) {
  userpwd <- make_up()
  tt <- GET("http://127.0.0.1:9200/_alias?ignore_unavailable=true", c(userpwd, ...))
  if(tt$status_code > 202) geterror(tt)
  jsonlite::fromJSON(content(tt, as = "text"), FALSE)
}
