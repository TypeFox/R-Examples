#
VideoSearch<-function(url="http://localhost:9200/selected_db/selected_db/_search") {
  str <- fromJSON ('
    {
    "facets": {},
    "from": 0,
    "query": {
        "bool": {
            "must": [
                {
                    "term": {
                        "FileName": "_video.json"
                    }
                }
            ],
            "must_not": [],
            "should": []
        }
    },
    "size": 500,
    "sort": []
  }
  ')
  res2 <- fromJSON(getURL(url, customrequest="GET",
      httpheader=c('Content-Type'='application/json'),
      postfields=toJSON(str)));
  vnam<-NULL
  if (  res2$hits$total > 0 ) {
    vnam<- sort(unlist(lapply(res2$hits$hits,
                function(x){return(x$`_source`$ev_father)})))
  }
  return(vnam)
}
#
