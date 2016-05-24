ExtractMotion<-function(father,pos=1,
		url="http://localhost:9200/selected_db/selected_db/_search") {
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
                },
                {
                    "text": {
                        "ev_father": "to_be_replaced"
                    }
                }
            ],
            "must_not": [],
            "should": []
        }
    },
    "size": 50,
    "sort": []
  }
  ')
  nm<-names(str$query$bool$must[[2]]$text);
  str$query$bool$must[[2]]$text <- father;
  names(str$query$bool$must[[2]]$text)<-nm;
  res2 <- fromJSON(getURL(url, customrequest="GET",
      httpheader=c('Content-Type'='application/json'),
      postfields=toJSON(str)))
  if (  res2$hits$total == 1 ) {
    kk<-res2$hits$hits[[1]]$`_source`$video_hash$frames
    return(as.numeric(unlist(lapply(kk,function(x){return(x[1])}))));
  } else {
    return(NULL);
  }
}
#
