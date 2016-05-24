# Extract frames where images were taken out 
ExtractImgPos<-function(father,
				url="http://localhost:9200/selected_db/selected_db/_search") {
  str <- fromJSON ('
    {
    "facets": {},
    "from": 0,
    "query": {
        "bool": {
            "must": [
                {
                    "text": {
                        "ev_father": "to_be_replaced"
                    }
                }
            ],
            "must_not": [
                {
                    "term": {
                        "FileName": "_video.json"
                    }
                }
            ],
            "should": []
        }
    },
    "size": 50,
    "sort": []
  }
  ');
  nm<-names(str$query$bool$must[[1]]$text);
  str$query$bool$must[[1]]$text <- father; 
  names(str$query$bool$must[[1]]$text)<-nm;
  res2 <- fromJSON(getURL(url, customrequest="GET",
      httpheader=c('Content-Type'='application/json'),
      postfields=toJSON(str)));
  if (  res2$hits$total > 0 ) {
    out<-unlist(lapply(res2$hits$hits,function(x){return(x$`_source`$FileName)}));
    return(unique(sort(unlist(lapply(strsplit(out,"[-_]"),function(x){return(x[3])})))));
  } else {
    return(NULL);
  }
}
#
