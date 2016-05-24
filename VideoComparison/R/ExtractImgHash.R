#
# Extract hashes from one specific frame where images were taken out 
ExtractImgHash<-function(pos,father,
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
                },
                {
                    "query_string": {
                        "default_field": "FileName",
                        "query": "to_be_replaced"
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
  nm<-names(str$query$bool$must[[2]]$query_string[2]);
  str$query$bool$must[[2]]$query_string[2]<-paste("*",pos,"_out*",sep="");
  names(str$query$bool$must[[2]]$query_string[2])<-nm;
  res2 <- fromJSON(getURL(url, customrequest="GET",
      httpheader=c('Content-Type'='application/json'),
      postfields=toJSON(str)));
  if (  res2$hits$total == 1 ) {
    out<-res2$hits$hits[[1]]$`_source`$Hash;
    rt<-list(dct=as.character(out[1]), 
             hstrada=as.numeric(unlist(strsplit(out[2],","))),
             mw=as.numeric(unlist(strsplit(out[3],","))),
             rd=as.numeric(unlist(strsplit(out[4],",")))
             )
    
    return(rt);
  } else {
    return(NULL);
  }
}
