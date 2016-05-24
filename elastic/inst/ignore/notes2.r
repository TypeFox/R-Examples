countUAlangEN <- '{
  "query": {
    "bool": {
      "must": [
         {"term": {"ctweet.place.country_code": "UA"}},
         {"term": {"ctweet.msg.lang":"EN"}}
      ],
      "must_not": [ ],
      "should": [ ]
    }
  },
  "sort": [ ],
  "facets": { }
}'
...
Search(index=nextElem(twitter.index), type="ctweet", search_type="count", body=countUAlangEN)$hits$total

Search(index="shakespeare")
Search(index="shakespeare", search_type = "query_then_fetch")
Search(index="shakespeare", search_type = "dfs_query_then_fetch")
Search(index="shakespeare", search_type = "count")
Search(index="shakespeare", search_type = "scan", scroll = "5m")
Search(index="shakespeare", search_type = "query_and_fetch")
Search(index="shakespeare", search_type = "dfs_query_and_fetch")

#For whatever reason, the country code and language must be lower case in this query (ua,en)
countUAlangEN <- '{
  "query": {
    "bool" : {
      "must" : [
        {"term": {"ctweet.place.country_code":"ua"}},
        {"term": {"ctweet.msg.lang":"en"}}
      ]
    }
  },
  "aggs": {
    "by_date": {
      "date_histogram": {
        "field": "created_at",
        "interval": "day",
        "format": "yyyy-MM-dd"
      }
    }
  }
}'
...
Search(index="twitter_all_geo", type="ctweet", body=countUAlangEN)
...
