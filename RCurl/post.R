library(RCurl)
u = "http://requestb.in/168944w1"

json = '{
      "project": 1,
          "event_id": "fc6d8c0c43fc4630ad850ee518f1b9d0",
          "culprit": "my.module.function_name",
          "timestamp": "2011-05-02T17:41:36",
          "message": "SyntaxError: Wattttt!"
          "sentry.interfaces.Exception": {
                    "type": "SyntaxError":
                              "value": "Wattttt!",
                              "module": "__builtins__"
                  }
    }'

hdr = c('X-Sentry-Auth' = 'Sentry sentry_version=2.0, sentry_signature=a3901c854752a61636560638937237c8d7a9561d, sentry_timestamp=1329096377, sentry_key=b70a31b3510c4cf793964a185cfe1fd0, sentry_client=raven-python/1.0')


httpPOST(u, httpheader = hdr, postfields = json, verbose = TRUE)
