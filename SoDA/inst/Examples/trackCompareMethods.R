setMethod("Compare", c("track", "track"),
          function(e1, e2) {
              cmpx <- callGeneric(e1@x, e2@x)
              cmpy <- callGeneric(e1@y, e2@y)
              ifelse(cmpx & cmpy, TRUE,
                     ifelse(cmpx | cmpy, NA, FALSE))
          })
