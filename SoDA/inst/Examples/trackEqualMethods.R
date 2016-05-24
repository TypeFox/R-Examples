setMethod("==", c("track", "track"),
          function(e1, e2) {
              e1@x == e2@x &
                e1@y == e2@y
          })

setMethod("!=", c("track", "track"),
          function(e1, e2) {
              e1@x != e2@x |
                e1@y != e2@y
          })
