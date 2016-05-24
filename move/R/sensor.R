setGeneric("sensor", function(this,...) standardGeneric("sensor"))
setMethod("sensor", ".MoveTrack",
          function(this,...) {
            this@sensor
          })
setMethod("sensor", ".unUsedRecords",
          function(this,...) {
            this@sensorUnUsedRecords
          })
