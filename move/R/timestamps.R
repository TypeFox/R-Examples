setGeneric("timestamps", function(this) standardGeneric("timestamps"))
setMethod("timestamps", ".MoveTrack",
          function(this) {
            this@timestamps
          })


setMethod("timestamps", ".unUsedRecords",
          function(this) {
            this@timestampsUnUsedRecords
          })

setMethod("timestamps", ".MoveTrackSingle",
          function(this) {
            this@timestamps
          })

setGeneric("timestamps<-", function(this, value) standardGeneric("timestamps<-"))
setReplaceMethod("timestamps", ".MoveTrack",
                 function(this, value) {
                   if (length(value)!=length(this@timestamps)) 
                     stop(paste("The number of timestamps does not match the original number of timestamps! (",length(value),":",length(this@timestamps),")"))
                   this@timestamps <- value
                   this
                 })
