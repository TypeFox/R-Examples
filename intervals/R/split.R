# Cause split.data.frame to be used. See "Methods for S3 Generic Functions" in
# help for "Methods", for guidance and examples.

split.Intervals_virtual <- split.data.frame

setMethod( "split", "Intervals_virtual", split.Intervals_virtual )
