setGeneric("burst", function(x, f, ...) {standardGeneric("burst")})
setMethod("burst", 
          signature=c(x = "Move", f = "factor"), 
          definition = function(x, f, ...) {
            levels(f)<-validNames(levels(f))
          new("MoveBurst", 
              as(x, ".MoveTrackSingle"), 
              as(x, ".MoveGeneral"), 
              burstId = f, 
              idData = x@idData)
          })

setMethod("burst", 
          signature=c(x = "Move", f = "numeric"), 
          definition = function(x, f, ...) {
            burst(x = x, f = as.factor(f))
          })

setMethod("burst", 
          signature=c(x = "Move", f = "character"), 
          definition = function(x, f, ...) {
            if (length(f) == 1) {
                burst(x = x, f = do.call("$", list(x, f))[-n.locs(x)])
            } else {
                burst(x = x, f = as.factor(f))
            }
          })

setMethod('burst',
	  signature=c(x='.MoveTrackSingleBurst', 'missing'),
	  definition=function(x,f,...){
		  return(slot(x,'burstId'))
	  })


setAs("MoveBurst", "Move", function(from) {
# last id can be different since that one is not telling anything about a segment in this obj
      if (length(unique(from@burstId)) != 1) 
        stop("Does not work with one burst id only")
      new("Move", 
          as(from,".MoveGeneral"), 
          as(from,".MoveTrackSingle"))
    }) 
setAs("MoveStack", "Move", function(from) {
# last id can be different since that one is not telling anything about a segment in this obj
      if (length(unique(from@trackId)) != 1) 
        stop("Does only work with one id")
      new("Move", 
          as(from,".MoveGeneral"), 
          as(from,".MoveTrack"),timestamps=timestamps(from), as(unUsedRecords(from),'.unUsedRecords')
	  )
    }) 
setAs("dBMvarianceBurst", "dBMvariance", function(from) {
      if (length(unique(from@burstId)) != 1) 
	      stop("Not one unique burst id method wont work")
      new("dBMvariance", as(from, "dBMvarianceTmp"), as(from, ".MoveTrackSingle"))
	 }) 
