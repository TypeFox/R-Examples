###create a list of Move objects from a Move Stack (hand over additional arguments!)
setGeneric("split") 
setMethod(f = "split",
	  signature = c(x="MoveStack", f="missing"),
	  definition = function(x, f, ...){
		 s<-split(as(x,'.MoveTrackStack'))
	  return(lapply(s, new, Class='Move', as(x, '.MoveGeneral')))
	  })
setMethod(f = "split",
	  signature = c(x=".MoveTrackStack", f="missing"),
	  definition = function(x, f, ...){
		  moveList <- list()
		  unUsed<-as(x,".unUsedRecordsStack")
		  for (ID in unique(x@trackId)) {
			  s<-x@trackId==ID
			  spdf <- SpatialPointsDataFrame(coords = matrix(x@coords[s,], ncol=2),
							 data=x@data[s,],
							 proj4string=x@proj4string)
			  mt <- new(Class=".MoveTrack",
				    spdf,
				    timestamps=x@timestamps[s],
				    sensor=x@sensor[s])
			  unUsedSub<-as(unUsed[unUsed@trackIdUnUsedRecords==ID,],'.unUsedRecords')
			  moveObj <- new(Class=".MoveTrackSingle", 
					 mt,
					 idData=x@idData[row.names(x@idData)==ID, ,drop=F],
					 unUsedSub)
			  moveList[[ID]]  <- moveObj
		  }
		  return(moveList)
	  }) 

setMethod(f = "split",
	  signature = c(x="DBBMMStack", f="missing"),
	  definition = function(x, f, ...){
		  DBBMMList <- list()
		  moveTrackStack<-split(as(x@DBMvar,'.MoveTrackStack'))
		  for (Id in as.character(unique(x@DBMvar@trackId))) { 
			  UD <- new(Class=".UD", 
				    method=x@method,
				    x[[Id]])
			  dbmv<-new('dBMvariance',moveTrackStack[[Id]],
				      window.size=x@DBMvar@window.size,
				      break.list=x@DBMvar@break.list,
				      interest=as.logical(x@DBMvar@interest[x@DBMvar@trackId==Id]),
				      means=x@DBMvar@means[x@DBMvar@trackId==Id],
				      in.windows=x@DBMvar@in.windows[x@DBMvar@trackId==Id],
				      margin=x@DBMvar@margin)
			  DBBMMObj <- new(Class="DBBMM",
					  UD,
					  ext=x@ext,
					  DBMvar=dbmv)
			  DBBMMList[[Id]]  <- DBBMMObj
		  }
		  return(DBBMMList)
	  })

setMethod(f = "split",
	  signature = c(x=".UDStack", f="missing"),
	  definition = function(x, f, ...){
		  xx<-lapply(unstack(x), new, Class=".UD")
		  names(xx)<-names(x)
		  return(xx)
	  })

setMethod("split", 
          signature = c(x = ".MoveTrackSingleBurst", f = "missing"), 
          definition = function(x, f, ...) {    
            f <- c(0, cumsum(diff(as.numeric(factor(paste(as.character(x@burstId), is.na(x@burstId)))))!=0))
            f <- c(f, max(f))
            res <- list()
            for (i in unique(f)) 
              res[[i + 1]] <- x[f == i | c(0, f[-n.locs(x)]) == i, ]
            names(res) <- as.character(unlist(lapply(lapply(res, slot, "burstId"), "[", 1)))
            res <- lapply(res, as, sub("Burst", "", class(x)))
            return(res)
          })


