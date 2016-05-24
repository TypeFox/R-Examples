setGeneric("move2ade", function(x){standardGeneric("move2ade")})
setMethod("move2ade", 
          signature=".MoveTrackSingle", 
          definition=function(x){ 
            SpatialPointsDataFrame(as(x,'SpatialPoints'),
                                   data=data.frame(id=rep(rownames(idData(x, drop=F)), n.locs(x)))) 
            })

setMethod("move2ade", 
          signature=".MoveTrackStack", 
          definition=function(x){ 
#            SpatialPointsDataFrame(coords=coordinates(x), 
 #                                  data=data.frame(id=as.character(trackId(x)))) 
	  SpatialPointsDataFrame(as(x, 'SpatialPoints'), data=data.frame(id=as.character(trackId(x))))
            })


# define ltrajs when neede
if(!isClass("ltraj"))
    setOldClass("ltraj")


setAs("Move", "ltraj", function(from){
      if(isLonLat(from))
	      warning('Converting a long lat projected object while the ltraj does not deal with long lat projected data')
      adehabitatLT::as.ltraj(as.data.frame(coordinates(from)),date=timestamps(from), id=rownames(from@idData), typeII=T, infolocs=data.frame(sensor=from@sensor,from@data))
})
setAs("MoveStack", "ltraj", function(from){
      if(isLonLat(from))
	      warning('Converting a long lat projected object while the ltraj does not deal with long lat projected data')
      adehabitatLT::as.ltraj(as.data.frame(coordinates(from)),date=timestamps(from), id=from@trackId, typeII=T, infolocs=data.frame(sensor=from@sensor,from@data))
})

setAs("ltraj", "Move", function(from) {
    if (!inherits(from, "ltraj"))
        stop("from should be of class \"ltraj\"")
    if(length(from)!=1)
	    stop("Can only convert one individual to a move object")
    if(!attr(from,"typeII"))
	    stop('Can only work on typeII objects')
    spdf<-adehabitatLT::ltraj2spdf(from)
    new("Move",data=(attr(from[[1]],'infolocs')), spdf, sensor=rep(factor("unknown"), nrow(spdf)), timestamps=spdf$date, idData=data.frame(row.names=paste0(attr(from[[1]], 'id'),'_', attr(from[[1]],'id')),burst=attr(from[[1]],'burst'), id=attr(from[[1]],'id')), sensorUnUsedRecords=factor(levels='unknown'))


})

setAs("ltraj", "MoveStack", function(from) {
    if (!inherits(from, "ltraj"))
        stop("from should be of class \"ltraj\"")
    res<-list()
    for(i in 1:length(from))
    {
	    res[[i]]<-as(from[i], 'Move')
    }
    moveStack(res)


})

