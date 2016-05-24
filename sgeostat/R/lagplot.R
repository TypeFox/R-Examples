# lagplot.s creates a lagplot for an object of type "point"
assign("lagplot",
function(point.obj,pair.obj,a1,a2,lag=1,std=FALSE,query.a=NULL,xlim=NULL,ylim=NULL) {

  if (!inherits(point.obj,"point")) stop('Point.obj must be of class, "point".\n')

  if (!inherits(pair.obj,"pair")) stop('Pair.obj must be of class, "pair".\n')

  if(missing(a1)) stop('Must enter at least one attribute.\n')
  if(missing(a2)) a2 <- a1

  att1 <- point.obj[[match(a1,names(point.obj))]]

  att2 <- point.obj[[match(a2,names(point.obj))]]

  if(std) {
    att1 <- (att1 - mean(att1,na.rm=TRUE))/var(att1[!is.na(att1)])
    att2 <- (att2 - mean(att2,na.rm=TRUE))/var(att2[!is.na(att2)])
  }

  plot((att1[pair.obj$from])[pair.obj$lags==lag],
       (att2[pair.obj$to])[pair.obj$lags==lag],
#       xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),
       xlab=paste(a1,'(s)',sep=''),
#       ylab=paste(a2,'(s+h)',sep=''),fty='s')
       ylab=paste(a2,'(s+h)',sep=''))
  title(paste(deparse(substitute(point.obj)),": lag " ,lag,sep=''))
  abline(0,1)

  if(!is.null(query.a)) {
    i <- (1:length(names(point.obj)))[names(point.obj) == query.a]
    query.att <- point.obj[[i]]
    cat('Identify "from" points...')
    identify((att1[pair.obj$from])[pair.obj$lags==lag],
             (att2[pair.obj$to])[pair.obj$lags==lag],
             (query.att[pair.obj$from])[pair.obj$lags==lag],col=2)
    cat('\nIdentify "to" points...')
    identify((att1[pair.obj$from])[pair.obj$lags==lag],
             (att2[pair.obj$to])[pair.obj$lags==lag],
             (query.att[pair.obj$to])[pair.obj$lags==lag],col=3)

  }
       
})
#

