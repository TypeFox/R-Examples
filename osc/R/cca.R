cca <- function(data, s=1, mode=3, count.cells=FALSE, count.max=ncol(data)*3, res.x=NULL, res.y=NULL, cell.class=1, unit="", compare=""){
	if(class(data) %in% c("data.frame", "RasterLayer")){
    if(unit=="m"){
      ccaRm(data=data, d=s,res.x=res.x, res.y=res.y, cell.class=cell.class, compare)
    } else {
      ccaRd(data=data, d=s, cell.class=cell.class, compare)
    }
	} else {
		ccaM(data=data, s=s,mode=mode, count.cells=count.cells, count.max=count.max)
	}

}
ccaM <- function(data, s=1, mode=3, count.cells=FALSE, count.max=ncol(data)*3){
	#do checks
	stopifnot(is.numeric(data))
	stopifnot(is.matrix(data))
	stopifnot(is.numeric(s))
	stopifnot(mode==1 | mode==2 | mode==3)
	the.data <- as.integer(t(data))
	clu <- as.integer(rep(0, ncol(data)*nrow(data)))
        count.max <- as.integer(count.max)
	count <- as.integer(rep(0, count.max))
	if(count.cells==TRUE){
	  out <- .C("callburn_count",  s=as.integer(s), xmax=nrow(data), ymax=ncol(data), mode=as.integer(mode)[1], data=the.data, clu=clu, count=count, count.max=count.max,CLASSES=c("integer", "integer", "integer", "integer", "integer","integer"))
	  return(list(clusters=matrix(out$clu, ncol=ncol(data), byrow=TRUE), cluster.count=out$count[1:max(out$clu)]))
	  }
	if(count.cells==FALSE){
	  out <- .C("callburn",  s=as.integer(s), xmax=nrow(data), ymax=ncol(data), mode=as.integer(mode)[1], data=the.data, clu=clu, count=count, count.max=count.max,CLASSES=c("integer", "integer", "integer", "integer", "integer","integer"))
	  return(clusters=matrix(out$clu, ncol=ncol(data), byrow=TRUE))
	  }
}

