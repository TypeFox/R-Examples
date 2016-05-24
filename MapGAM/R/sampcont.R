sampcont <-
function(rdata, type="stratified", regions=NULL, times=NULL, n=1, nrow=100, ncol=100) { 
  if (type=="stratified") {
    if (!is.null(times) && length(times)!=dim(rdata)[1]) 
	  stop("times argument must be NULL or a vector/factor of length = rows in rdata")
    if (is.null(regions)) {
      XYnames = names(rdata)[2:3]
      Xrange = range(rdata[,2])
      Yrange = range(rdata[,3])
      polyGrid = PBSmapping::makeGrid(x=seq(Xrange[1],Xrange[2],length.out=nrow),
      		y=seq(Yrange[1],Yrange[2],length.out=ncol))
      names(rdata)[2:3] = c("X","Y")
      rdata$EID = 1:length(rdata$X)
      idpolys = PBSmapping::findCells(PBSmapping::as.EventData(rdata),polyGrid)
      rdata = merge(idpolys,rdata)
      regions = paste((1:nrow)[rdata$PID],(1:ncol)[rdata$SID],sep=",")
      rdata = rdata[,!(names(rdata) %in% c("EID","PID","SID","Bdry"))]
      names(rdata)[2:3] = XYnames
    } else
    if (length(regions)!=dim(rdata)[1]) 
	  stop("regions argument must be NULL or a vector/factor of length = rows in rdata")
    strata = paste(regions,times,sep="")
    eligible = rdata[,1]==0		  	# eligible controls
    counts = table(strata[eligible])		# number of eligible controls in each stratum
    estrat = names(counts)			# strata with at least 1 eligible control
    ns = pmin(n,counts)	    			# number to sample from each stratum
    wts = counts/ns				# inverse probability weights
    dsamp = rdata[rdata[,1]==1,]		# take all cases  	
    ncases = dim(dsamp)[1]			# number of cases
    w = rep(1,ncases)				# cases get weight=1
    # sample controls from each stratum and add to dsamp data frame
    for (i in 1:length(counts)) {
      ind = sample(1:counts[i],ns[i])		# sample without replacement
      controls = rdata[strata==estrat[i] & eligible,][ind,]
      dsamp = rbind(dsamp,controls)
      w = c(w,rep(wts[i],dim(controls)[1]))	# keep track of sampling weights
    }
	n = length(w)-ncases		# redefine n as total numer of controls
    cat(paste(n,"controls selected from",sum(eligible),
	"eligibles in",length(estrat),"strata."),fill=T)
  }
  if (type=="simple") {
    eligible = rdata[,1]==0		  		# eligible controls
    if (n > sum(eligible)) stop(paste("rdata contains only ",n," eligible controls"))
    dsamp = rdata[rdata[,1]==1,]			# take all cases  	
    ind = sample(1:dim(rdata)[1],n,prob=eligible)	# simple random sample
	dsamp = rbind(dsamp,rdata[ind,])
	w = c(rep(1,sum(!eligible)),rep(n/sum(eligible),n))
  }
  return(list(rdata=dsamp,w=w,ncont=n))
}
