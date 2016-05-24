`pava` <-
function(x,w=rep(1,length(x)),block=weighted.mean){
  nblock<-n<-length(x); blocklist<-array(1:n,c(n,2)); blockvalues<-x; active<-1
  repeat{
  	if (!isUpSatisfied(blockvalues,active)) {
  		blockmerge<-mergeBlockup(blocklist,blockvalues,x,w,active,block)
  		blockvalues<-blockmerge$v; blocklist<-blockmerge$l
  		nblock<-nblock-1
  		while (!isDownSatisfied(blockvalues,active)) {
  			blockmerge<-mergeBlockup(blocklist,blockvalues,x,w,active-1,block)
  			blockvalues<-blockmerge$v; blocklist<-blockmerge$l;
  			nblock<-nblock-1; active<-active-1;
  			}
  		}
  	else if (active == nblock) break() else active<-active+1
  	}
  putBack(n,blocklist,blockvalues)
}

