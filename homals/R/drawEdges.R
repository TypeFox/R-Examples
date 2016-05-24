`drawEdges` <-
function(a,far=1000) {
  lns<-makeEdges(a)
  p<-dim(lns)[1]; n<-dim(a)[1]
  for (i in 1:p) {
  	ee<-lns[i,3:4]; ff<-lns[i,5:6]
  	xlw<-lns[i,7]; if (xlw == -Inf) xlw<--far
  	plw<-ee+xlw*ff
  	xup<-lns[i,8]; if (xup == Inf) xup<-far
  	pup<-ee+xup*ff
  	lines(rbind(plw,pup),col="RED")
  	}
}

