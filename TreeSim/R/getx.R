getx<-function(datatree,sersampling=0){
    datatree$node.label<-NULL	
	if (sersampling==1){
		#9.7.14
		if (length(datatree$tip.label)==2) {
			temp<-datatree$edge.length
			ttype <- c(0,0,1)
			temp2<-0
			if (length(temp)==3) {ttype<-c(ttype,1)
				temp2<-temp[1]
				temp<-temp[-1]
				}
			times<-ttype*0
			times[2]<-abs(temp[1]-temp[2])
			times[3]<-max(temp[1],temp[2])
			if (temp2>0) {times[4]<-temp2+times[3]}
			x<-cbind(times,ttype)
		} else {
		x<-get.times.polytomy(datatree)}} else {
	br<-branching.times(datatree)
	edges<-datatree$edge
	edges1<-edges[,1]
	edges2<-edges[,2]
	ord<-order(edges1)
	edges2<-cbind(edges1[ord],edges2[ord])
	br2<-vector()
	temp<-edges[1,1]
	index<-1
	for (j in 2:length(edges2[,1])){
		if (edges2[j,1] == temp){
		index <- index+1
		label<-paste(edges2[j,1],sep="")
		br2<-c(br2,br[label])
		} else {
		if (index==1){
		print("wierd")
		break
		}
		index<-1
		temp<-edges2[j,1]
		}
	}
	x<-sort(br2)}
	x
}