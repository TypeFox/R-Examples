getUnifiedMap <- function(seg.L.chr,seg.B.chr,map.chr,thr=5){
	#Merge breakpoints within thr by seg.L
	seg.new<-c()
	i<-1
	j<-1
	while(1){
		pos.B<-which(map.chr==seg.B.chr[j])
		pos.L<-which(map.chr==seg.L.chr[i])
		if(seg.L.chr[i]>=seg.B.chr[j]){
			if(i%%2==j%%2 & pos.L-pos.B<=thr){
				seg.new<-c(seg.new,seg.L.chr[i])
				i<-i+1
				j<-j+1
			}
			else{
				seg.new<-c(seg.new,seg.B.chr[j])
				j<-j+1
			}
		}
		else{
			if(i%%2==j%%2 & pos.B-pos.L<=thr){
				seg.new<-c(seg.new,seg.L.chr[i])
				i<-i+1
				j<-j+1
			}
			else{
				i<-i+1
			}
		}
		if(i>length(seg.L.chr)|j>length(seg.B.chr)){break}
	}
	return(seg.new)
}