segmentWrapper <-
function(CNAvec,chrom,maploc,ranseed,dtype,...){
	if(all(CNAvec==chrom)|all(CNAvec==maploc))return(NULL)
	if(!is.null(ranseed))set.seed(ranseed)
	CNAobj<-CNA(genomdat=CNAvec,chrom=chrom,maploc=maploc,data.type=dtype)
	seg<-segment(CNAobj,...)
	return(data.frame(seg$output,seg$segRows))
}
