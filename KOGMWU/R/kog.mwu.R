kog.mwu <-
function(data,gene2kog,Alternative="t") {
	rsq=data
	names(rsq)=c("seq","value")
	bads=which(rsq[,2]==Inf | rsq[,2]==(-Inf) | is.na(rsq[,2]))
	if (length(bads)>0) { rsq=rsq[-bads,]}
	kogs=gene2kog
	annotated=rsq[rsq[,1] %in% kogs[,1],]
	kogs=kogs[kogs[,1] %in% rsq[,1],]
	
	kogrows=match(kogs[,1],annotated[,1])
	annotated=annotated[kogrows,]
	annotated$term=as.character(kogs[,2])
	annotated$value=as.numeric(as.character(annotated[,2]))
	
	mwut.t=TRUE
	if (length(levels(as.factor(annotated[,2])))==2) {
		print("Binary classification detected; will perform Fisher's test");
		mwut.t=F
		rr=kog.ft(annotated)
	} else {
		print("Continuous measure of interest: will perform MWU test");		
		rr=kog.mwut(annotated,Alternative)
	}
	return(rr)
}
