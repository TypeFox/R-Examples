NormalizeLRR <- function(x,para){
	peak.int<-c(para$LRR_correction_del,para$LRR_correction_amp)
	print(peak.int)
	dat<-x
	delPos <- which(dat[,4]<0)
	ampPos <- which(dat[,4]>=0)
	dat[delPos,4] <- dat[delPos,4]*(1/para$LRR_correction_del)
	dat[ampPos,4] <- dat[ampPos,4]*(1/para$LRR_correction_amp)
	x <- dat
	return(x)
}