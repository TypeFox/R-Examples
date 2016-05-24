getPara <- function(){
	para<-NULL
	para$datafile<-''
	para$savefile<-''
	para$pngdir<-getwd()
	para$BAFfilter<-10	## filter original file to remove small noisy segments and improve computational efficacy
	para$thr.kmeans<-0.1
	para$thr.originsize<-500
	para$thr.CL<-0.05
	para$thr.CP<-0.05
	para$thr.penalty<-500
	para$std.BAF<-0.01
	para$std.LRR<-0.4
	para$exclude.chr<-NULL
	para$res.r<-0.8 # permutation resample ratio
	para$num.tracks<-4 # number of canonical tracks concerned
	para$model<-1 #1: contraction towards diploid; 2: contraction towards baseline
	para$is.normalize<-TRUE	#if data needs to be normalized. If T, the data will be normalized
	para$is.perm<-FALSE
	para$is.png<-TRUE
    para$is.plot<-TRUE
	para$is.seed<-TRUE
	para$is.bubble<-TRUE
    para$LRR_correction_del<-0.572
	para$LRR_correction_amp<-0.553
	para$save.origin<-TRUE
	para$thr.triploidy<-0.06 #BAF threshold for identifying a tri-ploid sample
	para$is.tri<-FALSE
	return(para)
}