getPara.sAGP <- function(){
	para<-NULL
	para$inputdata<-''	#input data, for segment purity calculation
	para$savedata<-''
	para$pngdir<-''
	para$purityfile<-''
	para$thr.kmeans<-0.1
	para$K<-12	## constant parameter for objective function F (see details in getSegPurity)
	para$BAFfilter<-10
	para$thr.originsize<-500
	para$is.plot<-FALSE
    para$is.normalize<-TRUE
	para$std.BAF<-0.025
	para$std.LRR<-0.1
	para$max.cn<-NULL	## maximum copy number examined. Default NULL, this number is determined by ploidy+4
	para$strictness<-1.1	## if 1, strict. larger than 1 indicating segments can have >1 AGP, due to noise
	para$is.seed<-TRUE
	para$res.r<-0.8 # permutation resample ratio
	para$is.LRRcorrection<-FALSE
	para$is.multicore<-FALSE
	para$is.scale<-FALSE	#When computing minimun distance from segment to canonical track, if rescale BAF, recommended F.
	para$LRR_correction_del<-0.572
	para$LRR_correction_amp<-0.553
	para$originFile<-'Origin.Rdata'
	para$sensitivity<-FALSE	## if perform sensitivity analysis
	para$ss.BAF<-0.02	## radius of random sampling region
	para$ss.LRR<-0.08
	para$ss.N<-10
	para$ss.dir<-''
	return(para)
}