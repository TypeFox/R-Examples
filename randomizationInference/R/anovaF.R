#calculate ANOVA F statistic for 1 factor, >2 treatments
#can be used for CRD, BRD, and Latin squares
#calcOptions=list(block=)
#calcOptions=list(row=,col=)

anovaF=function(y,w,calcOptions=NULL){
	#formatting
	if(sum(factor(w)==w)>0) w=data.matrix(w)
	if(ncol(w)==1) w=as.vector(w)
	#calculations
	#default to complete randomization
	if(is.null(calcOptions)){
		anova(lm(y~factor(w)))$F[1]
	}else if(is.null(calcOptions$block)==FALSE){
		anova(lm(y~factor(w)+factor(calcOptions$block)))$F[1]
	}else if(is.null(calcOptions$row)==FALSE){
		anova(lm(y~factor(w)+factor(calcOptions$row)+factor(calcOptions$col)))$F[1]
	}
}
