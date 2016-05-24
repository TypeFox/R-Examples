 # # # # # # # # # # # # # # # # # #
# # #           feedback        # # #
 # # #          functions       # # #
  # # # # # # # # # # # # # # # # # 

#####  computation of after-before differences

after.minus.before=function(data,operator){
    if(is.matrix(data)){
    	before.after=data
    } else {
    	if(is(data,"KDD") || is(data,"KDD.yearly.average")){
            before.after=data@before.after
    	} else {
    		stop("[after.minus.before] Error: data is neither a matrix nor a KDD / KDD.yearly.average object.")
    	}
    }
    if(nrow(before.after)/2==round(nrow(before.after)/2) || nrow(before.after)<3){
        stop("[after.minus.before] Error: the number of rows of matrix before.after should be at least 3 and should be odd.")
    }
    nb.days=floor(nrow(before.after)/2)
    before=before.after[nb.days:1,]
    after=before.after[nb.days+1+1:nb.days,]
    if(operator=="dmv"){
        difference=apply(after,2,mean)-apply(before,2,mean)
    } else {
        if(operator=="dmpiv"){
            difference=apply(after>0,2,mean)-apply(before>0,2,mean)
        } else {
            if(operator=="dmgiv"){
                difference=apply(after>before,2,mean)-apply(before>after,2,mean)
            } else {
                stop(paste("[after.minus.before] Error: operator -",operator,"- is not available.",sep=""))
                difference=NULL
            }
        }
    }
    return(difference)
}

## feedback statistics

feedback.stats=function(object,operator,turning.year=NULL,
    trend.correction=list(apply=FALSE,object2=NULL)){
    feedback.stats.original=function(object,operator,turning.year){
        diff=after.minus.before(object,operator=operator)
	stats=mean(diff)
	names="whole period"
	if(length(turning.year)>0){
            for(i in turning.year){
                stats=c(stats,mean(diff[object@year<i]),mean(diff[object@year>=i]))
                stats=c(stats,stats[length(stats)]-stats[length(stats)-1])
                names=c(names,paste("period <",i,sep=""),paste("period >=",i,sep=""),"change")
            }
	}
	names(stats)=names
	return(stats)
    }
    out=feedback.stats.original(object,operator,turning.year)
    if(trend.correction$apply){
        alldates=substring(trend.correction$object2@date,first=6)
        STATS=NULL
        for(i in 1:length(object@date)){
            keep=(substring(object@date,first=6)[i]==alldates)
            objecti=kdd(trend.correction$object2@before.after[,keep],
                trend.correction$object2@date[keep],
                trend.correction$object2@year[keep],
                trend.correction$object2@day[keep],
                trend.correction$object2@keyday.threshold)
            STATS=rbind(STATS,feedback.stats.original(objecti,
                operator,turning.year))
        }
        out=out-apply(STATS,2,mean)
    }
    return(out)
}


###################     NONSPATIAL FEEDBACK TESTS     ####################

feedback.test=function(object,test,operator,nb.rand,plots=TRUE,
    turning.year=NULL){
    if(test=="feedback"){
        permute.fct=function(object){
            kddstar=object
            permute=(rbinom(n,1,0.5)==1)
            kddstar@before.after[,permute]=kddstar@before.after[p:1,permute]
            return(kddstar)
        }
    }
    if(test=="change.in.feedback"){
        permute.fct=function(object){
            kddstar=object
            data1=(1:n)[kddstar@year<turning.year]
            permute=c(data1,sample(data1,n-length(data1),replace=TRUE))
            kddstar@before.after=kddstar@before.after[,permute]
            return(kddstar)
        }
    }
    diff=cumsum(after.minus.before(object@before.after,operator=operator))
    DIFF=diff
    n=ncol(object@before.after)
    p=nrow(object@before.after)
    for(j in 1:nb.rand){
        kddstar=permute.fct(object)
        diffstar=cumsum(after.minus.before(kddstar@before.after,
            operator=operator))
        DIFF=rbind(DIFF,diffstar)
    }
    quant=apply(DIFF,2,quantile,c(0.5,0.25,0.75,0.025,0.975))
    nbout=NULL
    for(j in 1:nrow(DIFF)){
        nbout=c(nbout,sum(!(DIFF[j,]>=quant[4,] & DIFF[j,]<=quant[5,])))
    }
    pvalue=mean(nbout>=nbout[1])
    if(plots){
        plot(object@year,diff,ylim=range(c(diff,DIFF)),type="l",col=2,
             xlab="Year",ylab="Cumul(After-Before)",
             cex.axis=1.5,cex.lab=1.5,lwd=1.5)
        lines(object@year,quant[1,],lwd=1.5)
        lines(object@year,quant[2,],lty="dashed",lwd=1.5)
        lines(object@year,quant[3,],lty="dashed",lwd=1.5)
        lines(object@year,quant[4,],lty="dotted",lwd=1.5)
        lines(object@year,quant[5,],lty="dotted",lwd=1.5)
        hist(nbout,breaks=0:max(c(1,nbout)),cex.axis=1.5,cex.lab=1.5,lwd=1.5,
             main=paste("p-value:",round(pvalue,digits=3)),
             xlab="Number of exits")
        abline(v=nbout[1],col=2,lwd=1.5)
    }
    return(pvalue)
}


