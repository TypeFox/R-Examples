 # # # # # # # # # # # # # # # # # #
# # #     Class  definition     # # #
 # # #          KT.output            # # #
  # # # # # # # # # # # # # # # # # 


.KT.output.validity <- function(object){
	return(TRUE)
}

setClass(
    Class="KT.output",
    slots=c(
        krige.output="list",
        subregion="list",
        averageKrigingPrediction.rand="numeric",
        averageKrigingPrediction.obs="numeric",
        alternative="character",
		p.value="numeric"
   ),
prototype=prototype(
        krige.output=list(),
        subregion=list(),
        averageKrigingPrediction.rand=numeric(),
        averageKrigingPrediction.obs=numeric(),
        alternative=character(),
        p.value=numeric()
   ),
    validity=.KT.output.validity
)


cat("### Show ###\n")
.KT.output.show <- function(object){
    cat("   ~~~ Class:",class(object),"~~~ ")
    cat("\n ~ Observed average kriging prediction: ",object@averageKrigingPrediction.obs)
    cat("\n ~ Test alternative: ",object@alternative)
    cat("\n ~ test p-value: ",object@p.value)
    cat("\n")
    return(invisible())
}
setMethod(f="show",signature="KT.output",definition=.KT.output.show)


cat("### Getteur ###\n")
.KT.output.get <- function(x,i,j,drop){
    switch(EXPR=i,
    "krige.output"={return(x@krige.output)},
    "subregion"={return(x@subregion)},
    "averageKrigingPrediction.rand"={return(x@averageKrigingPrediction.rand)},
    "averageKrigingPrediction.obs"={return(x@averageKrigingPrediction.obs)},
    "alternative"={return(x@alternative)},
    "p.value"={return(x@p.value)},
       stop("[KT.output:get] ",i," is not a 'KT.output' slot")
    )
    return(invisible())
}
setMethod(f="[",signature="KT.output",definition=.KT.output.get)


cat("### Setteur ###\n"	)
.KT.output.set <- function(x,i,j,value){
    switch(EXPR=i,
       "krige.output"={x@krige.outpu=value},
       "subregion"={x@subregion=value},
       "averageKrigingPrediction.rand"={x@averageKrigingPrediction.rand=value},
       "averageKrigingPrediction.obs"={x@averageKrigingPrediction.obs=value},
       "alternative"={x@alternative=value},
       "p.value"={x@p.value=value},
       stop("[KT.output:set] Error:",i," is not a 'KT.output' slot")
    )
    validObject(x)
    return(x)
}
setMethod(f="[<-",signature="KT.output",definition=.KT.output.set)

.KT.output.summary=function(object){
	cat("Krige.test output:\n")
    cat("\n ~ Observed average kriging prediction: ",object@averageKrigingPrediction.obs)
    cat("\n ~ Test alternative: ",object@alternative)
    cat("\n ~ test p-value: ",object@p.value)
    cat("\n")
	return(invisible())
}
setMethod(f="summary",signature="KT.output",definition=.KT.output.summary)

.KT.output.names=function(x){
	slotNames(x)
}
setMethod(f="names",signature="KT.output",definition=.KT.output.names)

.KT.output.plot=function(x, digits=list(predict=5,pvalue=5), breaks=12){
	plot(x@krige.output$MAP,type="l",axes=F,xlab="",ylab="",main="",asp=1)
	polygon(x@subregion$x,x@subregion$y,border=2)
	hist(x@averageKrigingPrediction.rand,
		sub=paste("Observed value =",round(x@averageKrigingPrediction.obs,
		digits=digits$predict)),
		main=paste("p-value =",round(x@p.value,digits=digits$pvalue)),
        xlab="Average interpolated value in the red region",
        xlim=range(c(x@averageKrigingPrediction.obs,
        	x@averageKrigingPrediction.rand)),
        breaks=breaks)
	box()
	abline(v=x@averageKrigingPrediction.obs,col=2)
    return(invisible())
}
setMethod(f="plot",signature="KT.output",definition=.KT.output.plot)





