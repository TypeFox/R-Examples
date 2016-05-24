"batons" <- 
function(..., waist = FALSE) {
	res <- boxplot(..., boxwex=0, medpch=19, medcex=1.5, boxlwd=4, whisklty=1, outpch=4, outcex=.5)
	if(waist){
		segments(seq(ncol(res$conf)), res$conf[1,], seq(ncol(res$conf)), res$conf[2,], col="white", lwd=1.2)
	}
}