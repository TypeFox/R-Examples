## ==============================
## Convert from STS to SRS format
## ==============================

STS_to_SRS <- function(seqdata,nrep) {

	nbseq <- seqdim(seqdata)[1]
	seql <- seqdim(seqdata)[2]

	if (is.null(nrep)) nrep <- 0

	out <- data.frame(id=1:nbseq,idx=(1+nrep),t(apply(seqdata,1,seqshift,(1+nrep))))

	names(out) <- c("id","idx",paste("T-",seq(seql-1,1),sep=""),"T")

	for (i in (2+nrep):seql) {
			tmp <- data.frame(id=seq(1:nbseq),idx=i,t(apply(seqdata,1,seqshift,i)))
			names(tmp) <- c("id","idx",paste("T-",seq(seql-1,1),sep=""),"T")
			out <- rbind(out,tmp)
		}

	out <- out[order(out$id,out$idx),]
	names(out) <- c("id","idx",paste("T-",seq(seql-1,1),sep=""),"T")

	return(out)
}

