## Impute missign values in sequence data using a PST

setMethod("impute", signature=c(object="PSTf", data="stslist"), 
	def=function(object, data, method="pmax") {
	
	ipdata <- as.data.frame(data)
	L <- summary(object)@depth
		
	A <- alphabet(data)
	if (!all.equal(alphabet(object), alphabet(data))) {
		stop(" [!] object and data do not have same alphabet")
	}

	nr <- attr(data, "nr")
	sl <- ncol(data)

	ismiss <- matrix(nrow=nrow(ipdata), ncol=ncol(ipdata))
	hasmiss <- vector("logical", length=nrow(data))

	for (i in 1:sl) {
		ismiss[,i] <- ipdata[,i]==nr
	}

	nbmiss <- rowSums(ismiss)
	hasmiss <- nbmiss>0

	message(" [>] found ", sum(hasmiss), " sequence(s) with missing values")
	message(" [>] ", sum(nbmiss), " missing states ")

	for (i in which(hasmiss)) {
		for (l in 1:sl) {
			if (ipdata[i,l]==nr) {
				if (l==1) {
					context <- "e"
				} else {
					context <- seqconc(ipdata[i, max(1, l-L):(l-1)])
				}

				p <- suppressMessages(query(object, context))

				## Imputation
				if (method=="pmax") {
					ipdata[i,l] <- A[which.max(p)]
				} else {
					ipdata[i,l] <- sample(A, size=1, prob=p)
				}
			}
		}
	}

	ipdata <- seqdef(ipdata, alphabet=A, missing=nr, labels=stlab(data), cpal=cpal(data), 
		weights=attr(data, "weights"))

	return(ipdata)		
}
)
		

	
