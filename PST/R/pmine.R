## Frequent patterns

setMethod("pmine", signature=c(object="PSTf", data="stslist"), 
	def=function(object, data, l, pmin=0, pmax=1, prefix, lag, 
		average=FALSE, output="sequences", with.prefix=TRUE, sorted=TRUE, decreasing=TRUE, score.norm=FALSE) {

	A <- alphabet(object)

	if (missing(lag)) { lag <- 0 }

	data <- unique(data)

	if (!missing(prefix)) {
		prefix.length <- length(seqdecomp(prefix))
		lag <- prefix.length
		sel <- which(seqconc(data[,1:lag])==prefix)
		prefix.SPS <- suppressMessages(seqformat(prefix, from="STS", to="SPS", compressed=TRUE,SPS.out=list(xfix="", sdsep="/")))

		if (length(sel)>0) {
			message(" [>] selecting ", length(sel)," sequence(s) starting with: ", prefix.SPS)  
			data <- data[sel,]
		} else {
			stop(" [!] no sequence starting with: ", prefix.SPS)
		}
	}

	prob <- predict(object, data, decomp=TRUE)
	prob <- log(prob, base=2)
	pmin <- log(pmin, base=2)
	pmax <- log(pmax, base=2)

	nbps <- rowSums(!is.na(prob))

	if (!average) { 
		prob.check <- prob>=pmin & prob<=pmax
	}
	
	select.seq <- vector(mode="logical", length=nrow(data))
	score <- vector(mode="numeric", length=nrow(data))
	score <- log(score)

	sl <- seqlength(data)

	patterns.list <- NULL

	if (missing(l)) { l <- sl-lag }

	for (p in (1+lag):(max(sl)-max(l)+1)) {
		score.tmp <- rowSums(prob[, p:(p+max(l)-1)])/sl
		if (average) {
			fp <- score.tmp>=pmin & score.tmp<=pmax
		} else {
			fp <- rowSums(prob.check[, p:(p+max(l)-1), drop=FALSE], na.rm=TRUE)==l
		}

		## Tag as selected
		select.seq[fp] <- TRUE

		## Updating score  
		score.update <- which(score.tmp>score) 
		score[score.update] <- score.tmp[score.update] 

		if (output=="patterns" & sum(fp)>0) {
			pstart <- if (with.prefix) { 1 } else { p } 
			tmp <- seqconc(data[fp,  pstart:(p+max(l)-1)])
			patterns.list <- c(patterns.list, unique(tmp[!tmp %in% patterns.list]))
		}
	}

	if (sum(select.seq)>0) {
		if (output=="patterns") {
			message(" [>] ", length(patterns.list), " distinct pattern(s) found")

			nr <- if ("*" %in% A) { "#" } else { "*" }
			res <- seqdef(patterns.list, alphabet=A, labels=object@labels, cpal=cpal(object), nr=nr)
		} else {
			message(" [>] ", sum(select.seq), " sequence(s) selected")

			score <- score[select.seq]
			score <- 2^score
			if (!score.norm) { 
				sl <- sl[select.seq]
				score <- score^sl 
			}

			data <- data[select.seq,]

			if (sorted) {
				score.sort <- order(score, decreasing=decreasing)
				data <- data[score.sort,]
				score <- score[score.sort]
			}

			attr(data, "weights") <- score
			res <- data
		}
	} else {
		message(" [>] no pattern found")
		res <- NULL
	}	

	return(res)
}
)




