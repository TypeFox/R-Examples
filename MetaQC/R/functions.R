GetPVal <- function(dat, resp.type=c("Twoclass", "Multiclass", "Survival")) {
	
	t.test2 <- function(x,y) {
		varr <- function(x, meanx=NULL){
			n <- ncol(x)
			p <- nrow(x)
			Y <-matrix(1,nrow=n,ncol=1)
			if(is.null(meanx)){   meanx <- rowMeans(x)}
			ans<- rep(1, p)
			xdif <- x - meanx %*% t(Y)
			ans <- (xdif^2) %*% rep(1/(n - 1), n)
			ans <- drop(ans)
			return(ans)
		}
		n1 <- sum(y==levels(y)[1])
		n2 <- sum(y==levels(y)[2])
		
		m1 <- rowMeans(x[,y==levels(y)[1],drop=F])
		m2 <- rowMeans(x[,y==levels(y)[2],drop=F])
		
		df <- n1+n2-2
		sd <- sqrt( ((n2-1) * varr(x[, y==levels(y)[2]], meanx=m2) + (n1-1) * varr(x[, y==levels(y)[1]], meanx=m1) )*(1/n1+1/n2)/df )
		
		tstat <-  (m2 - m1) / sd
		pval <- 2 * pt(-abs(tstat), df)
		return(pval)
	}
	
	if(resp.type == "Twoclass") {
		group <- factor(dat$y)
		if(any(is.na(dat$x))) {
			p.value <- apply(dat$x, 1, function(x) {
						t.test(x[which(group==levels(group)[1])], x[which(group==levels(group)[2])], var.equal=T)$p.value 
					}) 
		} else {
			p.value <- t.test2(dat$x, group)
		}
	} else if(resp.type == "Multiclass") {
		p.value <- apply(dat$x, 1, function(xx) summary(lm(dat$y ~ xx))$coefficients[2,4])
	} else if(resp.type == "Survival") {
		suppressPackageStartupMessages(require(survival))
		p.value <- apply(dat$x, 1, function(xx) summary(coxph(Surv(dat$y, dat$censoring.status) ~ xx))$logtest['pvalue'])
	}
	
	return(p.value)
}

GetEWPval <- function(PDat) {
	pchisq(-2*rowSums(log(PDat),na.rm=TRUE), df=2*rowSums(!is.na(PDat)), lower.tail=F)
}

printLog <- function(msg, verbose) {
	if(verbose)
		print(sprintf("[%s] %s", Sys.time(), msg))
}

#From gtools
combinations <- function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) {
	if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) != 
			0) 
		stop("bad value of n")
	if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) != 
			0) 
		stop("bad value of r")
	if (!is.atomic(v) || length(v) < n) 
		stop("v is either non-atomic or too short")
	if ((r > n) & repeats.allowed == FALSE) 
		stop("r > n and repeats.allowed=FALSE")
	if (set) {
		v <- unique(sort(v))
		if (length(v) < n) 
			stop("too few different elements")
	}
	v0 <- vector(mode(v), 0)
	if (repeats.allowed) 
		sub <- function(n, r, v) {
			if (r == 0) 
				v0
			else if (r == 1) 
				matrix(v, n, 1)
			else if (n == 1) 
				matrix(v, 1, r)
			else rbind(cbind(v[1], Recall(n, r - 1, v)), Recall(n - 
										1, r, v[-1]))
		}
	else sub <- function(n, r, v) {
			if (r == 0) 
				v0
			else if (r == 1) 
				matrix(v, n, 1)
			else if (r == n) 
				matrix(v, 1, n)
			else rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), 
						Recall(n - 1, r, v[-1]))
		}
	sub(n, r, v[1:n])
}

#From matrixStats
rowVars <- function (x, center = NULL, ...) {
	n <- !is.na(x)
	n <- rowSums(n)
	n[n <= 1] <- NA
	if (is.null(center)) {
		center <- rowMeans(x, ...)
	}
	x <- x - center
	x <- x * x
	x <- rowSums(x, ...)
	x <- x/(n - 1)
	x
}

GMT2List <- function(f, saveAs=NULL) {
	trim <- function(x) {
		sub(" +$", "", sub("^ +", "", x))
	}
	
	con <- file(f, "r", blocking = FALSE)
	GList <- readLines(con) # empty
	GList <- lapply(GList, function(x) unlist(strsplit(x, "\t", fixed=T)))
	GList <- foreach(x=iter(GList)) %dopar% sapply(strsplit(x, "///", fixed=T), function(xx) trim(xx[[1]]))
	GList.names <- unlist(foreach(x=iter(GList)) %dopar% x[1])
	GList <- foreach(x=iter(GList)) %dopar% x[-(1:2)]
	names(GList) <- GList.names
	close(con)
	if(!is.null(saveAs))
		save(GList, file=saveAs)
	return(GList)
}

getFileExt <- function(x) {
	sub(".*[.]", "", x)
}

getFileName <- function(x) {
	sub("(.+)[.][^.]+$", "\\1", basename(x))
}

union.rec <- function(.list, ...){
	if(length(.list)==1) return(.list[[1]])
	Recall(c(list(union(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
}

Download <- function(pkg, fn) {
	url <- paste("http://cloud.github.com/downloads/donkang75", pkg, fn, sep="/")
	res <- try(download.file(url, fn, mode = "wb"), silent=TRUE)
	return(res)
}
