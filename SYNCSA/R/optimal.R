optimal<-function (comm, envir, traits, subset.min = 2, subset.max = 3, pattern = "tcap", dist = "euclidean", method = "pearson", scale = TRUE, scale.envir = TRUE , na.rm = FALSE, notification = TRUE, progressbar=FALSE) 
{
    part.cor <- function(rxy, rxz, ryz) {
        (rxy - rxz * ryz)/sqrt(1 - rxz * rxz)/sqrt(1 - ryz * ryz)
    }
    comm<-as.matrix(comm)
    envir<-as.matrix(envir)
    traits<-as.matrix(traits)
    if(notification==TRUE){
    	c.NA <- apply(comm, 2, is.na)
    	if(length(which(unique(as.vector(c.NA))==TRUE))>0){
			warning("Warning: NA in community data",call.=FALSE)	
    	}
		t.NA <- apply(traits, 2, is.na)
    	if(length(which(unique(as.vector(t.NA))==TRUE))>0){
			warning("Warning: NA in traits matrix",call.=FALSE)	
		}
		e.NA <- apply(envir, 2, is.na)
		if(length(which(unique(as.vector(e.NA))==TRUE))>0){
			warning("Warning: NA in environmental data",call.=FALSE)	
    	}
    }
    colnames(traits) <- colnames(traits, do.NULL = FALSE, prefix = "T")
    if (scale.envir == "TRUE") {
        envir <- cent.norm(envir,na.rm = na.rm)
    }
    dist.y <- vegdist(envir, method = dist, na.rm = na.rm)
    m <- dim(traits)[2]
    if (subset.max > m) {
        stop("\n Subset must be lower than the number of traits\n")
    }
    PATTERNS <- c("tcap", "tdap", "tcap.tdap")
    pattern <- pmatch(pattern, PATTERNS)
    if (length(pattern) > 1) {
        stop("\n Only one argument is accepted in pattern \n")
    }
    if (is.na(pattern)) {
        stop("\n Invalid pattern \n")
    }
    p <- 1:subset.max
    bin <- factorial(m)/(factorial(p) * factorial(m - p))
    nT<-sum(bin[subset.min:subset.max])
    comb <- matrix(NA, nrow = sum(bin[subset.min:subset.max]), ncol = 1)
    n=0
    for (i in subset.min:subset.max) {
        combinations <- combn(colnames(traits), i, simplify = TRUE)
        for (j in 1:bin[i]) {
        	n=n+1
            comb[n, 1] <- paste(combinations[,j], collapse = " ")
        }
    }
    n=0
    correlation <- matrix(NA, nrow = sum(bin[subset.min:subset.max]), ncol = 1)
    for (i in subset.min:subset.max) {
        combinations1 <- combn(colnames(traits), i, simplify = TRUE)
        for (j in 1:bin[i]) {
            if (pattern == 1) {
            	n=n+1
                T <- matrix.t(comm, as.matrix(traits[, combinations1[, j]]), scale = scale, notification = FALSE)
                correlation[n, 1] <- cor(vegdist(as.matrix(T$matrix.T), method = dist, na.rm = na.rm), dist.y, method = method)
                if(progressbar){
					ProgressBAR(n,nT,style=3)
				}
            }
            if (pattern == 2) {
            	n=n+1
                T <- matrix.t(comm, as.matrix(traits[, combinations1[, j]]), scale = scale, notification = FALSE)
                X <- matrix.x(comm, as.matrix(traits[, combinations1[, j]]), scale = scale, notification = FALSE)
                dist.x <- vegdist(X$matrix.X, method = dist, na.rm = na.rm)
                dist.z <- vegdist(T$matrix.T, method = dist, na.rm = na.rm)
                rxy <- cor(dist.x, dist.y, method = method)
                rxz <- cor(dist.x, dist.z, method = method)
                ryz <- cor(dist.y, dist.z, method = method)
                correlation[n, 1] <- part.cor(rxy, rxz, ryz)
                if(!is.na(rxz)){
                	if(!is.na(ryz)){
                		if(!is.na(rxy)){
                			if ((rxz == 1 | ryz == 1) == TRUE) {
	                  			correlation[n, 1] <- 0
	                  		}
        	        	}
    	            }
                }
				if(progressbar){
					ProgressBAR(n,nT,style=3)
				}
            }
            if (pattern == 3) {
            	n=n+1
                X <- matrix.x(comm, as.matrix(traits[, combinations1[, j]]), scale = scale, notification = FALSE)
                correlation[n, 1] <- cor(vegdist(as.matrix(X$matrix.X), method = dist, na.rm = na.rm), dist.y, method = method)
				if(progressbar){
					ProgressBAR(n,nT,style=3)
				}
            }
        }
    }
    result <- data.frame(Subset = comb, ro = correlation, stringsAsFactors = FALSE)
    return(result[order(result[, 2], decreasing = TRUE), ])
}