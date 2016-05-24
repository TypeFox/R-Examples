UCexactTest = function(object, pair = 1:2, dispersion = "auto", prior.count = 0.125, upper=5e4, by=100)
{
    if (!is(object, "DGEList")) 
        stop("Currently only supports DGEList objects as the object argument.")
    if (length(pair) != 2) 
        stop("Pair must be of length 2.")
    group <- as.factor(object$samples$group)
    levs.group <- levels(group)
    if (is.numeric(pair)) 
        pair <- levs.group[pair]
    else pair <- as.character(pair)
    if (!all(pair %in% levs.group)) 
        stop("At least one element of given pair is not a group.\n Groups are: ", 
            paste(levs.group, collapse = " "))
    if (is.null(dispersion)) 
        dispersion <- "auto"
    if (is.character(dispersion)) {
        dispersion <- match.arg(dispersion, c("auto", "common", 
            "trended", "tagwise"))
        dispersion <- switch(dispersion, common = object$common.dispersion, 
            trended = object$trended.dispersion, tagwise = object$tagwise.dispersion, 
            auto = getDispersion(object))
        if (is.null(dispersion)) 
            stop("specified dispersion not found in object")
        if (is.na(dispersion[1])) 
            stop("dispersion is NA")
    }
    ldisp <- length(dispersion)
    ntags <- nrow(object$counts)
    if (ldisp != 1 && ldisp != ntags) 
        stop("Dispersion provided by user must have length either 1 or the number of tags in the DGEList object.")
    if (ldisp == 1) 
        dispersion <- rep(dispersion, ntags)
    group <- as.character(group)
    j <- group %in% pair
    y <- object$counts[, j, drop = FALSE]
    			lib.size <- object$samples$lib.size[j]
    norm.factors <- object$samples$norm.factors[j]
    group <- group[j]
    if (is.null(rownames(y))) 
        rownames(y) <- paste("tag", 1:ntags, sep = ".")
    			lib.size <- lib.size * norm.factors
    			offset <- log(lib.size)
    			lib.size.average <- exp(mean(offset))
    			prior.count <- prior.count * lib.size/mean(lib.size)
    			offset.aug <- log(lib.size + 2 * prior.count)
    j1 <- group == pair[1]
    n1 <- sum(j1)
    if (n1 == 0) 
        stop("No libraries for", pair[1])
    y1 <- y[, j1, drop = FALSE]
    abundance1 <- mglmOneGroup(y1 + matrix(prior.count[j1], ntags, 
        n1, byrow = TRUE), offset = offset.aug[j1], dispersion = dispersion)
    j2 <- group == pair[2]
    n2 <- sum(j2)
    if (n1 == 0) 
        stop("No libraries for", pair[2])
    y2 <- y[, j2, drop = FALSE]
    abundance2 <- mglmOneGroup(y2 + matrix(prior.count[j2], ntags, 
        n2, byrow = TRUE), offset = offset.aug[j2], dispersion = dispersion)
    logFC <- (abundance2 - abundance1)/log(2)
    abundance <- mglmOneGroup(y, dispersion = dispersion, offset = offset)
    e <- exp(abundance)
    input.mean <- matrix(e, ntags, n1)
    output.mean <- input.mean * lib.size.average
    input.mean <- t(t(input.mean) * lib.size[j1])
    y1 <- q2qnbinom(y1, input.mean = input.mean, output.mean = output.mean, 
        dispersion = dispersion)
    input.mean <- matrix(e, ntags, n2)
    output.mean <- input.mean * lib.size.average
    input.mean <- t(t(input.mean) * lib.size[j2])
    y2 <- q2qnbinom(y2, input.mean = input.mean, output.mean = output.mean, 
        dispersion = dispersion)
      
      y1[y1<0] = 0 #after quantile adjustment, sometimes values are neg.
	y2[y2<0] = 0
	n1 = ncol(y1)
	n2 = ncol(y2)
	if(n1==n2) {
		s1 = rowSums(y1) 
		s2 = rowSums(y2)
		exact.pvals <- pvalue(s1=s1,s2=s2,phi=as.numeric(dispersion),n1=n1,n2=n2,upper=upper,by=by)
	}
	if(n1!=n2) {
		n1.bigger = n1 > n2
		if(n1.bigger) {s1 = rowSums(y1);s2 = rowSums(y2)/n2 * n1}
		if(!n1.bigger) {s1 = rowSums(y1)/n1 * n2;s2 = rowSums(y2)}
		exact.pvals <- pvalue(s1=s1,s2=s2,phi=dispersion,n1=max(n1,n2),n2=max(n1,n2),upper=upper,by=by)
	}
	
    AveLogCPM <- object$AveLogCPM
    if (is.null(AveLogCPM)) 
        AveLogCPM <- aveLogCPM(object)
    de.out <- data.frame(logFC = logFC, logCPM = AveLogCPM, PValue = exact.pvals)
    rn <- rownames(object$counts)
    if (!is.null(rn)) 
        rownames(de.out) <- make.unique(rn)
    new("DGEExact", list(table = de.out, comparison = pair, genes = object$genes))
}