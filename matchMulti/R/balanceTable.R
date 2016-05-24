balanceTable <-
function(df.orig, df.match = NULL, treatment, #cat.vars = NULL,
 treat.wts = NULL, ctrl.wts = NULL, mt.wts = NULL, mc.wts = NULL, verbose = FALSE){
	#if(is.null(cat.vars)) cat.vars <- rep(FALSE, ncol(df.orig))
	
	#Validate input
	if (!(treatment %in% colnames(df.orig))) {
		stop(paste('Treatment variable', treatment,'not found'))
	}	
	if(!is.null(df.match) & !all(colnames(df.orig) == colnames(df.match))){ 
		stop('df.orig and df.match must have identical column names')
	}
	#if(!is.null(cat.vars) && length(cat.vars) != ncol(df.orig)){
	#	stop('cat.vars must have exactly one entry for each column in df.orig')
	#}
	if(any(is.na(df.orig[[treatment]])) || (!is.null(df.match) && any(is.na(df.match[[treatment]])))){
		stop('NAs are present in the treatment variable')
	}
	
	#non.numeric <- colnames(df.orig)[laply(df.orig, inherits, what = c('character','factor'))]
	
	cov.orig <- handleNA(df.orig, verbose = verbose)
	if(!is.null(df.match)){		
		cov.match <- handleNA(df.match, verbose = verbose)
		cov.match <- resolve.cols(cov.match, cov.orig)
		cov.orig <- resolve.cols(cov.orig, cov.match)
	} else {
		cov.match <- NULL
	}
	
	#search through and find all binaries.
	binary.ind <- laply(cov.orig, is.binary)
	treat.ind <- colnames(cov.orig) == treatment
	
	sdiff.out <- aaply(colnames(cov.orig)[which(!treat.ind)], 1, sdiff, treatment = treatment, orig.data = cov.orig, match.data = cov.match, treat.wts = treat.wts, ctrl.wts = ctrl.wts, mt.wts = mt.wts, mc.wts = mc.wts)
	rownames(sdiff.out) <- colnames(cov.orig)[which(!treat.ind)]

	t.test.out <- as.matrix(aaply(colnames(cov.orig)[which(!treat.ind)], 1, ttest.balance, treatment = treatment, orig.data = cov.orig, match.data = cov.match, treat.wts = treat.wts, ctrl.wts = ctrl.wts, mt.wts = mt.wts, mc.wts = mc.wts)) 
	
	#TODO: figure out whether to keep weight arguments and incorporate them or drop them
	
	if(any(!binary.ind)){
		wilc.test.out <- as.matrix(aaply(colnames(cov.orig)[which(!binary.ind & !treat.ind)], 1, wilc.balance, treatment = treatment, orig.data = cov.orig, match.data = cov.match, treat.wts = treat.wts, ctrl.wts = ctrl.wts, mt.wts = mt.wts, mc.wts = mc.wts))
		if (sum(!binary.ind & !treat.ind) == 1) wilc.test.out <- t(wilc.test.out)
		rownames(wilc.test.out) <- colnames(cov.orig)[which(!binary.ind & !treat.ind)]
	}

	if(any(binary.ind)){
		fisher.test.out <- as.matrix(aaply(colnames(cov.orig)[which(binary.ind & !treat.ind)], 1, fisher.balance, treatment = treatment, orig.data = cov.orig, match.data = cov.match, treat.wts = treat.wts, ctrl.wts = ctrl.wts, mt.wts = mt.wts, mc.wts = mc.wts)) 
		rownames(fisher.test.out) <- colnames(cov.orig)[which(binary.ind & !treat.ind)]
	}
		
	if (all(binary.ind)){
		test.out <- fisher.test.out
	} else if (all(!binary.ind[-which(treat.ind)])) {
		test.out <- wilc.test.out	
	} else {
		test.out <- as.data.frame(matrix(nrow = nrow(sdiff.out), ncol = ncol(fisher.test.out)))
		rownames(test.out) <- rownames(sdiff.out)
		if (ncol(test.out) == 2){
			colnames(test.out) <- c('Fisher/Wilcox Pvalue Before', 'Fisher/Wilcox Pvalue After')
		} else {
			colnames(test.out) <- 'Fisher/Wilcox Pvalue'
		}
		test.out[which(rownames(sdiff.out) %in% rownames(fisher.test.out)),] <- fisher.test.out
		test.out[which(rownames(sdiff.out) %in% rownames(wilc.test.out)),] <- wilc.test.out
	}
	
	out.tab <- cbind(sdiff.out,t.test.out, test.out)
	out.tab
}
