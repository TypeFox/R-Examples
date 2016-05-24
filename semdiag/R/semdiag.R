## SEM Diagnostics based on EQS
semdiag.combinations<-
function (n, r) 
{
    v<-1:n
    v0 <- vector(mode(v), 0)
    sub <- function(n, r, v) {
            if (r == 0) 
                v0
            else if (r == 1) 
                matrix(v, n, 1)
            else if (n == 1) 
                matrix(v, 1, r)
            else rbind(cbind(v[1], Recall(n, r - 1, v)), Recall(n - 1, r, v[-1]))
        }
    sub(n, r, v[1:n])
}

## Modified functions from REQS
semdiag.read.eqs<-function (file) 
{
    file.ets <- file
    file.split <- strsplit(file.ets, "\\.")
    if (length(file.split[[1]]) > 2) 
        stop("File name or folders should not contain '.'")
    if (file.split[[1]][2] != "ets") 
        stop("File should be of the form 'xxxxxx.ets'")
    file.cbk <- paste(file.split[[1]][1], ".CBK", sep = "")
    file.etp <- paste(file.split[[1]][1], ".ETP", sep = "")
    cbk.info1 <- scan(file.cbk, skip = 2, nlines = 2, quiet = TRUE)
    cbk.info2 <- scan(file.cbk, skip = 4, nlines = 2, quiet = TRUE)
    endfile <- scan(file.cbk, nlines = 1, quiet = TRUE)
    cbk.info.mat <- cbind(cbk.info2, cbk.info1)
    rownames(cbk.info.mat) <- c("Parameter estimates", "Standard errors", 
        "Robust standard errors", "Corrected standard errors", 
        "Gradients", "Sample covariance matrix", "Model Covariance Matrix (Sigma hat)", 
        "Inverted Information matrix", "Robust inverted information matrix", 
        "Corrected inverted information matrix", "First derivatives", 
        "4th Moment weight matrix", "Standardized Elements", 
        "R-squares", "Factor means", "Univariate statistics (means)", 
        "Univariate statistics (standard deviations)", "Univariate statistics (skewness)", 
        "Univariate statistics (kurtosis)", "Univariate statistics (sample size)", 
        "Dependent variable standardization vector", "Independent variable standardization vector")
    colnames(cbk.info.mat) <- c("Line Number", "Number of Elements")
    cbk.base <- read.fwf(file.cbk, widths = c(13, 3), skip = 6, 
        col.names = c("variable", "line"), buffersize = 1, n = 98)
    nminfo <- length(which(cbk.base[, 2] == 2))
    minfo.val <- scan(file.ets, skip = 1, nlines = 1, quiet = TRUE)
    minfo.dframe <- data.frame(minfo.val)
    colnames(minfo.dframe) <- "values"
    rownames(minfo.dframe) <- cbk.base[1:nminfo, 1]
    start.cbk <- nminfo + 1
    ntprobs <- length(c(which(cbk.base[, 2] == 3), which(cbk.base[, 
        2] == 4)))
    probs.val <- scan(file.ets, skip = 2, nlines = 2, quiet = TRUE)
    probs.val[which(probs.val == -1)] <- NA
    probs.dframe <- data.frame(probs.val, row.names = cbk.base[start.cbk:(start.cbk + 
        ntprobs - 1), 1])
    colnames(probs.dframe) <- "p-values"
    start.cbk <- start.cbk + ntprobs
    nfit <- 60
    fit.val <- scan(file.ets, skip = 4, nlines = 6, quiet = TRUE)
    fit.val[which(fit.val == -9)] <- NA
    fit.dframe <- data.frame(fit.val, row.names = cbk.base[start.cbk:(start.cbk + 
        nfit - 1), 1])
    colnames(fit.dframe) <- "fit values"
    start.cbk <- start.cbk + nfit
    ndesc <- 9
    desc.val <- scan(file.ets, skip = 11 - 1, nlines = 1, quiet = TRUE)
    desc.dframe <- data.frame(desc.val, row.names = cbk.base[start.cbk:(start.cbk + 
        ndesc - 1), 1])
    colnames(desc.dframe) <- "values"
    n.ind <- desc.dframe[8, 1]
    n.dep <- desc.dframe[9, 1]
    n.fac <- desc.dframe[3, 1]
    n.tot <- n.ind + n.dep
    if (n.ind%%32 == 0) {
        skiplines <- n.ind/32
    }
    else {
        skiplines <- trunc(n.ind/32) + 1
    }
    if (n.dep%%32 == 0) {
        skiplines <- skiplines + n.dep/32
    }
    else {
        skiplines <- skiplines + trunc(n.dep/32) + 1
    }
    parindvec <- scan(file.etp, skip = skiplines, quiet = TRUE)
    varnames.string <- readLines(file.etp, n = skiplines, warn = FALSE)
    varnames.chvec <- unlist(strsplit(varnames.string, split = " "))
    varnames.vec <- varnames.chvec[which(varnames.chvec != "")]
    nout <- dim(cbk.info.mat)[1]
    model.list <- as.list(rep(NA, nout))
    for (i in 1:nout) {
        startline <- cbk.info.mat[i, 1]
        if (i != nout) {
            endlinevec <- cbk.info.mat[(i + 1):nout, 1]
            endline <- (endlinevec[endlinevec > 0])[1]
            nlines <- endline - startline
        }
        else {
            if ((cbk.info.mat[i, 2]) > 0) 
                nlines <- endfile - startline
            else nlines <- 0
        }
        if (startline != 0) {
            vals <- scan(file.ets, skip = startline - 1, nlines = nlines, 
                quiet = TRUE)
        }
        else {
            vals <- NA
        }
        model.list[[i]] <- vals
    }
    par.val <- model.list[[1]]
    par.pos <- which(parindvec > 0)
    phi.dim <- n.ind * n.ind
    gamma.dim <- n.dep * n.ind
    indpos1 <- (par.pos > (phi.dim)) + (par.pos < (phi.dim + 
        gamma.dim))
    cumpos1 <- par.pos[indpos1 == 2]
    cumpos2 <- par.pos[par.pos > (phi.dim + gamma.dim)]
    if (length(cumpos1) > 0) 
        parindvec[cumpos1] <- parindvec[cumpos1] + max(parindvec)
    if (length(cumpos2) > 0) 
        parindvec[cumpos2] <- parindvec[cumpos2] + max(parindvec)
    negpos <- which(parindvec == -1)
    parindvec[parindvec <= 0] <- NA
    parvec <- par.val[parindvec]
    parvec[negpos] <- -1
    parvec[is.na(parvec)] <- 0
    cuts <- c(n.ind * n.ind, n.dep * n.ind, n.dep * n.dep)
    dimlist <- list(c(n.ind, n.ind), c(n.dep, n.ind), c(n.dep, 
        n.dep))
    cutfac <- rep(1:3, cuts)
    parlist <- split(parvec, cutfac)
    parmat <- mapply(function(xx, dd) {
        matrix(xx, nrow = dd[1], ncol = dd[2], byrow = TRUE)
    }, parlist, dimlist)
    names(parmat) <- c("Phi", "Gamma", "Beta")
    colnames(parmat$Phi) <- rownames(parmat$Phi) <- colnames(parmat$Gamma) <- varnames.vec[1:n.ind]
    rownames(parmat$Gamma) <- rownames(parmat$Beta) <- colnames(parmat$Beta) <- varnames.vec[(n.ind + 
        1):length(varnames.vec)]
    parse.mat <- NULL
    for (i in 1:5) parse.mat <- cbind(parse.mat, model.list[[i]])
    colnames(parse.mat) <- c("Parameter", "SE", "RSE", "CSE", 
        "Gradient")
    npar <- dim(parse.mat)[1]
    namesvec <- NULL
    partable <- NULL
    for (i in 1:3) {
        if (i == 1) {
            combmat <- semdiag.combinations(dim(parmat[[i]])[1], 
                2)
            comb.names <- apply(combmat, 2, function(rn) rownames(parmat[[i]])[rn])
            par.val0 <- parmat[[i]][lower.tri(parmat[[i]], diag = TRUE)]
            partable <- rbind(partable, cbind(comb.names, par.val0))
        }
        if (i == 2) {
            comb.names <- as.matrix(expand.grid(colnames(parmat[[i]]),rownames(parmat[[i]]) 
                ))
            par.val0 <- as.vector(t(parmat[[i]]))
            partable <- rbind(partable, cbind(comb.names, par.val0))
        }
        if (i == 3) {
            comb.names <- as.matrix(expand.grid(rownames(parmat[[i]]), 
                colnames(parmat[[i]])))
            par.val0 <- as.vector(t(parmat[[i]]))
            partable <- rbind(partable, cbind(comb.names, par.val0))
        }
        par.val.ind <- which(((par.val0 != 0) + (par.val0 != 
            -1)) == 2)
        names.mat <- rbind(comb.names[par.val.ind, ])
        
        names <- apply(names.mat, 1, function(ss) paste("(", 
                ss[2], ",", ss[1], ")", sep = ""))       
        namesvec <- c(namesvec, names)
    }
    if ((dim(parse.mat)[1]) != (length(namesvec))) {
        parse.mat <- parse.mat[parse.mat[, 1] != 0, ]
        if ((dim(parse.mat)[1]) == (length(namesvec))) 
            rownames(parse.mat) <- namesvec
    }else {
        rownames(parse.mat) <- namesvec
    }
    meanjn <- scan(file.cbk, skip = 1, nlines = 1, quiet = TRUE)[3]
    Vcheckstr <- colnames(parmat$Phi)
    compstr <- paste("V", 1:999, sep = "")
    TFVcheck <- Vcheckstr %in% compstr
    if (any(TFVcheck)) 
        depnames.add <- Vcheckstr[TFVcheck]
    else depnames.add <- NULL
    VBcheckstr <- colnames(parmat$Beta)
    TFVBcheck <- VBcheckstr %in% compstr
    if (any(TFVBcheck)) 
        depnames.addB <- VBcheckstr[TFVBcheck]
    else depnames.addB <- NULL
    depnames <- c(depnames.addB, depnames.add)
    rm(compstr)
    if (meanjn == 0) 
        p <- n.dep
    else p <- n.dep + 1
    cov.list <- as.list(rep(NA, 5))
    names(cov.list) <- c("sample.cov", "sigma.hat", "inv.infmat", 
        "rinv.infmat", "cinv.infmat")
    for (i in 6:10) {
        if (length(model.list[[i]]) > 1) {
            cov.list[[i - 5]] <- matrix(model.list[[i]], nrow = sqrt(length(model.list[[i]])))
            if (i <= 7) {
                dimnames(cov.list[[i - 5]]) <- list(depnames[1:dim(cov.list[[i - 
                  5]])[1]], depnames[1:dim(cov.list[[i - 5]])[2]])
                order.V <- order(depnames)
                cov.list[[i - 5]] <- cov.list[[i - 5]][order.V, 
                  order.V]
            }
            if (i >= 8) 
                dimnames(cov.list[[i - 5]]) <- list(namesvec, 
                  namesvec)
        }
    }
    pstar <- p * (p + 1)/2
    if (length(model.list[[11]]) > 1) {
        deriv1 <- matrix(model.list[[11]], nrow = npar, ncol = pstar)
    }
    else {
        deriv1 <- NA
    }
    if (length(model.list[[12]]) > 1) {
        moment4 <- matrix(model.list[[12]], nrow = pstar, ncol = pstar)
    }
    else {
        moment4 <- NA
    }
    ustatmat <- cbind(model.list[[16]], model.list[[17]], model.list[[18]], 
        model.list[[19]], model.list[[20]])
    if (dim(ustatmat)[1] == 1) 
        ustatmat <- NA
    else colnames(ustatmat) <- c("means", "sd", "skewness", "kurtosis", 
        "n")
    result <- c(list(model.info = minfo.dframe), list(pval = probs.dframe), 
        list(fit.indices = fit.dframe), list(model.desc = desc.dframe), 
        parmat, list(par.table = parse.mat), cov.list, list(derivatives = deriv1), 
        list(moment4 = moment4), list(ssolution = model.list[[13]]), 
        list(Rsquared = model.list[[14]]), list(fac.means = model.list[[15]]), 
        list(var.desc = ustatmat), list(depstd = model.list[[21]]), 
        list(indstd = model.list[[22]]))
    result
}


semdiag.run.eqs<-function (EQSpgm, EQSmodel, serial, Rmatrix = NA, datname = NA, 
    LEN = 2e+06) 
{
    res <- semdiag.call.eqs(EQSpgm = EQSpgm, EQSmodel = EQSmodel, serial = serial, 
        Rmatrix = Rmatrix, datname = datname, LEN = LEN)
    if (!res) 
        warning("EQS estimation not successful!")
    filedir.split <- strsplit(EQSmodel, "/")[[1]]
    n <- length(filedir.split)
    etsname <- strsplit(filedir.split[n], "\\.")[[1]][1]
    etsfile <- paste(etsname, ".ets", sep = "")
    reslist <- semdiag.read.eqs(etsfile)
    return(c(list(success = res), reslist))
}

semdiag.call.eqs<-function (EQSpgm, EQSmodel, serial, Rmatrix = NA, datname = NA, 
    LEN = 2e+06) 
{
    if (!file.exists(EQSmodel)) 
        stop("The .eqs file not found in the current folder!")
    filedir.split <- strsplit(EQSmodel, "/")[[1]]
    n <- length(filedir.split)
    filedir <- paste(filedir.split[1:(n - 1)], collapse = "/")
    if (n > 1) 
        setwd(filedir)
    outname <- strsplit(filedir.split[n], "\\.")[[1]][1]
    file.out <- paste(outname, ".out", sep = "")
    lenstring <- paste("LEN=", as.integer(LEN), sep = "")
    filepathin <- paste("IN=", EQSmodel, sep = "")
    fileout <- paste("OUT=", file.out, sep = "")
    serstring <- paste("SER=", serial, "\n", sep = "")
    if (length(Rmatrix) > 1) {
        if (is.na(datname)) {
            warning(paste("No filename for data specified! ", 
                outname, ".dat is used", sep = ""))
            datname <- paste(outname, ".dat", sep = "")
        }
        write.table(as.matrix(Rmatrix), file = datname, col.names = FALSE, 
            row.names = FALSE)
    }
    EQScmd <- paste(deparse(EQSpgm), filepathin, fileout, lenstring, 
        serstring)
    RetCode <- system(EQScmd, intern = FALSE, ignore.stderr = TRUE, 
        wait = TRUE, input = NULL)
    if (RetCode == 0) {
        success <- TRUE
    }
    else {
        success <- FALSE
    }
    return(success = success)
}

semdiag.input.model<-function(file = "")
{
    eqs <- scan(file = file, what = "", sep = "\n", strip.white = TRUE, 
        comment.char = "#", fill = TRUE)
    eqs
}



## Parse the EQS input file based on model specification
semdiag.parse<-function(eqs){
	## check whether this is a vector
	if (length(eqs)>1){
		eqsinput<-eqs
		keywordlno<-grep('/', eqsinput, ignore.case=TRUE)
		d1<-NULL
		d2<-NULL
		d3<-NULL
		d4<-NULL
		d5<-NULL
		d6<-NULL
		count<-0
		## parse equations for fixed parameters
		eqslno<-grep('/equations', eqsinput, ignore.case=TRUE)
		eqslend<-keywordlno[which(keywordlno==eqslno)+1]
		for (i in (eqslno+1):(eqslend-1)){
			temp<-gsub(' +', '', eqsinput[i])
			temp<-gsub(';', '', temp)
			temp<-gsub('-', '+-', temp)
			if (nchar(temp)>1){
				
				
				temp<-unlist(strsplit(temp, '='))
				left<-temp[1]
				if (substr(temp[2], 1, 1)=='+') temp[2]<-substring(temp[2], 2)
				right<-unlist(strsplit(temp[2], '+', fixed=TRUE))
				m<-length(right)
				for (j in 1:m){
					if (length(grep('*',right[j], fixed=TRUE))==0){
						str.ind<-unlist(strsplit(right[j],''))
						pos<-ifelse('F' %in% str.ind, which(str.ind=='F'), 0)
						pos<-pos + ifelse('V' %in% str.ind, which(str.ind=='V'), 0)
						pos<-pos + ifelse('E' %in% str.ind, which(str.ind=='E'), 0)
						pos<-pos + ifelse('D' %in% str.ind, which(str.ind=='D'), 0)
						d1<-c(d1, left)
					
						l<-length(str.ind)
						d2<-c(d2, substr(right[j], pos, l))
						if (pos==1) {d3<-c(d3, '1')}else{
						d3<-c(d3, substr(right[j], 1, (pos-1)))}
						d4<-c(d4, 0)
						d5<-c(d5, 1)
						count<-count+1
					}else{
						str.ind<-unlist(strsplit(right[j],'*', fixed=TRUE))
						temp.value<-'0'
						if (nchar(str.ind[1])>0){ temp.value = str.ind[1] }
						d1<-c(d1, left)

						d2<-c(d2, str.ind[2])
						
						d3<-c(d3, temp.value)
						d4<-c(d4, 1)
						d5<-c(d5, 1)
						count<-count+1	
						}
					}
				d6<-c(d6, count)
				}			
			} ## finish parse equations
		## parse variances for fixed parameters
		varlno<-grep('/variances', eqsinput, ignore.case=TRUE)
		varend<-length(eqsinput)+1
		if (length(keywordlno)>2) varend<-keywordlno[which(keywordlno==varlno)+1]
		for (i in (varlno+1):(varend-1)){
			temp<-gsub(' +', '', eqsinput[i])
			temp<-gsub(';', '', temp)

			if (nchar(temp)>1){
				temp<-unlist(strsplit(temp, '='))
				left<-temp[1]
			
				if (length(grep('*',temp[2], fixed=TRUE))==0){
					d1<-c(d1, left)
					d2<-c(d2, left)
					d3<-c(d3, temp[2])
					d4<-c(d4, 0)
					d5<-c(d5, 2)
				}else{
					str.ind<-unlist(strsplit(temp[2],'*', fixed=TRUE))
					d1<-c(d1, left)
					d2<-c(d2, left)
					temp.value<-'0'
					if (nchar(str.ind)>1) temp.value<-str.ind
					d3<-c(d3, temp.value)
					d4<-c(d4, 1)
					d5<-c(d5, 2)
					}	
				}		
			} ## finish parse variances
			
		## parse covariances for fixed parameters
		covlno<-grep('/covariances', eqsinput, ignore.case=TRUE)
		if (length(covlno)>0){
			covend<-length(eqsinput)+1
			if (length(keywordlno)>3) covend<-keywordlno[which(keywordlno==covlno)+1]
			for (i in (covlno+1):(covend-1)){
				temp<-gsub(' +', '', eqsinput[i])
				temp<-gsub(';', '', temp)

				if (nchar(temp)>1){
					temp<-unlist(strsplit(temp, '='))
					left<-temp[1]
			
					if (length(grep('*',temp[2], fixed=TRUE))==0){
						tempvar<-unlist(strsplit(temp[1], ','))
						d1<-c(d1, tempvar[1])
						d2<-c(d2, tempvar[2])
						d3<-c(d3, temp[2])
						d4<-c(d4, 0)
						d5<-c(d5, 3)
					}else{
						tempvar<-unlist(strsplit(temp[1], ','))
						d1<-c(d1, tempvar[1])
						d2<-c(d2, tempvar[2])
						str.ind<-unlist(strsplit(temp[2],'*', fixed=TRUE))
					
						temp.value<-'0'
						if (nchar(str.ind)>1) temp.value<-str.ind
						d3<-c(d3, temp.value)
						d4<-c(d4, 1)
						d5<-c(d5, 3)
						}	
				}
				}		
			}else{ covend<-varend-1} ## finish parse covariance
		
		row.names<-paste('(',d1,',',d2,')',sep='')
		ram<-cbind(d1,d2,d3,d4,d5)
		rownames(ram)<-row.names
		colnames(ram)<-c('end','start','value','free','type')
		eqsstart<-NA
		if (covend==length(eqsinput)) {eqsend=NA} else{ eqsend<-eqsinput[covend:length(eqsinput)]}

		}else{

	eqsinput<-readLines(eqs,warn = FALSE)	
	keywordlno<-grep('/', eqsinput, ignore.case=TRUE)
	method<-grep('/model', eqsinput, ignore.case=TRUE)
	
	d1<-NULL
	d2<-NULL
	d3<-NULL
	d4<-NULL
	d5<-NULL
	d6<-NULL
	count<-0
	
	if (length(method)==0){
		## parse equations for fixed parameters
		eqslno<-grep('/equations', eqsinput, ignore.case=TRUE)
		eqsend<-keywordlno[which(keywordlno==eqslno)+1]
		for (i in (eqslno+1):(eqsend-1)){
			temp<-gsub(' +', '', eqsinput[i])
			temp<-gsub(';', '', temp)
			temp<-gsub('-', '+-', temp)
			if (nchar(temp)>1){
				
				
				temp<-unlist(strsplit(temp, '='))
				left<-temp[1]
				if (substr(temp[2], 1, 1)=='+') temp[2]<-substring(temp[2], 2)
				right<-unlist(strsplit(temp[2], '+', fixed=TRUE))
				m<-length(right)
				for (j in 1:m){
					if (length(grep('*',right[j], fixed=TRUE))==0){
						str.ind<-unlist(strsplit(right[j],''))
						pos<-ifelse('F' %in% str.ind, which(str.ind=='F'), 0)
						pos<-pos + ifelse('V' %in% str.ind, which(str.ind=='V'), 0)
						pos<-pos + ifelse('E' %in% str.ind, which(str.ind=='E'), 0)
						pos<-pos + ifelse('D' %in% str.ind, which(str.ind=='D'), 0)
						d1<-c(d1, left)
					
						l<-length(str.ind)
						d2<-c(d2, substr(right[j], pos, l))
						if (pos==1) {d3<-c(d3, '1')}else{
						d3<-c(d3, substr(right[j], 1, (pos-1)))}
						d4<-c(d4, 0)
						d5<-c(d5, 1)
						count<-count+1
					}else{
						str.ind<-unlist(strsplit(right[j],'*', fixed=TRUE))
						temp.value<-'0'
						if (nchar(str.ind[1])>0){ temp.value = str.ind[1] }
						d1<-c(d1, left)

						d2<-c(d2, str.ind[2])
						
						d3<-c(d3, temp.value)
						d4<-c(d4, 1)
						d5<-c(d5, 1)
						count<-count+1	
						}
					}
				d6<-c(d6, count)
				}			
			} ## finish parse equations
		## parse variances for fixed parameters
		varlno<-grep('/variances', eqsinput, ignore.case=TRUE)
		varend<-keywordlno[which(keywordlno==varlno)+1]
		for (i in (varlno+1):(varend-1)){
			temp<-gsub(' +', '', eqsinput[i])
			temp<-gsub(';', '', temp)

			if (nchar(temp)>1){
				temp<-unlist(strsplit(temp, '='))
				left<-temp[1]
			
				if (length(grep('*',temp[2], fixed=TRUE))==0){
					d1<-c(d1, left)
					d2<-c(d2, left)
					d3<-c(d3, temp[2])
					d4<-c(d4, 0)
					d5<-c(d5, 2)
				}else{
					str.ind<-unlist(strsplit(temp[2],'*', fixed=TRUE))
					d1<-c(d1, left)
					d2<-c(d2, left)
					temp.value<-'0'
					if (nchar(str.ind)>1) temp.value<-str.ind
					d3<-c(d3, temp.value)
					d4<-c(d4, 1)
					d5<-c(d5, 2)
					}	
				}		
			} ## finish parse variances
			
		## parse covariances for fixed parameters
		covlno<-grep('/covariances', eqsinput, ignore.case=TRUE)
		if (length(covlno)>0){
		covend<-keywordlno[which(keywordlno==covlno)+1]
		for (i in (covlno+1):(covend-1)){
			temp<-gsub(' +', '', eqsinput[i])
			temp<-gsub(';', '', temp)

			if (nchar(temp)>1){
				temp<-unlist(strsplit(temp, '='))
				left<-temp[1]
			
				if (length(grep('*',temp[2], fixed=TRUE))==0){
					tempvar<-unlist(strsplit(temp[1], ','))
					d1<-c(d1, tempvar[1])
					d2<-c(d2, tempvar[2])
					d3<-c(d3, temp[2])
					d4<-c(d4, 0)
					d5<-c(d5, 3)
				}else{
					tempvar<-unlist(strsplit(temp[1], ','))
					d1<-c(d1, tempvar[1])
					d2<-c(d2, tempvar[2])
					str.ind<-unlist(strsplit(temp[2],'*', fixed=TRUE))
					
					temp.value<-'0'
					if (nchar(str.ind)>1) temp.value<-str.ind
					d3<-c(d3, temp.value)
					d4<-c(d4, 1)
					d5<-c(d5, 3)
					}	
				}
				}		
			}else{ covend<-varend} ## finish parse covariance
		}
		row.names<-paste('(',d1,',',d2,')',sep='')
		ram<-cbind(d1,d2,d3,d4,d5)
		rownames(ram)<-row.names
		colnames(ram)<-c('end','start','value','free','type')
		eqsstart<-eqsinput[1:(eqslno-1)]
		eqsend<-eqsinput[covend:length(eqsinput)]
		
		## change the output file name
		eqsname <- strsplit(eqs, "\\.")[[1]][1]
		lno<-grep(eqsname, eqsend, ignore.case=TRUE)
     temp.str<-eqsend[lno]
     temp.str.new<-sub(eqsname, 'eqsonce', temp.str)
     eqsend[lno]<-temp.str.new
     
     ## change the iteration number to 1
     lno<-grep('/end', eqsend, ignore.case=TRUE)
     eqsend[lno]<-''
     eqsend<-c(eqsend, '/TECHNICAL', 'ITR=1;', '/END')
     }
		list(eqsstart=eqsstart, eqsend=eqsend, ram=ram, count=d6, eqsinput=eqs)
		
	}

## Function to write EQS input file
semdiag.write.eqs<-function(eqs, par, N, P){
	#eqs<-semdiag.parse(eqsinput)
	filename<-'eqsonce.eqs'
	if (length(eqs$eqsinput)==1){
	write(eqs$eqsstart, file=filename)
	ram<-eqs$ram
	m<-nrow(ram)
	row.names<-rownames(ram)
	index<-cumsum(table(ram[,5]))
	## output equations
	cat('/EQUATIONS\n', file=filename, append=TRUE)
	## output equation
	eq.no<-length(eqs$count)
	count<-eqs$count
	## for first equation
	cat(c(' ', ram[count[1],1], ' = '), file=filename, append=TRUE)
	for (j in 1:(count[1]-1)){
		plus<-'+'
		
		if (ram[j, 4]=='0'){ 
			cat(paste(ram[j,3], ram[j,2], plus, sep=''), file=filename, append=TRUE)
		}else{ 
		   cat(paste(par[row.names[j],1], '*', ram[j,2], plus, sep=''), file=filename, append=TRUE)}
		}
		
		if (ram[count[1], 4]=='0'){ cat(paste(ram[count[1],3], ram[count[1],2], ";\n", sep=''), file=filename, append=TRUE)
		   }else{ cat(paste(par[row.names[count[1]],1], '*', ram[count[1],2], ';\n', sep=''), file=filename, append=TRUE)}
		   
	for (i in 2:eq.no){
		cat(c(' ', ram[count[i],1], ' = '), file=filename, append=TRUE)
		for (j in (count[i-1]+1):(count[i]-1)){
			plus<-'+'

			if (ram[j, 4]=='0'){ 
				cat(paste(ram[j,3], ram[j,2], plus, sep=''), file=filename, append=TRUE)
		  		}else{ 
		  			cat(paste(par[row.names[j],1], '*', ram[j,2], plus, sep=''), file=filename, append=TRUE)}
			}
			
			if (ram[count[i], 4]=='0'){ cat(paste(ram[count[i],3], ram[count[i],2], ";\n", sep=''), file=filename, append=TRUE)
		  		}else{ cat(paste(par[row.names[count[i]],1], '*', ram[count[i],2], ';\n', sep=''), file=filename, append=TRUE)}

		}
	## output variances
	cat('/VARIANCES\n', file=filename, append=TRUE)
	for (i in (index[1]+1):index[2]){
		cat(c(ram[i,1], ' = '), file=filename, append=TRUE)
		if (ram[i, 4]=='0'){ cat(c(ram[i,3], ';\n'), file=filename, append=TRUE)
			}else{ 
				par.value<-as.numeric(par[row.names[i],1])
				#if (par.value<0) par.value<-0
				cat(c(par.value, '*;\n'), file=filename, append=TRUE)
				}
			}
		
	if (length(index)>2){	
		cat('/COVARIANCES\n', file=filename, append=TRUE)	
		for (i in (index[2]+1):index[3]){
			cat(c(' ', ram[i,1], ',', ram[i,2], ' = '), file=filename, append=TRUE)
			if (ram[i, 4]=='0'){ cat(c(ram[i,3], ';\n'), file=filename, append=TRUE)
				}else{ cat(c(par[row.names[i],1], '*;\n'), file=filename, append=TRUE)}
			}
		}
	
	## nonlinear constraints for variance parameters
	cat('/INEQUALITIES\n', file=filename, append=TRUE)
	for (i in (index[1]+1):index[2]){
		cat(c('(', ram[i,1], ',', ram[i,1],') '), file=filename, append=TRUE)
		if (ram[i, 4]=='1') cat('GT -99999999;\n', file=filename, append=TRUE)		}
		
	write(eqs$eqsend, file=filename, append=TRUE)
	
	## replace +- with -
	temp.eqs<-readLines('eqsonce.eqs', warn=FALSE)
	for (i in 1:length(temp.eqs)) temp.eqs[i]<-gsub('+-', '-', temp.eqs[i], fixed=TRUE)
	write(temp.eqs, file=filename)
	}else{
		## for model input in R
		cat('/TITLE\n', file=filename)
		cat('EQS input generated by the semdiag package;\n', file=filename, append=TRUE)
		cat('/SPECIFICATIONS\n', file=filename, append=TRUE)
		if (missing(N)) stop('No sample size is supplied')
		cat(c('CASES=',N,'; VARIABLES=', P, '; METHOD=ML;\n'), file=filename, append=TRUE)
		cat('MATRIX=COVARIANCE; ANALYSIS=COVARIANCE;\n', file=filename, append=TRUE)
		
	ram<-eqs$ram
	m<-nrow(ram)
	row.names<-rownames(ram)
	index<-cumsum(table(ram[,5]))
	## output equations
	cat('/EQUATIONS\n', file=filename, append=TRUE)
	## output equation
	eq.no<-length(eqs$count)
	count<-eqs$count
	## for first equation
	cat(c(' ', ram[count[1],1], ' = '), file=filename, append=TRUE)
	for (j in 1:(count[1]-1)){
		if (ram[j, 4]=='0'){ cat(paste(ram[j,3], ram[j,2], '+', sep=''), file=filename, append=TRUE)
		   }else{ cat(paste(par[row.names[j],1], '*', ram[j,2], '+', sep=''), file=filename, append=TRUE)}
		}
		if (ram[count[1], 4]=='0'){ cat(paste(ram[count[1],3], ram[count[1],2], ";\n", sep=''), file=filename, append=TRUE)
		   }else{ cat(paste(par[row.names[count[1]],1], '*', ram[count[1],2], ';\n', sep=''), file=filename, append=TRUE)}
		   
	for (i in 2:eq.no){
		cat(c(' ', ram[count[i],1], ' = '), file=filename, append=TRUE)
		for (j in (count[i-1]+1):(count[i]-1)){
			if (ram[j, 4]=='0'){ cat(paste(ram[j,3], ram[j,2], '+', sep=''), file=filename, append=TRUE)
		  		}else{ cat(paste(par[row.names[j],1], '*', ram[j,2], '+', sep=''), file=filename, append=TRUE)}
			}
			if (ram[count[i], 4]=='0'){ cat(paste(ram[count[i],3], ram[count[i],2], ";\n", sep=''), file=filename, append=TRUE)
		  		}else{ cat(paste(par[row.names[count[i]],1], '*', ram[count[i],2], ';\n', sep=''), file=filename, append=TRUE)}

		}
	## output variances
	cat('/VARIANCES\n', file=filename, append=TRUE)
	for (i in (index[1]+1):index[2]){
		cat(c(ram[i,1], ' = '), file=filename, append=TRUE)
		if (ram[i, 4]=='0'){ cat(c(ram[i,3], ';\n'), file=filename, append=TRUE)
			}else{ 
				par.value<-as.numeric(par[row.names[i],1])
				if (par.value<0) par.value=0
				cat(c(par.value, '*;\n'), file=filename, append=TRUE)}
				}
		
	if (length(index)>2){	
		cat('/COVARIANCES\n', file=filename, append=TRUE)	
		for (i in (index[2]+1):index[3]){
			cat(c(' ', ram[i,1], ',', ram[i,2], ' = '), file=filename, append=TRUE)
			if (ram[i, 4]=='0'){ cat(c(ram[i,3], ';\n'), file=filename, append=TRUE)
				}else{ cat(c(par[row.names[i],1], '*;\n'), file=filename, append=TRUE)}
			}
		}
	if (!is.na(eqs$eqsend)) write(eqs$eqsend, file=filename, append=TRUE)
	cat('/TECHNICAL\n', file=filename, append=TRUE)
	cat('ITR=1;\n', file=filename, append=TRUE)
	cat('/OUTPUT\n', file=filename, append=TRUE)
	cat('CODEBOOK; DATA=\"eqsonce\".ETS; PARAMETERS ESTIMATES; STANDARD ERRORS; LISTING;\n \\END', file=filename, append=TRUE)
		}
	}
	
	

## construct lisrel model matrices from EQS output

## model: EQS output

semdiag.eqs.lisrel<-function(model){

	p<-model$model.desc$values[2] ## number of observed varialbes

	## identify x, y, eta, xi

	var.ind<-colnames(model$Phi)

	## take out the first letter

	var.ind.sub<-substr(var.ind,1,1)

	

	if ('D' %in% var.ind.sub){

		## 

		xi<-var.ind[var.ind.sub=='F']    ## independent factors

		d.phi<-var.ind[var.ind.sub=='D']   ## dependent factors

		e.phi<-var.ind[var.ind.sub=='E']   ## error for observed variables



		## identify eta

		var.ind.gamma<-rownames(model$Gamma)

		var.ind.gamma.sub<-substr(var.ind.gamma, 1, 1)

		eta<-var.ind.gamma[var.ind.gamma.sub=='F']



		## identify x

		xy<-var.ind.gamma[var.ind.gamma.sub=='V']



		## take out a submatrix

		gamma<-model$Gamma[xy, xi]

		gamma<-as.matrix(gamma)

		sum.gamma<-apply(gamma,1,sum)

		x<-xy[!(sum.gamma==0)]

		y<-xy[(sum.gamma==0)]



		nx<-length(x)

		ny<-length(y)

		## get ex and ey

		ex<-e.phi[1:nx]

		ey<-e.phi[(nx+1):p]



		## build the lisrel matrix

		lambda.x<-model$Gamma[x, xi]

		lambda.x<-as.matrix(lambda.x)

		lambda.x[lambda.x==-1]<-1



		lambda.y<-model$Beta[y,eta]

		lambda.y<-as.matrix(lambda.y)

		lambda.y[lambda.y==-1]<-1



		beta<-model$Beta[eta, eta]

		beta<-as.matrix(beta)

		diag(beta)<-0



		gamma<-model$Gamma[eta, xi]

		gamma<-as.matrix(gamma)



		theta.x<-model$Phi[ex,ex]

		theta.y<-model$Phi[ey,ey]



		phi<-model$Phi[xi,xi]

		phi<-as.matrix(phi)

		psi<-model$Phi[d.phi, d.phi]

		psi<-as.matrix(psi)

		

		A<-diag(nrow(beta))

		AB<-A-beta

		

		## computing Sigma

		ABin<-solve(AB)

		sigyy<-lambda.y%*%ABin%*%(gamma%*%phi%*%t(gamma)+psi)%*%t(ABin)%*%t(lambda.y)+theta.y

		sigxy<-lambda.x%*%phi%*%t(gamma)%*%t(ABin)%*%t(lambda.y)

		sigxx<-lambda.x%*%phi%*%t(lambda.x)+theta.x

		sigma0<-rbind(cbind(sigxx,sigxy),cbind(t(sigxy),sigyy))

	

		## put all matrices in a list for future use

		list(lambda.x=lambda.x, lambda.y=lambda.y, beta=beta, gamma=gamma, phi=phi, psi=psi, theta.x=theta.x, theta.y=theta.y, sigma0=sigma0, eqsresults=model, model=2)

	}else{

		## factor model

		xi<-var.ind[var.ind.sub=='F']    ## independent factors

	   e.phi<-var.ind[var.ind.sub=='E']   ## error for observed variables



		## identify eta

		var.ind.gamma<-rownames(model$Gamma)

		var.ind.gamma.sub<-substr(var.ind.gamma, 1, 1)

		x<-var.ind.gamma[var.ind.gamma.sub=='V']



		nx<-length(x)

		ex<-e.phi[1:nx]



		## build the lisrel matrix

		lambda.x<-model$Gamma[x, xi]

		lambda.x<-as.matrix(lambda.x)

		lambda.x[lambda.x==-1]<-1



		theta.x<-model$Phi[ex,ex]

	

		phi<-model$Phi[xi,xi]

		phi<-as.matrix(phi)

		

		sigma0<-lambda.x%*%phi%*%t(lambda.x)+theta.x

	

		## put all matrices in a list for future use

		list(lambda.x=lambda.x, phi=phi, theta.x=theta.x, sigma0=sigma0, eqsresults=model, model=1)

	}	

}





## Function estimate the saturated meand and covariance using Huber weight

## x: data

## varphi: proportion of downweight

semdiag.musig<-function(x,varphi, max_it=1000){

	if(!is.matrix(x)) x<-data.matrix(x)

	ep<-.000001

	n<-nrow(x)

	p<-ncol(x)	

	#Tunning parameters in Huber weight;

	prob<-1-varphi

	chip<-qchisq(prob,p)

	ck<-sqrt(chip) # this is rho, the d start to downweight 

	cbeta<-(p*pchisq(chip,p+2)+chip*(1-prob))/p



   #initial values 

	mu0<-colMeans(x) 

	sigma0<-cov(x)

	n_it<-1

	dt<-1

	while (dt>ep && n_it <= max_it){ 

		sigin<-solve(sigma0)

		sumw1<-0.0

		mu<-rep(0,p) # a p vector 

		sigma<-matrix(0,p,p)	

    	

   	for (i in 1:n) {

			xii<-x[i,] 

  			xi0<-xii-mu0 

  			di2<-xi0%*%sigin%*%xi0

  			di<-sqrt(di2)  		

			#Huber weight functions

  			if (di<=ck) {

     			wi1<-1.0

     			wi2<-1.0/cbeta

			}else {

     			wi1<-ck/di

     			wi2<-wi1^2/cbeta

 			}

  			sumw1<-sumw1+wi1

  			mu<-mu+c(wi1)*xii

  			sigma<-sigma+c(wi2)*xi0%*%t(xi0)

   } # end of loop n 

   

	mu1<-mu/sumw1

	sigma1<-sigma/n

    

	dt<-max(c(max(abs(mu1-mu0)), max(abs(sigma1-sigma0))));

	n_it<-n_it+1

	mu0<-mu1

	sigma0<-sigma1      	   

	} 

	list(mu=mu1,sigma=sigma1, n_it=n_it)

}

	



## Function to calculate M-distance for factors based on Bartlett-factor score

## x: data matrix

## mu: saturated mean from the robust estimation

## model: eqs output from running eqs

semdiag.mdist.f<-function(x, mu, lisrel){

	x<-data.matrix(x)

	n<-nrow(x)

	p<-ncol(x)

	mu<-matrix(mu, 1, length(mu))



	if (lisrel$model==2){

		## SEM model

		lambx<-lisrel$lambda.x

		lamby<-lisrel$lambda.y

		B<-lisrel$beta

		A<-diag(nrow(B))

		AB<-A-B

		gamma<-lisrel$gamma

		phi<-lisrel$phi

		psi<-lisrel$psi

		thdx<-lisrel$theta.x

		they<-lisrel$theta.y

		

		## computing Sigma

		ABin<-solve(AB)

		sigyy<-lamby%*%ABin%*%(gamma%*%phi%*%t(gamma)+psi)%*%t(ABin)%*%t(lamby)+they

		sigxy<-lambx%*%phi%*%t(gamma)%*%t(ABin)%*%t(lamby)

		sigxx<-lambx%*%phi%*%t(lambx)+thdx

		sigma0<-rbind(cbind(sigxx,sigxy),cbind(t(sigxy),sigyy))

		

		## constructing distance based on factor model notation

		lamb<-rbind(lambx, lamby%*%ABin%*%gamma) ## factor loading matrix

		px<-nrow(thdx)

		py<-nrow(they)

		psi_new<-rbind(cbind(thdx,matrix(0,px,py) ),cbind(matrix(0,py,px),they+lamby%*%ABin%*%psi%*%t(lamby%*%ABin)))

		psiin<-solve(psi_new);

		lpsiin<-t(lamb)%*%psiin;

		sig_f<-lpsiin%*%sigma0%*%t(lpsiin)

		sig_fin<-solve(sig_f)

		Md_f<-matrix(0,n,1)

		for (i in 1:n){

			xi0<-x[i,]-mu

			fi0<-lpsiin%*%t(xi0)

			d_fi2<-t(fi0)%*%sig_fin%*%fi0

			Md_f[i,1]<-sqrt(d_fi2)

      	}		

		}

		if (lisrel$model==1){

			## for a factor model

			lamb<-lisrel$lambda.x

			phi<-lisrel$phi

			psi<-lisrel$theta.x

			

			sigma0<-lamb%*%phi%*%t(lamb)+psi

			psiin<-solve(psi)

			lpsiin<-t(lamb)%*%psiin;

			sig_f<-lpsiin%*%sigma0%*%t(lpsiin)

			sig_fin<-solve(sig_f)

			Md_f<-matrix(0,n,1);

			for (i in 1:n){

				xi0<-x[i,]-mu

				fi0<-lpsiin%*%t(xi0)

				d_fi2<-t(fi0)%*%sig_fin%*%fi0

				Md_f[i,1]<-sqrt(d_fi2)

      		}			

		}

		return(Md_f)

	}



## Function to calculate M-distance for factors based on Bartlett-factor score

## Based on all latent factors

## x: data matrix

## mu: saturated mean from the robust estimation

## model: eqs output from running eqs

semdiag.mdist.f1<-function(x, mu, lisrel){

	x<-data.matrix(x)

	n<-nrow(x)

	p<-ncol(x)

	mu<-matrix(mu, 1, length(mu))

	

	if (lisrel$model==2){

## SEM model

		lambx<-lisrel$lambda.x

		lamby<-lisrel$lambda.y

		B<-lisrel$beta

		A<-diag(nrow(B))

		AB<-A-B

		gamma<-lisrel$gamma

		phi<-lisrel$phi

		psi<-lisrel$psi

		thdx<-lisrel$theta.x

		they<-lisrel$theta.y

		

## computing Sigma

		ABin<-solve(AB)

		sigyy<-lamby%*%ABin%*%(gamma%*%phi%*%t(gamma)+psi)%*%t(ABin)%*%t(lamby)+they

		sigxy<-lambx%*%phi%*%t(gamma)%*%t(ABin)%*%t(lamby)

		sigxx<-lambx%*%phi%*%t(lambx)+thdx

		sigma0<-rbind(cbind(sigxx,sigxy),cbind(t(sigxy),sigyy))

		

## constructing distance based on factor model notation

		px<-nrow(thdx)

		pb<-nrow(B)

		py<-nrow(they)

		lamb<-rbind(cbind(lambx,matrix(0,px,pb)),cbind(lamby%*%ABin%*%gamma, lamby%*%ABin) ) ## factor loading matrix		

		psi_new<-rbind(cbind(thdx,matrix(0,px,py) ),cbind(matrix(0,py,px),they ))

		psiin<-solve(psi_new);

		lpsiin<-t(lamb)%*%psiin;

		sig_f<-lpsiin%*%sigma0%*%t(lpsiin)

		sig_fin<-solve(sig_f)

		Md_f<-matrix(0,n,1)

		for (i in 1:n){

			xi0<-x[i,]-mu

			fi0<-lpsiin%*%t(xi0)

			d_fi2<-t(fi0)%*%sig_fin%*%fi0

			Md_f[i,1]<-sqrt(d_fi2);

      	}		

	}

	if (lisrel$model==1){

## for a factor model

		lamb<-lisrel$lambda.x

		phi<-lisrel$phi

		psi<-lisrel$theta.x

		

		sigma0<-lamb%*%phi%*%t(lamb)+psi

		psiin<-solve(psi)

		lpsiin<-t(lamb)%*%psiin;

		sig_f<-lpsiin%*%sigma0%*%t(lpsiin)

		sig_fin<-solve(sig_f)

		Md_f<-matrix(0,n,1);

		for (i in 1:n){

			xi0<-x[i,]-mu

			fi0<-lpsiin%*%t(xi0)

			d_fi2<-t(fi0)%*%sig_fin%*%fi0

			Md_f[i,1]<-sqrt(d_fi2);

		}			

	}

	return(Md_f)

}



#--------------------------------------------------------------------------*;

#### Ac is an orthogonal complement of A used in calculating d_r;

#--------------------------------------------------------------------------*;

semdiag.orthog<-function(A){

	p<-nrow(A)

	q<-ncol(A)

	pmq<-p-q

	At<-t(A)

	Pa<-A%*%solve(At%*%A)%*%At

  I_p<-diag(p)

	Qa<-I_p-Pa

	Qv<-eigen(Qa)

	Qvec<-Qv$vec

	Qvec[,1:pmq]

}





## Function to calculate M-distance for residuals

## x: data matrix

## mu: saturated mean from the robust estimation

## model: eqs output from running eqs

semdiag.mdist.r<-function(x, mu, lisrel){

	if(!is.matrix(x)) x<-data.matrix(x)

	n<-nrow(x)

	p<-ncol(x)

	mu<-matrix(mu, 1, length(mu))

	if (lisrel$model==2){

		## SEM model

		lambx<-lisrel$lambda.x

		lamby<-lisrel$lambda.y

		B<-lisrel$beta

		A<-diag(nrow(B))

		AB<-A-B

		gamma<-lisrel$gamma

		phi<-lisrel$phi

		psi<-lisrel$psi

		thdx<-lisrel$theta.x

		they<-lisrel$theta.y

		

		## computing Sigma

		ABin<-solve(AB)

		sigyy<-lamby%*%ABin%*%(gamma%*%phi%*%t(gamma)+psi)%*%t(ABin)%*%t(lamby)+they

		sigxy<-lambx%*%phi%*%t(gamma)%*%t(ABin)%*%t(lamby)

		sigxx<-lambx%*%phi%*%t(lambx)+thdx

		sigma<-rbind(cbind(sigxx,sigxy),cbind(t(sigxy),sigyy))

		

		## constructing distance based on factor model notation

		lamb<-rbind(lambx, lamby%*%ABin%*%gamma) ## factor loading matrix

		lambt<-t(lamb)

		px<-nrow(thdx)

		py<-nrow(they)

		psi_new<-rbind(cbind(thdx,matrix(0,px,py) ),cbind(matrix(0,py,px),they+lamby%*%ABin%*%psi%*%t(ABin)%*%t(lamby)))

		psiin<-solve(psi_new);

		Pmat<-diag(p)-lamb%*%solve(lambt%*%psiin%*%lamb)%*%lambt%*%psiin

		Covy<-Pmat%*%sigma%*%t(Pmat)

		Av<-psiin%*%lamb

		Avc<-semdiag.orthog(Av)

		cove<-t(Avc)%*%Covy%*%Avc

		covein<-solve(cove)



		dr2<-matrix(0,n,1)



		for (i in 1:n){

	 		xi<-x[i,]-mu

			yi<-Pmat%*%t(xi) # yi is a col vector

			ei<-t(Avc)%*%yi # a col vector 

			dr2[i]<-t(ei)%*%covein%*%ei

			}

		}

	

	if (lisrel$model==1){

		## for a factor model

		lamb<-lisrel$lambda.x

		phi<-lisrel$phi

		psi<-lisrel$theta.x

			

		psiin<-solve(psi)

		lambt<-t(lamb)

		sigma<-lamb%*%phi%*%lambt+psi

		I_p<-diag(p)

		pmat<-I_p-lamb%*%solve(lambt%*%psiin%*%lamb)%*%lambt%*%psiin

		covy<-pmat%*%sigma%*%t(pmat)

		Av<-psiin%*%lamb

		Avc<-semdiag.orthog(Av)

		Avct<-t(Avc)

		cove<-Avct%*%covy%*%Avc

		covein<-solve(cove)



		dr2<-matrix(0,n,1)

		for (i in 1:n){

			xi<-x[i,]-mu

			yi<-pmat%*%t(xi)

			ei<-Avct%*%yi

			dr2[i]<-t(ei)%*%covein%*%ei

			}			

		}

	dr2

	}



## Function to calculate M-distance for residuals using all factors

## x: data matrix

## mu: saturated mean from the robust estimation

## model: eqs output from running eqs

semdiag.mdist.r1<-function(x, mu, lisrel){

	if(!is.matrix(x)) x<-data.matrix(x)

	n<-nrow(x)

	p<-ncol(x)

	mu<-matrix(mu, 1, length(mu))

	if (lisrel$model==2){

## SEM model

		lambx<-lisrel$lambda.x

		lamby<-lisrel$lambda.y

		B<-lisrel$beta

		A<-diag(nrow(B))

		AB<-A-B

		gamma<-lisrel$gamma

		phi<-lisrel$phi

		psi<-lisrel$psi

		thdx<-lisrel$theta.x

		they<-lisrel$theta.y

		

## computing Sigma

		ABin<-solve(AB)

		sigyy<-lamby%*%ABin%*%(gamma%*%phi%*%t(gamma)+psi)%*%t(ABin)%*%t(lamby)+they

		sigxy<-lambx%*%phi%*%t(gamma)%*%t(ABin)%*%t(lamby)

		sigxx<-lambx%*%phi%*%t(lambx)+thdx

		sigma<-rbind(cbind(sigxx,sigxy),cbind(t(sigxy),sigyy))

		

## constructing distance based on factor model notation

		px<-nrow(thdx)

		pb<-nrow(B)

		py<-nrow(they)

		lamb<-rbind(cbind(lambx,matrix(0,px,pb)),cbind(lamby%*%ABin%*%gamma, lamby%*%ABin) ) ## factor loading matrix		

		psi_new<-rbind(cbind(thdx,matrix(0,px,py) ),cbind(matrix(0,py,px),they ))

		lambt<-t(lamb)

		psiin<-solve(psi_new);

		Pmat<-diag(p)-lamb%*%solve(lambt%*%psiin%*%lamb)%*%lambt%*%psiin

		Covy<-Pmat%*%sigma%*%t(Pmat)

		Av<-psiin%*%lamb

		Avc<-semdiag.orthog(Av)

		cove<-t(Avc)%*%Covy%*%Avc

		covein<-solve(cove)

		

		dr2<-matrix(0,n,1)

		

		for (i in 1:n){

	 		xi<-x[i,]-mu

			yi<-Pmat%*%t(xi) # yi is a col vector

			ei<-t(Avc)%*%yi # a col vector 

			dr2[i]<-t(ei)%*%covein%*%ei

		}

	}

	

	if (lisrel$model==1){

## for a factor model

		lamb<-lisrel$lambda.x

		phi<-lisrel$phi

		psi<-lisrel$theta.x

		

		psiin<-solve(psi)

		lambt<-t(lamb)

		sigma<-lamb%*%phi%*%lambt+psi

		I_p<-diag(p)

		pmat<-I_p-lamb%*%solve(lambt%*%psiin%*%lamb)%*%lambt%*%psiin

		covy<-pmat%*%sigma%*%t(pmat)

		Av<-psiin%*%lamb

		Avc<-semdiag.orthog(Av)

		Avct<-t(Avc)

		cove<-Avct%*%covy%*%Avc

		covein<-solve(cove)

		

		dr2<-matrix(0,n,1)

		for (i in 1:n){

			xi<-x[i,]-mu

			yi<-pmat%*%t(xi)

			ei<-Avct%*%yi

			dr2[i]<-t(ei)%*%covein%*%ei

		}			

	}

	dr2

}



##############################################

## Stacking lower triange of a matrix to    ##

##   a vector                               ##

## Function rsem_vech                       ##

##############################################

semdiag.vech<-function(x){

	t(x[!upper.tri(x)])

}



##############################################

## Generate a duplication matrix            ##

## Function rsem_DP                         ##

##############################################

semdiag.DP <- function(x){ 	

	mat <- diag(x)

 	index <- seq(x*(x+1)/2)

 	mat[ lower.tri( mat , TRUE ) ] <- index

 	mat[ upper.tri( mat ) ] <- t( mat )[ upper.tri( mat ) ]

 	outer(c(mat), index , function( x , y ) ifelse(x==y, 1, 0 ) )

}





#--------------------------------------------------------------------------*;

# parameter estimation by H_r(varphi) using iteratively reweighted least squares;

##Huber-type weight

#--------------------------------------------------------------------------*;

## theta0: starting values

## x: data

## varphi: 

semdiag.robfit<-function(lisrel0, x, q, varphi, EQSmodel, EQSdata, max_it=1000, EQSprog='C:/Progra~1/EQS61/WINEQS',serial=1234){
	## parse EQS input file
	eqs<-semdiag.parse(EQSmodel)
	

	ep<-.0001

	n<-nrow(x)

	p<-ncol(x)

	p1<-p-q 

	x<-data.matrix(x)

	#Tunning parameters in Huber weight;

    prob<-1-varphi

	chip1<-qchisq(prob, p1)

	ck<-sqrt(chip1)

	cbeta<-( p1*pchisq(chip1,p1+2)+ chip1*(1-prob) )/p1



	mu0<-colMeans(x)

	theta0<-lisrel0$eqsresults$par.table[,1]

	sigma0<-lisrel0$sigma0



	Md_r<-matrix(0,n,3)	



	n_it<-1

	dt<-1

	while (dt>ep && n_it <= max_it){

		dr2<-semdiag.mdist.r(x,mu0,lisrel0)

		sumw1<-0.0

		mu<-rep(0,p)

		sigma<-matrix(0,p,p)



		for (i in 1:n) {

			xii<-x[i,] 

  			xi0<-xii-mu0 

  			di2<-dr2[i]

  			di<-sqrt(di2)

  		

			#Huber weight functions

  			if (di<=ck) {

     			wi1<-1.0

     			wi2<-1.0/cbeta

			}else {

     			wi1<-ck/di

     			wi2<-wi1^2/cbeta

 			}



			Md_r[i,1]<-i

  			Md_r[i,2]<-di

  			Md_r[i,3]<-wi2

  			sumw1<-sumw1+wi1

  			mu<-mu+c(wi1)*xii

  			sigma<-sigma+c(wi2)*xi0%*%t(xi0)  		

		} # end of loop n 



		mu1<-mu/sumw1

		sigma1<-sigma/n

	

		## write data and run EQS for parameter estimates

		write.table(sigma1, EQSdata, row.names=FALSE, col.names=FALSE)
		
		## write EQS input file with updated starting values
		semdiag.write.eqs(eqs, lisrel0$eqsresults$par.table)

		model1<-semdiag.run.eqs(EQSprog, 'eqsonce.eqs', serial=serial)

	

		## estimated model parameters

		theta1<-model1$par.table[,1]

	

		dt<-max(c(max(abs(mu1-mu0)), max(abs(theta1-theta0))))	

		n_it<-n_it+1

		mu0<-mu1

		lisrel0<-semdiag.eqs.lisrel(model1)

		theta0<-lisrel0$eqsresults$par.table[,1]

		# print(dt)

  } # end of loop k              

 list(Md_r=Md_r,mu1=mu1,eqs=model1, n_it=n_it)                  

}	



#--------------------------------------------------------------------------*;

# parameter estimation by H_r(varphi) using iteratively reweighted least squares;

##Huber-type weight

#--------------------------------------------------------------------------*;

## theta0: starting values

## x: data

## varphi: 

semdiag.robfit1<-function(lisrel0, x, q, varphi, EQSmodel, EQSdata, max_it=1000, EQSprog='C:/Progra~1/EQS61/WINEQS',serial=1234){

	## parse EQS input file
	eqs<-semdiag.parse(EQSmodel)
	
	ep<-.0001

	n<-nrow(x)

	p<-ncol(x)

	p1<-p-q 

	x<-data.matrix(x)

#Tunning parameters in Huber weight;

    prob<-1-varphi

	chip1<-qchisq(prob, p1)

	ck<-sqrt(chip1)

	cbeta<-( p1*pchisq(chip1,p1+2)+ chip1*(1-prob) )/p1

	

	mu0<-colMeans(x)

	theta0<-lisrel0$eqsresults$par.table[,1]

	sigma0<-lisrel0$sigma0

	

	Md_r<-matrix(0,n,3)	

	

	n_it<-1

	dt<-1

	while (dt>ep && n_it <= max_it){

		dr2<-semdiag.mdist.r1(x,mu0,lisrel0)

		sumw1<-0.0

		mu<-rep(0,p)

		sigma<-matrix(0,p,p)

		

		for (i in 1:n) {

			xii<-x[i,] 

  			xi0<-xii-mu0 

  			di2<-dr2[i]

  			di<-sqrt(di2)

			

#Huber weight functions

  			if (di<=ck) {

     			wi1<-1.0

     			wi2<-1.0/cbeta

			}else {

     			wi1<-ck/di

     			wi2<-wi1^2/cbeta

 			}

			

			Md_r[i,1]<-i

  			Md_r[i,2]<-di

  			Md_r[i,3]<-wi2

  			sumw1<-sumw1+wi1

  			mu<-mu+c(wi1)*xii

  			sigma<-sigma+c(wi2)*xi0%*%t(xi0)  		

		} # end of loop n 

		

		mu1<-mu/sumw1

		sigma1<-sigma/n

		

## write data and run EQS for parameter estimates

		write.table(sigma1, EQSdata, row.names=FALSE, col.names=FALSE)

		## write EQS input file with updated starting values
		semdiag.write.eqs(eqs, lisrel0$eqsresults$par.table)
		model1<-semdiag.run.eqs(EQSprog, 'eqsonce.eqs', serial=serial)	

##model1<-run.eqs("C:/Progra~1/EQS61/WINEQS", 'c:/rsem/lisrel.eqs', serial=1234)

		

## estimated model parameters

		theta1<-model1$par.table[,1]

		

		dt<-max(c(max(abs(mu1-mu0)), max(abs(theta1-theta0))))	

		n_it<-n_it+1

		mu0<-mu1

		lisrel0<-semdiag.eqs.lisrel(model1)

		theta0<-lisrel0$eqsresults$par.table[,1]

	} # end of loop k              

	list(Md_r=Md_r,mu1=mu1,eqs=model1, n_it=n_it)                  

}	



semdiag.eqs<-function(res){

	pval <- res$pval['PVAL',1]

	chi <- res$fit.indices['CHI',1]

	df <- res$model.info['DF ',1]   

	z <- res$par.table[, 1]/res$par.table[, 2]

	par.est <- cbind(res$par.table[, c(1, 2)], z)

	list(chi=chi, df=df, pval=pval,par=par.est)

	}

semdiag.lisrel<-function(model){
	var.names<-model$var.names
	nvar<-length(var.names)  ## total number of variables
	p<-model$n
	ram<-model$ram
	obs.var<-rownames(model$C)
	latent.names<-var.names[(p+1):nvar]
	obs.names<-var.names[1:p]
	q<-length(latent.names)
	
	## single-headed arrows to identify ind and dep variables
	ram.s<-ram[ram[,1]==1, ]
	dep.var<-unique(ram.s[,2])
	
	## For latent variables
	dep.lat.var<-dep.var[dep.var>p]-p
	
	if (length(dep.lat.var)>0){
		dep.lat.names<-latent.names[dep.lat.var]
		index<-1:q
		index.ind<-index[-dep.lat.var]
		ind.lat.names<-latent.names[index.ind]
	
		## Observed variables related to depedent and independent latent variables
		## index for independent latent variables
		index.ind.lat<-index.ind+p
		ram.l<-ram[ram[,1]==1 & ram[,2]<=p, ]
		ram.l.temp<-ram.l[ ram.l[,3] %in% index.ind.lat, ]
		ind.obs.var<-unique(ram.l.temp[,2])
		ind.obs.names<-obs.var[ind.obs.var]
		## observed variables related to the dependent latent variables
		index<-1:p
		index.ind<-index[-ind.obs.var]
		dep.obs.names<-obs.var[index.ind]
	
		lambda.x<-as.matrix(model$A[ind.obs.names, ind.lat.names])
		lambda.y<-as.matrix(model$A[dep.obs.names, dep.lat.names])
	
		beta<-as.matrix(model$A[dep.lat.names,dep.lat.names])
		gamma<-as.matrix(model$A[dep.lat.names,ind.lat.names])
	
		theta.x<-as.matrix(model$P[ind.obs.names,ind.obs.names])
		theta.y<-as.matrix(model$P[dep.obs.names,dep.obs.names])
	
		psi<-as.matrix(model$P[dep.lat.names,dep.lat.names])
		phi<-as.matrix(model$P[ind.lat.names,ind.lat.names])
		list(lambda.x=lambda.x, lambda.y=lambda.y, beta=beta, gamma=gamma, phi=phi, psi=psi, theta.x=theta.x, theta.y=theta.y, sigma0=model$C, semoject=model, coeff=model$coeff)
	}else{
		## for factor model
		lambda.x<-as.matrix(model$A[obs.names, latent.names])
		theta.x<-as.matrix(model$P[obs.names, obs.names])
		phi<-as.matrix(model$P[latent.names,latent.names])
		list(lambda.x=lambda.x, phi=phi, theta.x=theta.x, sigma0=model$C, semoject=model, coeff=model$coeff)
	}	
	}

## Function to calculate M-distance for factors based on Bartlett-factor score
## x: data matrix
## mu: saturated mean from the robust estimation
## lisrel: lisrel notation model
semdiag.mdist.f.r<-function(x, mu, lisrel){
	x<-data.matrix(x)
	n<-nrow(x)
	p<-ncol(x)
	mu<-t(mu)

	if (!is.null(lisrel$lambda.y)){
		## SEM model
		lambx<-lisrel$lambda.x
		lamby<-lisrel$lambda.y
		B<-lisrel$beta
		A<-diag(nrow(B))
		AB<-A-B
		gamma<-lisrel$gamma
		phi<-lisrel$phi
		psi<-lisrel$psi
		thdx<-lisrel$theta.x
		they<-lisrel$theta.y
		
		## computing Sigma
		sigma0<-lisrel$sigma0
		ABin<-solve(AB)
		## constructing distance based on factor model notation
		lamb<-rbind(lambx, lamby%*%ABin%*%gamma) ## factor loading matrix
		px<-nrow(thdx)
		py<-nrow(they)
		psi_new<-rbind(cbind(thdx,matrix(0,px,py) ),cbind(matrix(0,py,px),they+lamby%*%ABin%*%psi%*%t(lamby%*%ABin)))
		psiin<-solve(psi_new);
		lpsiin<-t(lamb)%*%psiin;
		sig_f<-lpsiin%*%sigma0%*%t(lpsiin)
		sig_fin<-solve(sig_f)
		Md_f<-matrix(0,n,1)
		for (i in 1:n){
			xi0<-x[i,]-mu
			fi0<-lpsiin%*%t(xi0)
			d_fi2<-t(fi0)%*%sig_fin%*%fi0
			Md_f[i,1]<-sqrt(d_fi2)
      	}		
	}else{
			## for a factor model
			lamb<-lisrel$lambda.x
			phi<-lisrel$phi
			psi<-lisrel$theta.x
			
			sigma0<-sigma0<-lisrel$sigma0
			psiin<-solve(psi)
			lpsiin<-t(lamb)%*%psiin;
			sig_f<-lpsiin%*%sigma0%*%t(lpsiin)
			sig_fin<-solve(sig_f)
			Md_f<-matrix(0,n,1);
			for (i in 1:n){
				xi0<-x[i,]-mu
				fi0<-lpsiin%*%t(xi0)
				d_fi2<-t(fi0)%*%sig_fin%*%fi0
				Md_f[i,1]<-sqrt(d_fi2)
      		}			
		}
		return(Md_f)
	}


## Function to calculate M-distance for factors based on Bartlett-factor score
## Based on all latent factors
## x: data matrix
## mu: saturated mean from the robust estimation
## lisrel: lisrel notation model
semdiag.mdist.f1.r<-function(x, mu, lisrel){
	x<-data.matrix(x)
	n<-nrow(x)
	p<-ncol(x)
	mu<-t(mu)
	
	if (!is.null(lisrel$lambda.y)){
## SEM model
		lambx<-lisrel$lambda.x
		lamby<-lisrel$lambda.y
		B<-lisrel$beta
		A<-diag(nrow(B))
		AB<-A-B
		gamma<-lisrel$gamma
		phi<-lisrel$phi
		psi<-lisrel$psi
		thdx<-lisrel$theta.x
		they<-lisrel$theta.y
		
		sigma0<-lisrel$sigma0
		ABin<-solve(AB)		
## constructing distance based on factor model notation
		px<-nrow(thdx)
		pb<-nrow(B)
		py<-nrow(they)
		lamb<-rbind(cbind(lambx,matrix(0,px,pb)),cbind(lamby%*%ABin%*%gamma, lamby%*%ABin) ) ## factor loading matrix		
		psi_new<-rbind(cbind(thdx,matrix(0,px,py) ),cbind(matrix(0,py,px),they ))
		psiin<-solve(psi_new);
		lpsiin<-t(lamb)%*%psiin;
		sig_f<-lpsiin%*%sigma0%*%t(lpsiin)
		sig_fin<-solve(sig_f)
		Md_f<-matrix(0,n,1)
		for (i in 1:n){
			xi0<-x[i,]-mu
			fi0<-lpsiin%*%t(xi0)
			d_fi2<-t(fi0)%*%sig_fin%*%fi0
			Md_f[i,1]<-sqrt(d_fi2);
      	}		
	}else{
## for a factor model
		lamb<-lisrel$lambda.x
		phi<-lisrel$phi
		psi<-lisrel$theta.x
		
		sigma0<-lisrel$sigma0
		psiin<-solve(psi)
		lpsiin<-t(lamb)%*%psiin;
		sig_f<-lpsiin%*%sigma0%*%t(lpsiin)
		sig_fin<-solve(sig_f)
		Md_f<-matrix(0,n,1);
		for (i in 1:n){
			xi0<-x[i,]-mu
			fi0<-lpsiin%*%t(xi0)
			d_fi2<-t(fi0)%*%sig_fin%*%fi0
			Md_f[i,1]<-sqrt(d_fi2);
		}			
	}
	return(Md_f)
}



## Function to calculate M-distance for residuals
## x: data matrix
## mu: saturated mean from the robust estimation
## lisrel: lisrel notation model
semdiag.mdist.r.r<-function(x, mu, lisrel){
	if(!is.matrix(x)) x<-data.matrix(x)
	n<-nrow(x)
	p<-ncol(x)
	mu<-matrix(mu, 1, length(mu))
	if (!is.null(lisrel$lambda.y)){
		## SEM model
		lambx<-lisrel$lambda.x
		lamby<-lisrel$lambda.y
		B<-lisrel$beta
		A<-diag(nrow(B))
		AB<-A-B
		gamma<-lisrel$gamma
		phi<-lisrel$phi
		psi<-lisrel$psi
		thdx<-lisrel$theta.x
		they<-lisrel$theta.y
		
		## computing Sigma
		ABin<-solve(AB)
		
		sigma<-lisrel$sigma0
		
		## constructing distance based on factor model notation
		lamb<-rbind(lambx, lamby%*%ABin%*%gamma) ## factor loading matrix
		lambt<-t(lamb)
		px<-nrow(thdx)
		py<-nrow(they)
		psi_new<-rbind(cbind(thdx,matrix(0,px,py) ),cbind(matrix(0,py,px),they+lamby%*%ABin%*%psi%*%t(ABin)%*%t(lamby)))
		psiin<-solve(psi_new);
		Pmat<-diag(p)-lamb%*%solve(lambt%*%psiin%*%lamb)%*%lambt%*%psiin
		Covy<-Pmat%*%sigma%*%t(Pmat)
		Av<-psiin%*%lamb
		Avc<-semdiag.orthog(Av)
		cove<-t(Avc)%*%Covy%*%Avc
		covein<-solve(cove)

		dr2<-matrix(0,n,1)

		for (i in 1:n){
	 		xi<-x[i,]-mu
			yi<-Pmat%*%t(xi) # yi is a col vector
			ei<-t(Avc)%*%yi # a col vector 
			dr2[i]<-t(ei)%*%covein%*%ei
			}
		}else{
		## for a factor model
		lamb<-lisrel$lambda.x
		phi<-lisrel$phi
		psi<-lisrel$theta.x
			
		psiin<-solve(psi)
		lambt<-t(lamb)
		sigma<-lisrel$sigma0
		I_p<-diag(p)
		pmat<-I_p-lamb%*%solve(lambt%*%psiin%*%lamb)%*%lambt%*%psiin
		covy<-pmat%*%sigma%*%t(pmat)
		Av<-psiin%*%lamb
		Avc<-semdiag.orthog(Av)
		Avct<-t(Avc)
		cove<-Avct%*%covy%*%Avc
		covein<-solve(cove)

		dr2<-matrix(0,n,1)
		for (i in 1:n){
			xi<-x[i,]-mu
			yi<-pmat%*%t(xi)
			ei<-Avct%*%yi
			dr2[i]<-t(ei)%*%covein%*%ei
			}			
		}
	dr2
	}

## Function to calculate M-distance for residuals using all factors
## x: data matrix
## mu: saturated mean from the robust estimation
## lisrel: lisrel notation model
semdiag.mdist.r1.r<-function(x, mu, lisrel){
	if(!is.matrix(x)) x<-data.matrix(x)
	n<-nrow(x)
	p<-ncol(x)
	mu<-matrix(mu, 1, length(mu))
	if (!is.null(lisrel$lambda.y)){
## SEM model
		lambx<-lisrel$lambda.x
		lamby<-lisrel$lambda.y
		B<-lisrel$beta
		A<-diag(nrow(B))
		AB<-A-B
		gamma<-lisrel$gamma
		phi<-lisrel$phi
		psi<-lisrel$psi
		thdx<-lisrel$theta.x
		they<-lisrel$theta.y
		
## computing Sigma
		ABin<-solve(AB)
		sigma<-lisrel$sigma0
		
## constructing distance based on factor model notation
		px<-nrow(thdx)
		pb<-nrow(B)
		py<-nrow(they)
		lamb<-rbind(cbind(lambx,matrix(0,px,pb)),cbind(lamby%*%ABin%*%gamma, lamby%*%ABin) ) ## factor loading matrix		
		psi_new<-rbind(cbind(thdx,matrix(0,px,py) ),cbind(matrix(0,py,px),they ))
		lambt<-t(lamb)
		psiin<-solve(psi_new);
		Pmat<-diag(p)-lamb%*%solve(lambt%*%psiin%*%lamb)%*%lambt%*%psiin
		Covy<-Pmat%*%sigma%*%t(Pmat)
		Av<-psiin%*%lamb
		Avc<-semdiag.orthog(Av)
		cove<-t(Avc)%*%Covy%*%Avc
		covein<-solve(cove)
		
		dr2<-matrix(0,n,1)
		
		for (i in 1:n){
	 		xi<-x[i,]-mu
			yi<-Pmat%*%t(xi) # yi is a col vector
			ei<-t(Avc)%*%yi # a col vector 
			dr2[i]<-t(ei)%*%covein%*%ei
		}
	}else{
## for a factor model
		lamb<-lisrel$lambda.x
		phi<-lisrel$phi
		psi<-lisrel$theta.x
		
		psiin<-solve(psi)
		lambt<-t(lamb)
		sigma<-lisrel$sigma0
		I_p<-diag(p)
		pmat<-I_p-lamb%*%solve(lambt%*%psiin%*%lamb)%*%lambt%*%psiin
		covy<-pmat%*%sigma%*%t(pmat)
		Av<-psiin%*%lamb
		Avc<-semdiag.orthog(Av)
		Avct<-t(Avc)
		cove<-Avct%*%covy%*%Avc
		covein<-solve(cove)
		
		dr2<-matrix(0,n,1)
		for (i in 1:n){
			xi<-x[i,]-mu
			yi<-pmat%*%t(xi)
			ei<-Avct%*%yi
			dr2[i]<-t(ei)%*%covein%*%ei
		}			
	}
	dr2
}

## A function to assign starting values to ram path
semdiag.start<-function(ram.path, coeff){
	n<-nrow(ram.path)
	for (i in 1:n){
		if (!is.na(ram.path[i,2])) ram.path[i,3]<-coeff[ram.path[i,2]]
		}
	ram.path	
}

#--------------------------------------------------------------------------*;
# parameter estimation by H_r(varphi) by iteratively calling sem package;
##Huber-type weight
#--------------------------------------------------------------------------*;
## theta0: starting values
## x: data
## varphi: 
semdiag.robfit.r<-function(lisrel0, x, q, varphi, ram.path, max_it=1000){
	ep<-.00001
	n<-nrow(x)
	p<-ncol(x)
	p1<-p-q 
	x<-data.matrix(x)
	#Tunning parameters in Huber weight;
    prob<-1-varphi
	chip1<-qchisq(prob, p1)
	ck<-sqrt(chip1)
	cbeta<-( p1*pchisq(chip1,p1+2)+ chip1*(1-prob) )/p1

	mu0<-colMeans(x)
	theta0<-lisrel0$coeff
	sigma0<-lisrel0$sigma0

	Md_r<-matrix(0,n,3)	

	n_it<-1
	dt<-1
	while (dt>ep && n_it <= max_it){
		dr2<-semdiag.mdist.r.r(x,mu0,lisrel0)
		sumw1<-0.0
		mu<-rep(0,p)
		sigma<-matrix(0,p,p)

		for (i in 1:n) {
			xii<-x[i,] 
  			xi0<-xii-mu0 
  			di2<-dr2[i]
  			di<-sqrt(di2)
  		
			#Huber weight functions
  			if (di<=ck) {
     			wi1<-1.0
     			wi2<-1.0/cbeta
			}else {
     			wi1<-ck/di
     			wi2<-wi1^2/cbeta
 			}

			Md_r[i,1]<-i
  			Md_r[i,2]<-di
  			Md_r[i,3]<-wi2
  			sumw1<-sumw1+wi1
  			mu<-mu+c(wi1)*xii
  			sigma<-sigma+c(wi2)*xi0%*%t(xi0)  		
		} # end of loop n 

		mu1<-mu/sumw1
		sigma1<-sigma/n
		
		rownames(sigma1)<-colnames(sigma1)
		
		ram.path<-semdiag.start(ram.path, theta0)
	
		model1<-sem(ram.path, sigma1, N=n, maxiter=1, analytic.gradient=TRUE, warn=FALSE)
		##model1<-semdiag.gradient(ram.path, sigma1, theta0)
	    ## dtheta<-model1$gradient
	
		## estimated model parameters
		theta1<-model1$coeff
	
		dt<-max(c(max(abs(mu1-mu0)), max(abs(theta1-theta0))))
		#print(dt)	
		n_it<-n_it+1
		mu0<-mu1
		
		lisrel0<-semdiag.lisrel(model1)
		theta0<-theta1			
  } # end of loop k              
 list(Md_r=Md_r,mu1=mu1, results=model1, n_it=n_it)                  
}	
	
#--------------------------------------------------------------------------*;
# parameter estimation by H_r(varphi) using iteratively reweighted least squares;
##Huber-type weight
#--------------------------------------------------------------------------*;
## theta0: starting values
## x: data
## varphi: 
semdiag.robfit1.r<-function(lisrel0, x, q, varphi, ram.path, max_it=1000){
	ep<-.0001
	n<-nrow(x)
	p<-ncol(x)
	p1<-p-q 
	x<-data.matrix(x)
#Tunning parameters in Huber weight;
    prob<-1-varphi
	chip1<-qchisq(prob, p1)
	ck<-sqrt(chip1)
	cbeta<-( p1*pchisq(chip1,p1+2)+ chip1*(1-prob) )/p1
	
	mu0<-colMeans(x)
	theta0<-lisrel0$coeff
	sigma0<-lisrel0$sigma0
	
	Md_r<-matrix(0,n,3)	
	
	n_it<-1
	dt<-1
	while (dt>ep && n_it <= max_it){
		dr2<-semdiag.mdist.r1.r(x,mu0,lisrel0)
		sumw1<-0.0
		mu<-rep(0,p)
		sigma<-matrix(0,p,p)
		
		for (i in 1:n) {
			xii<-x[i,] 
  			xi0<-xii-mu0 
  			di2<-dr2[i]
  			di<-sqrt(di2)
			
#Huber weight functions
  			if (di<=ck) {
     			wi1<-1.0
     			wi2<-1.0/cbeta
			}else {
     			wi1<-ck/di
     			wi2<-wi1^2/cbeta
 			}
			
			Md_r[i,1]<-i
  			Md_r[i,2]<-di
  			Md_r[i,3]<-wi2
  			sumw1<-sumw1+wi1
  			mu<-mu+c(wi1)*xii
  			sigma<-sigma+c(wi2)*xi0%*%t(xi0)  		
		} # end of loop n 
		
		mu1<-mu/sumw1
		sigma1<-sigma/n
		rownames(sigma1)<-colnames(sigma1)
		ram.path<-semdiag.start(ram.path, theta0)
		model1<-sem(ram.path, sigma1, N=n, maxiter=1, analytic.gradient=TRUE, warn=FALSE)
		
## estimated model parameters
		theta1<-model1$coeff
		
		dt<-max(c(max(abs(mu1-mu0)), max(abs(theta1-theta0))))	
		n_it<-n_it+1
		mu0<-mu1
		lisrel0<-semdiag.lisrel(model1)
		theta0<-theta1
	} # end of loop k              
	list(Md_r=Md_r,mu1=mu1,results=model1, n_it=n_it)                  
}	

## The main function for analysis

semdiag<-function(x,  EQSmodel, varphi=0.1, EQSdata='data.txt', D='E', delete=integer(0), max_it=1000, EQSprog='C:/Progra~1/EQS61/WINEQS',serial="111111 222222 333333", ram.path, software='EQS'){

	if (!is.matrix(x)) x<-data.matrix(x)
	if (D=='E'){
		  model=0
		}else{
		  model=1	
		}
	

	## Check whether cases need to be deleted

	N<-nrow(x)

	index<-1:N

	rownames(x)<-index

## check for errors
if (software=='EQS' && missing(EQSmodel)) stop ("No EQS input file is supplied!")
if (software=='sem' && missing(ram.path)) stop ("No ram path is specified for the sem package!")

	if (software=='EQS'){	

	if (length(delete) > 0){

		index<-index[-delete]

		x<-x[index, ]

		eqsinput<-readLines(EQSmodel, warn=FALSE)



     ## find the line with cases

     lno<-grep('cases', eqsinput, ignore.case=TRUE)

     if (length(lno)==0 | length(lno)>1) stop("The number of cases is not supplied in your EQS input file!")

     N.new<-N-length(delete)

     temp.str<-eqsinput[lno]

     temp.str.new<-sub(N, N.new, temp.str)

     eqsinput[lno]<-temp.str.new

     

     if (length(grep('/', EQSmodel))>0){

     	filedir.split <- strsplit(EQSmodel, "/")[[1]]

     }else{

     	filedir.split <- strsplit(EQSmodel, "\\\\")[[1]]

     }

     n <- length(filedir.split)

     eqsname <- strsplit(filedir.split[n], "\\.")[[1]][1]

     

     lno<-grep(eqsname, eqsinput, ignore.case=TRUE)

     temp.str<-eqsinput[lno]

     temp.str.new<-sub(eqsname, 'deqs', temp.str)

     eqsinput[lno]<-temp.str.new

     

     ## write the new input file

     eqsinput[lno]<-temp.str.new

     writeLines(text=eqsinput, con='deqs.eqs')

		EQSmodel<-'deqs.eqs'		

		}

	

	## run an ML

	ml.cov<-cov(x)

	write.table(ml.cov, EQSdata, row.names=FALSE, col.names=FALSE)

	ml.temp<-semdiag.run.eqs(EQSprog, EQSmodel, serial=serial)

	nml.res<-semdiag.eqs(ml.temp)

	## copy the results to a file called nml.out
	nml.out<-unlist(strsplit(EQSmodel, '.', fixed=TRUE))[1]
	nml.out<-paste(nml.out, '.out', sep='')
	file.rename(nml.out, 'nml.out')

	p<-ncol(x)

## robust estimation for mu and sigma

	robmusig<-semdiag.musig(x, varphi)

## write data into data.txt and run EQS

	write.table(robmusig$sigma, EQSdata, row.names=FALSE, col.names=FALSE)
	model0<-semdiag.run.eqs(EQSprog, EQSmodel, serial=serial)
	tsr.res<-semdiag.eqs(model0)
	file.rename(nml.out, 'tsr.out')

## lisrel model 	
	model.lisrel<-semdiag.eqs.lisrel(model0)
	q<-nrow(model.lisrel$phi)
	if (model==1) q<-q+nrow(model.lisrel$psi)
## Calculate M-distance for f
	if (model==0){
		d_f<-semdiag.mdist.f(x, robmusig$mu, model.lisrel)
	}else{
		d_f<-semdiag.mdist.f1(x, robmusig$mu, model.lisrel)
	}
	
## iterative run EQS to calculate M-distance for residuals

	if (model==0){
		d_r<-semdiag.robfit(model.lisrel, x, q, varphi, EQSmodel, EQSdata, max_it=max_it, EQSprog=EQSprog, serial=serial)
	}else{
		d_r<-semdiag.robfit1(model.lisrel, x, q, varphi, EQSmodel, EQSdata, max_it=max_it, EQSprog=EQSprog, serial=serial)
	}
	dr.res<-semdiag.eqs(d_r$eqs)
	file.rename('eqsonce.out', 'dr.out')
	list(d_f=d_f, d_r=d_r$Md_r[,1:2], mu=robmusig$mu, p=p, q=q, res=list(nml=nml.res, tsr=tsr.res, dr=dr.res), eqs=list(nml=ml.temp, tsr=model0, dr=d_r$eqs), x=x, software=software)
	}else{
		## start the use of R
		if (length(delete) > 0){
		index<-index[-delete]
		x<-x[index, ]	
		}
	
	## run an ML
	ml.cov<-cov(x)	
	ml.temp<-sem(ram.path, ml.cov, N=N, warn=FALSE)
	nml.res<-summary(ml.temp)
		
	p<-ncol(x)
## robust estimation for mu and sigma
	robmusig<-semdiag.musig(x, varphi)
	tsr.cov<-robmusig$sigma
	rownames(tsr.cov)<-colnames(tsr.cov)
	model0<-sem(ram.path, tsr.cov, N=N)
	tsr.res<-summary(model0)
## lisrel model 	
	model.lisrel<-semdiag.lisrel(model0)
	q<-nrow(model.lisrel$phi)
	if (model==1) q<-q+nrow(model.lisrel$psi)
## Calculate M-distance for f
	if (model==0){
		d_f<-semdiag.mdist.f.r(x, robmusig$mu, model.lisrel)
	}else{
		d_f<-semdiag.mdist.f1.r(x, robmusig$mu, model.lisrel)
	}
	
## iterative run EQS to calculate M-distance for residuals
	if (model==0){
		d_r<-semdiag.robfit.r(model.lisrel, x, q, varphi, ram.path, max_it=max_it)
	}else{
		d_r<-semdiag.robfit1.r(model.lisrel, x, q, varphi, ram.path, max_it=max_it)
	}
	dr.res<-summary(d_r$results)
	list(d_f=d_f, d_r=d_r$Md_r[,1:2], mu=d_r$mu, p=p, q=q, res=list(nml=nml.res, tsr=tsr.res, dr=dr.res), x=x, software=software)
		
		}
}



semdiag.summary<-function(d, alpha=.01, digits=2){
	prob<-1-alpha
	p<-d$p
	q<-d$q
	critq<-sqrt( qchisq(prob,q) )
	critp1<-sqrt( qchisq(prob,p-q) ) 
	dfi<-d$d_f
	dri<-d$d_r[,2]
	index<-d$d_r[,1]
	x<-d$x

	l.o.index<-index[(dri>critp1) & (dfi>critq)]
	o.index<-index[(dri>critp1) & (dfi<=critq)]
	l.index<-index[(dri<=critp1) & (dfi>critq)]
	
	if (length(l.o.index)>=1){
		cat("Leverage observations and outliers =",  rownames(x)[l.o.index], "\n")
	}else{
		cat("No case is both outlier and leverage observation. \n")
	}

	if (length(l.index)>=1){
		cat("Leverage observations not outliers =",  rownames(x)[l.index], "\n")
	}else{
		cat("No case is just leverage observation. \n")
	}

	if (length(o.index)>=1){
		cat("Outliers not leverage observations =",  rownames(x)[o.index], "\n")
	}else{
		cat("No case is just outliers. \n")
	}

	

## Model information and parameter estimates

## Model fit information
if (d$software=='EQS'){
	fit<-rbind(cbind(d$res$nml$chi, d$res$tsr$chi, d$res$dr$chi),cbind(d$res$nml$df, d$res$tsr$df, d$res$dr$df), cbind(d$res$nml$pval, d$res$tsr$pval, d$res$dr$pval))
	colnames(fit)<-c('NML','TSR', 'DR')
	rownames(fit)<-c('Statistics','df','p-value')
## parameter estimates
	estimates<-cbind(d$res$nml$par, d$res$tsr$par, d$res$dr$par)
	estimates<-round(estimates, digits=digits)
	estimates<-format(estimates, digits=digits)
	row0<-c('','NML','','','TSR','','', 'DR','')
	row1<-c('Est.','S.E.','z','Est.','S.E.','z','Est.','S.E.','z')
	estimates<-rbind(row1, estimates)
	rownames(estimates)[1]<-c('Label')
	colnames(estimates)<-row0<-c('','NML','','','TSR','','', 'DR','')
	cat('\nModel fit comparison\n')
	print(fit)
	cat('\nParameter estimates\n')
	print(estimates, quote=FALSE, right=TRUE)
	
	cat('\nNote.\n')
	cat('  NML = Normal ML\n')
	cat('  TSR = Two-stage robust method\n')
	cat('  DR  = Direct robust method\n')
	cat('  Est.= Parameter estimates\n')
	cat('  S.E.= Standard error\n')
	cat('  z   = Z-score\n')
	invisible(list(fit=fit, estimates=estimates))
	}else{
		## Model fit information

	chisq1<-d$res$nml$chisq
	chisq2<-d$res$tsr$chisq
	chisq3<-d$res$dr$chisq
	df1<-d$res$nml$df
	df2<-d$res$tsr$df
	df3<-d$res$dr$df
	
	pvalue1<-1-pchisq(chisq1, df1)
	pvalue2<-1-pchisq(chisq2, df2)
	pvalue3<-1-pchisq(chisq3, df3)
	

	fit<-rbind(cbind(chisq1,chisq2,chisq3),cbind(df1,df2,df3), cbind(pvalue1,pvalue2,pvalue3))
	fit<-round(fit, 4)
	colnames(fit)<-c('NML','TSR', 'DR')
	rownames(fit)<-c('Statistics','df','p-value')
## parameter estimates
	estimates<-cbind(d$res$nml$coeff[,1:3], d$res$tsr$coeff[,1:3], d$res$dr$coeff[,1:3])
	estimates<-round(estimates, digits=digits)
	estimates<-format(estimates, digits=digits)
	estimates<-as.matrix(estimates)
	row0<-c('','NML','','','TSR','','', 'DR','')
	row1<-c('Est.','S.E.','z','Est.','S.E.','z','Est.','S.E.','z')
	estimates<-rbind(row1, estimates)
	rownames(estimates)[1]<-c('Label')
	colnames(estimates)<-row0<-c('','NML','','','TSR','','', 'DR','')
	cat('\nModel fit comparison\n')
	print(fit)
	cat('\nParameter estimates\n')
	print(estimates, quote=FALSE, right=TRUE)
	
	cat('\nNote.\n')
	cat('  NML = Normal ML\n')
	cat('  TSR = Two-stage robust method\n')
	cat('  DR  = Direct robust method\n')
	cat('  Est.= Parameter estimates\n')
	cat('  S.E.= Standard error\n')
	cat('  z   = Z-score\n')
	invisible(list(fit=fit, estimates=estimates))
		}
}



semdiag.plot<-function(d, alpha=.01, label=0, cex=1){

	prob<-1-alpha

	p<-d$p

	q<-d$q

	critq<-sqrt( qchisq(prob,q) )

	critp1<-sqrt( qchisq(prob,p-q) ) 

	dfi<-d$d_f

	dri<-d$d_r[,2]

	index<-d$d_r[,1]

	x<-d$x

	

	l.o.index<-index[(dri>critp1) & (dfi>critq)]

	o.index<-index[(dri>critp1) & (dfi<=critq)]

	l.index<-index[(dri<=critp1) & (dfi>critq)]

	

	if (length(l.o.index)>=1){

			cat("Leverage observations and outliers =",  rownames(x)[l.o.index], "\n")

		}else{

			cat("No case is both outlier and leverage observation. \n")

		}

	if (length(l.index)>=1){

			cat("Leverage observations not outliers =",  rownames(x)[l.index], "\n")

		}else{

			cat("No case is just leverage observation. \n")

		}

	if (length(o.index)>=1){

			cat("Outliers not leverage observations =",  rownames(x)[o.index], "\n")

		}else{

			cat("No case is just outliers. \n")

		}



		## plot of distance

		par(mfrow=c(1,2))

		plot(dfi, dri, xlab='M-distance d_f (leverage obs.)', ylab='M-distance d_r (potential outliers)', main='A', pch=20)

		abline(v=critq, h=critp1, lty=2)

				

		## label leverage points and outliers

		textlabel<-unique(c(o.index,l.index,l.o.index))

		if (label == 1) identify(dfi, dri, rownames(x), cex=cex)

		if (label == 0) if (length(textlabel>0)) text(dfi[textlabel], dri[textlabel], rownames(x)[textlabel], cex=cex)

		if (label==2) {
		    if (length(textlabel>0)) text(dfi[textlabel], dri[textlabel], rownames(x)[textlabel], cex=cex)
			identify(dfi, dri, rownames(x), cex=cex)
			}
		

		## Q-Q plot

		n<-length(dfi)

		Qchi<-sqrt(qchisq((1:n)/(n+1), p-q))

		dr<-d$d_r

		norder<-order(dr[,2])

		dr.o<-dr[norder, ]

		data<-cbind(dr.o, Qchi)

		plot(data[,3], data[,2], xlab='Quantiles of Chi-distribution', ylab='M-distance d_r', main='B', pch=20)

		textlabel<-unique(c(o.index,l.o.index))

		if (label ==1) identify(data[,3], data[,2], rownames(x)[norder], cex=cex)

		if (label ==0){ 
		    m<-length(textlabel)
			if (m>0){
				for (j in 1:m){
					s<-data[,1]==textlabel[j]
					text(data[s,3], data[s,2], rownames(x)[textlabel[j]], cex=cex)
				}				
			}
		}
		if (label==2) {
			m<-length(textlabel)
			if (m>0){
				for (j in 1:m){
					s<-data[,1]==textlabel[j]
					text(data[s,3], data[s,2], rownames(x)[textlabel[j]], cex=cex)
				}				
			}	
			identify(data[,3], data[,2], rownames(x)[norder], cex=cex)
			}

	par(mfrow=c(1,1))

	}

semdiag.cpp<-function(d, cases){
    x<-d$x
    mu<-d$mu
	p<-ncol(x)
	n<-nrow(x)
	xmean<-matrix(mu, n, p, byrow=TRUE)
	xcenter<-x-xmean	
	ylim<-round(max(abs(max(xcenter)), abs(min(xcenter))),1)
	plot(1:p, rep(0,p), ylim=c(-ylim, ylim), type='n', axes=FALSE, xlab='Order of variables', ylab=expression(v[ij]-hat(mu)[j]))
	abline(h=0)
	axis(2, at=c(-ylim, 0, ylim))
	axis(1, at=1:p, labels=colnames(x))
	for (i in cases){
		lines(1:p, xcenter[i, ])
		text(1:p, xcenter[i, ], i)
	}
}	

