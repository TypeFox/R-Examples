setClass("matrixCorr", representation(K="numeric",Run1="numeric",Run2="numeric",CorrMatrix="matrix",Pvalues="matrix"))
matrixCorr <- function(K, Run1, Run2, CorrMatrix, Pvalues = matrix(NA)) {
	res <- new("matrixCorr")
	res@K <- K
	res@Run1 <- Run1
	res@Run2 <- Run2
	res@CorrMatrix <- CorrMatrix
	res@Pvalues <- Pvalues
	return(res)
}
setClass("rowncolMatrix", representation(K="numeric", filtermatrix = "table"))
rowncolMatrix <- function(K, filtermatrix) {
	res <- new("rowncolMatrix")
	res@K <- K
	res@filtermatrix <- filtermatrix
	return(res)
}
setClass("QmatrixFilt", representation(rowncol="list",avmaxcorr="table",rawcorr="list"))
QmatrixFilt <- function(rowncol = list(""), avmaxcorr = as.table(matrix(NA)), rawcorr = list("")) {
	res <- new("QmatrixFilt")
	res@rowncol <- rowncol
	res@avmaxcorr <- avmaxcorr
	res@rawcorr <- rawcorr
	return(res)
}
read.struct <- function(filepath = "./", instruct = FALSE) {
	options(warn = -9)
	files <- list.files(path = filepath, pattern = "_f$")
	maxK <- 1
	for (i in files) {
		lin <- paste(filepath, i, sep = "")
		x <- readLines(lin)
		if (instruct == TRUE) {
			K <- x[c(grep("Population number assumed=", as.character(x), perl = TRUE))]
			K <- as.numeric(substring(K, 31, nchar(K)))}
		else {
			K <- x[c(grep("populations assumed", as.character(x), perl = TRUE))]
			K <- as.numeric(substring(K, 4, nchar(K)-20))}
		if (K > maxK) {maxK <- K}
	}
	data <- "NULL"
	for (i in files) {
		lin <- paste(filepath, i, sep = "")
		x <- readLines(lin)
		if (instruct == TRUE) {
			K <- x[c(grep("Population number assumed=", as.character(x), perl = TRUE))]
			lnPD <- x[c(grep("Posterior Mean = ", as.character(x), perl = TRUE))]
			K <- as.numeric(substring(K, 31, nchar(K)))}
		else {
			K <- x[c(grep("populations assumed", as.character(x), perl = TRUE))]
			lnPD <- x[c(grep("Estimated Ln Prob of Data", as.character(x), perl = TRUE))]
			K <- as.numeric(substring(K, 4, nchar(K)-20))}
		Fsts <- c()
		for (j in 1:K) {
			line <- paste("Mean value of Fst_", as.character(j), " ", sep = "")
			Fst <- x[c(grep(line, as.character(x), perl = TRUE))]
			Fst <- as.numeric(substring(Fst, nchar(Fst)-6, nchar(Fst)))
			Fsts <- c(Fsts, Fst)
		}
		if (length(Fsts) < maxK) {
			for (l in 1:(maxK-length(Fsts))) {
				Fsts <- c(Fsts, NA)
			}
		}
		if (instruct == TRUE) {
			lnPD <-as.numeric(substring(lnPD, 22, 31))}
		else {
			lnPD <- as.numeric(substring(lnPD, 31, 37))}
		dat <- matrix(c(K, lnPD, Fsts), ncol = (maxK + 2),byrow = TRUE)
		if (is(data, "character")) {data <- dat} else {data <- rbind(data, dat)}
	}
	nam <- c("K", "lnPD")
	for (i in 1:maxK) {
		nam <- c(nam, paste("Fst_",as.character(i),sep= ""))
	}
	colnames(data) <- nam
	data <- data.frame(data)
	data <- data[order(data$K),]
	return(data)
}
corr.Qmatrix <- function(filepath  = "./", instruct = FALSE, rowncol = TRUE, avmax = TRUE, pvalue = FALSE, raw = TRUE, r = 0.99, p = 0.05) {
	options(warn = -9)
	files <- list.files(path = filepath, pattern = "_f$")
	maxK <- 1
	for (i in files) {
		lin <- paste(filepath, i, sep = "")
		x <- readLines(lin)
		if (instruct == TRUE) {
			K <- x[c(grep("Population number assumed=", as.character(x), perl = TRUE))]
			K <- as.numeric(substring(K, 31, nchar(K)))}
		else {
			K <- x[c(grep("populations assumed", as.character(x), perl = TRUE))]
			K <- as.numeric(substring(K, 4, nchar(K)-20))}
		if (K > maxK) {maxK <- K}
	}
	Qmatrices <- "NULL"
	NoColK <- c()
	for (j in 1:maxK) {
		Kfiles <- c()
		for (i in files) {
			lin <- paste(filepath, i, sep = "")
			x <- readLines(lin)
			if (instruct == TRUE) {
				K <- x[c(grep("Population number assumed=", as.character(x), perl = TRUE))]
				K <- as.numeric(substring(K, 31, nchar(K)))}
			else {
				K <- x[c(grep("populations assumed", as.character(x), perl = TRUE))]
				K <- as.numeric(substring(K, 4, nchar(K)-20))}
			if (j == K) {Kfiles <- c(Kfiles, lin)}
		}
		NoColK[j] <- length(Kfiles)
		for (i in Kfiles) {
			Qmatrix <- c()
			x <- readLines(i)
			Qmat <- x[c(grep(" : ", as.character(x), perl = TRUE))]
			Qmat <- Qmat[c(grep("cluster", as.character(Qmat), invert = TRUE, perl = TRUE))]
			Qmat <- Qmat[c(grep("Cluster", as.character(Qmat), invert = TRUE, perl = TRUE))]
			Qmat <- Qmat[c(grep("Locus", as.character(Qmat), invert = TRUE, perl = TRUE))]
			Qmat <- gsub("inf", "0.000", Qmat, perl = TRUE)
			for (l in 1:j) {
				if (instruct == TRUE && l == 1) {beg <-nchar(Qmat)-4} else {beg <- nchar(Qmat)-5}
				Qmat2 <- as.numeric(substring(Qmat, beg, beg+4))
				Qmat <- substring(Qmat, 1, beg - 1)
				Qmatrix <- c(Qmat2, Qmatrix)
			}
			Qmatrix <- matrix(Qmatrix, nrow = j, byrow = TRUE)
			Qmatrix <- t(Qmatrix)
			if (is(Qmatrices, "character")) {Qmatrices <- Qmatrix} else {Qmatrices <- cbind(Qmatrices, Qmatrix)}
		}
	}
	Qmatrices <- data.frame(Qmatrices)
	rawmatrix <- c()
	avmaxcorr <- c()
	rowcolALL <- c()
	rowcolcorr <- c()
	end <- 0
	for (i in 1:maxK) {	
		start <- end + 1
		end <- end + (NoColK[i] * i)
		if (NoColK[i] > 1) {
			corr <- cor(Qmatrices[,start:end], Qmatrices[,start:end])
			if (pvalue == TRUE) {
				sub <- Qmatrices[,start:end]
                        len  <- i * NoColK[i]
                        pvar <- c()
                        for (j in 1:len) {
                                xt <- c(sub[,j])
                                for (k in 1:len) {
                                        yt <- c(sub[,k])
                                        ptmp <- cor.test(xt, yt, method = "pearson")
                                        pvar <- c(pvar, ptmp$p.value)
                                }
                        }
                        pvar <- matrix(pvar, ncol = len, byrow = TRUE)
			}
			maxcorr <- c()
			rowcolcorr <- c()
			for (j in 1:NoColK[i]) {
				ya <- (j-1)*i + 1
				yb <- j*i
				min <- j+1
				if (min < NoColK[i] + 1) {
					for (k in min:NoColK[i]) {
						xa <- (k-1)*i+1
						xb <- k*i
						rmat <- matrix(corr[ya:yb,xa:xb], ncol = i, byrow = TRUE)
						if (pvalue == TRUE) {
							pv <- matrix(pvar[ya:yb,xa:xb], ncol = i, byrow = TRUE)
							mat <- matrixCorr(K = i, Run1 = j, Run2 = k, CorrMatrix = rmat, Pvalues = pv)
							if (i > 1) {
								for (l in 1:i) {
									for (m in 1:i) {
                                        if (is.na(pv[l,m])) {rmat[l,m] <- NA} else {
                                            if (pv[l,m] > p) {rmat[l,m] <- 0}
                                        }
									}
								}
						}
						} else {mat <- matrixCorr(K = i, Run1 = j, Run2 = k, CorrMatrix = rmat)}
						rawmatrix <- c(rawmatrix, mat)
                        mx <- max(abs(rmat))
                        maxcorr <- c(maxcorr, mx)
						if (rowncol == TRUE) {
							if (i == 1) {rowcolcorr <- c(rowcolcorr, "?")} else {
								sig = TRUE
                                corflag = FALSE
								for (m in 1:i) {
                                    if (NA %in% rmat[m,]) {corflag <- TRUE} else {
                                        if (max(abs(rmat[m,])) < r) {sig <- FALSE}
                                    }
								}
								for (m in 1:i) {
                                    if (NA %in% rmat[,m]) {corflag <- TRUE} else {
                                        if (max(abs(rmat[,m])) < r) {sig <- FALSE}
                                    }
								}
								if (corflag == TRUE) {rowcolcorr <-c(rowcolcorr, "?")} else {
                                    if (sig == TRUE) {rowcolcorr <- c(rowcolcorr, "Y")} else {rowcolcorr <- c(rowcolcorr, "N")}
                                }
							}
						}
					}
				}
			}
			if (rowncol == TRUE) {
				if (NoColK[i] > 2) {
					for (z in 1:(NoColK[i]-2)) {
						ca <- rowcolcorr[1:(NoColK[i]*z-z)]
						for (y in 1:z) {ca <- c(ca, NA)}
						rowcolcorr <- c(ca, rowcolcorr[(NoColK[i]*z-z + 1):length(rowcolcorr)])
					}
				}
				rowcolcorr <- as.table(matrix(rowcolcorr, ncol = (NoColK[i] - 1), byrow = TRUE))
				rownames(rowcolcorr) <- c(2:NoColK[i])
				colnames(rowcolcorr) <- c(1:(NoColK[i]-1))
				rowcolcorr <- rowncolMatrix(K = i, filtermatrix = rowcolcorr)
				rowcolALL <- c(rowcolALL, rowcolcorr)
			}
			avmaxcorr <- c(avmaxcorr, mean(maxcorr))
		}else {
			avmaxcorr <- c(avmaxcorr, NA)
			rowcolcorr <- c(rowcolcorr, table(NA))
			}
	}
	results <- QmatrixFilt()
	if (avmax == TRUE) {
		avmaxfilter <- c()
		for (i in 1:maxK){
			if (is.na(avmaxcorr[i]) == TRUE) 
				{avmaxfilter <- c(avmaxfilter, "?")}
			else	
				{if (avmaxcorr[i] < r) {avmaxfilter <- c(avmaxfilter, "N")} else {avmaxfilter <- c(avmaxfilter, "Y")}}			
		}	
		avmaxcorr <- rbind(c(1:maxK),avmaxcorr,avmaxfilter)
		avmaxcorr <- as.table(matrix(avmaxcorr, nrow = 3, byrow = FALSE))
		rownames(avmaxcorr) <- c("K", "Average Max Corr", "Significantly Correlated")
		results@avmaxcorr <- avmaxcorr
	}
	if (raw == TRUE) {results@rawcorr <- rawmatrix}
	if (rowncol == TRUE) {results@rowncol <- rowcolALL}
	return(results)
}
summarise.lnPD <- function(input) {
	lnPD <- "NULL"
	for (i in min(input$K):max(input$K)) {
		data <- c()
		for (j in 1:nrow(input)){
			if (input$K[j]==i) {data <- c(data, input$lnPD[j])}
		}
		dat <- matrix(c(i, mean(data), sd(data)), ncol = 3, byrow = TRUE)
		if (is(lnPD, "character")) {lnPD <- dat} else {lnPD <- rbind(lnPD, dat)}
	}
	colnames(lnPD) <- c("K", "Mean", "StDev")
	lnPD <- data.frame(lnPD)
	return(lnPD)
}
summarize.lnPD <- function(input) {return(summarise.lnPD(input))}
summarise.Fst <- function(input, stdevopt = 1) {
	if (stdevopt == 3) {
		cat("\n", "Enter raw data file path", "\n")
		filepath <- scan(n=1, what = "character", quiet = TRUE)
	}
	Fst <- "NULL"
	sdmeans <- c()
	meansds <- c()
	sdssds <- c()
	fmeans <- c()
	fsds <- c()
	for (i in min(input$K):max(input$K)) {
		data <- c()
		for (j in 1:nrow(input)){			
			if (input$K[j]==i) {
				for (k in 1:i) {data <- c(data, input[,k+2][j])}	
			}
		}
		tmp <- matrix(data, ncol = i, byrow = TRUE)
		means <- c()
		sds <- c()
		if (stdevopt == 2) {
			for (kc in 1:nrow(tmp)){
				tmp2 <- t(tmp[kc,1:ncol(tmp)])
				tmp2 <- t(tmp2[order(tmp2[1,])])
				tmp[kc, 1:ncol(tmp2)] <- tmp2
			}
		}
		if (stdevopt == 3 && i > 1) {
			options(warn = -9)
			files <- list.files(path = filepath, pattern = "_f$")
			Qref <- "NULL"
			Qmatrices <- "NULL"
			Kfiles <- c()
			for (fil in files) {
				lin <- paste(filepath, fil, sep = "")
				x <- readLines(lin)
				K <- x[c(grep("populations assumed", as.character(x), perl = TRUE))]
				K <- as.numeric(substring(K, 4, nchar(K)-20))
				if (i == K) {Kfiles <- c(Kfiles, lin)}
			}
			for (fil in Kfiles) {
				Qmatrix <- c()
				x <- readLines(fil)
				Qmat <- x[c(grep(" : ", as.character(x), perl = TRUE))]
				Qmat <- Qmat[c(grep("cluster", as.character(Qmat), invert = TRUE, perl = TRUE))]
				Qmat <- Qmat[c(grep("Locus", as.character(Qmat), invert = TRUE, perl = TRUE))]
				for (l in 1:i) {
					beg <- nchar(Qmat)-5
					Qmat2 <- as.numeric(substring(Qmat, beg, beg+4))
					Qmat <- substring(Qmat, 1, beg-1)
					Qmatrix <- c(Qmat2, Qmatrix)	
				}
				
				Qmatrix <- matrix(Qmatrix, nrow = i, byrow = TRUE)
				Qmatrix <- t(Qmatrix)
				if (is(Qmatrices, "character") && is(Qref, "character")) {Qref <- Qmatrix}
					else if (is(Qmatrices, "character")) {Qmatrices <- Qmatrix} else {Qmatrices <- cbind(Qmatrices, Qmatrix)}
			}
			xb <- 0
			for (jy in 2:length(Kfiles)) {
				xa <- xb +1
				xb <- xb + i
				tmp2 <- c()
				rmat <- cor(Qref, Qmatrices[,xa:xb])
				for (ord in 1:i){
					rv <- -1.0
					sp <- 1
					for (a in 1:i) {
						if (rmat[ord,a] > rv) {
							rv <- rmat[ord,a]
							sp <- a}
					}
					tmp2 <- c(tmp2, tmp[jy,sp])
					rmat[,sp] <- -2.0
					}		
				tmp[jy,] <- tmp2
			}
		}
		for (j in 1:i) {
			means <- c(means, mean(tmp[,j]))
			sds <- c(sds, sd(tmp[,j]))
		}
		sdmeans <- c(sdmeans, sd(means))
		meansds <- c(meansds, mean(sds))
		sdssds <- c(sdssds, sd(sds))
		fmeans <- c(fmeans, means)
		fsds <- c(fsds, sds)
		if (i < max(input$K)) {
			for (j in 1:(max(input$K)-i)) {
				fmeans <- c(fmeans, NA)
				fsds <- c(fsds, NA)
			}
		}
		dat <- matrix(c(i, mean(data), sd(data)), ncol = 3, byrow = TRUE)
		if (is(Fst, "character")) {Fst <- dat} else {Fst <- rbind(Fst, dat)}
	}
	sdmeans <- matrix(sdmeans, ncol = 1, byrow = TRUE)
	meansds <- matrix(meansds, ncol = 1, byrow = TRUE)
	sdssds <- matrix(sdssds, ncol = 1, byrow = TRUE)
	fmeans <- matrix(fmeans, ncol = max(input$K), byrow = TRUE)
	fsds <- matrix(fsds, ncol = max(input$K), byrow = TRUE)
	Fst <- cbind(Fst, sdmeans, meansds, sdssds, fmeans, fsds)
	nam <- c()
	for (i in 1:max(input$K)) {
		nam <- c(nam, paste("Mean_Fst_",as.character(i),sep= ""))
	}
	for (i in 1:max(input$K)) {
		nam <- c(nam, paste("StDev_Fst_",as.character(i),sep= ""))
	}
	colnames(Fst) <- c("K", "Overall_Mean", "Overall_StDev", "StDev_of_Means", "Mean_StDev", "StDev_of_StDev", nam)
	Fst <- data.frame(Fst)
	return(Fst)
}
summarize.Fst <- function(input, stdevopt = 1) {return(summarise.Fst(input, stdevopt))}
calc.delta <- function(input, Fst = FALSE) {
	flag <- FALSE
	j <- min(input$K)
	for (i in min(input$K):(max(input$K)-1)) {
		if (i+1 != input$K[i-j+2]) {flag <- TRUE}
	}
	if (flag == TRUE) {print("Discontinuous data. Cannot process.")} else {
		LprK <- c(NA)
		LdprK <- c()
		deltaK <- c()
		for (i in 2:nrow(input)) {
			if (Fst ==FALSE) {lpr <- input$Mean[i] - input$Mean[i-1]} else {lpr <-input$Overall_Mean[i] - input$Overall_Mean[i-1]}
			LprK <- c(LprK, lpr)
		}
		for (i in 2:nrow(input)-1) {
			ldpr <- abs(LprK[i+1] - LprK[i])
			LdprK <- c(LdprK, ldpr)
			if (Fst ==FALSE) {dK <- ldpr/input$StDev[i]} else {dK <- ldpr/input$Overall_StDev[i]}
			deltaK <- c(deltaK, dK)
		}
		LdprK <- c(LdprK, NA)
		deltaK <- c(deltaK, NA)
		output <- matrix(c(input$K,LprK,LdprK,deltaK),ncol=nrow(input),byrow =TRUE)
		output <- t(output)
		if (Fst == FALSE) {colnames(output) <- c("K","LprK","LdprK","deltaK")} else {colnames(output) <- c("K","FprK","FdprK","deltaFst")}
		output <- data.frame(output)
		return(output)
	}
}
