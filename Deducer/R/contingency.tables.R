



##Corrected the formulas from M.J. Smithson
##http://psychology.anu.edu.au/people/smithson/details/CIstuff/CI.html
chi.noncentral.conf<-function (chival, df, conf, prec = 1e-05) 
{
	result <- NA
	ulim <- 1 - (1 - conf)/2
	lc <- c(0.001, chival/2, chival)
	while (pchisq(chival, df, lc[1]) < ulim) {
		if (pchisq(chival, df) < ulim) {
			result <- (c(0, pchisq(chival, df)))
			break
		}
		lc <- c(lc[1]/4, lc[1], lc[3])
	}
	diff <- 1
	if (all(is.na(result))) {
		while (diff > prec) {
			if (pchisq(chival, df, lc[2]) < ulim) 
				lc <- c(lc[1], (lc[1] + lc[2])/2, lc[2])
			else lc <- c(lc[2], (lc[2] + lc[3])/2, lc[3])
			diff <- abs(pchisq(chival, df, lc[2]) - ulim)
			ucdf <- pchisq(chival, df, lc[2])
		}
		result <- c(lc[2], ucdf)
	}
	uc <- c(chival, 2 * chival, 3 * chival)
	llim <- (1 - conf)/2

	while (pchisq(chival, df, uc[1]) < llim && uc>prec) {
		uc <- c(uc[1]/4, uc[1], uc[3])
	}
	while (pchisq(chival, df, uc[3]) > llim) {
		uc <- c(uc[1], uc[3], uc[3] + chival)
	}
	diff <- 1
	count<-0
	while (diff > prec) {
		if (pchisq(chival, df, uc[2]) < llim) 
			uc <- c(uc[1], (uc[1] + uc[2])/2, uc[2])
		else uc <- c(uc[2], (uc[2] + uc[3])/2, uc[3])
		diff <- abs(pchisq(chival, df, uc[2]) - llim)
		lcdf <- pchisq(chival, df, uc[2])
		count<-count+1
		if(count>1000){
			warning("Convergence not reached in chi.noncentral.conf .")
			uc[2]<-Inf
			lcdf<-0
			break
		}
	}
	result <- rbind(result, c(uc[2], lcdf))
	rownames(result) <- c("Lower", "Upper")
	colnames(result) <- c("Non-Central", "%")
	result
}



#Adapted from
# g.test V3.3 Pete Hurd Sept 29 2001. phurd@ualberta.ca
likelihood.test <- function(x, y = NULL, conservative=FALSE)
{
	DNAME <- deparse(substitute(x))
	if (is.data.frame(x)) x <- as.matrix(x)
	if (is.matrix(x)) {
		if (min(dim(x)) == 1) 
		x <- as.vector(x)
	}
	if (!is.matrix(x) && !is.null(y)) {
		if (length(x) != length(y)) 
			stop("x and y must have the same length")
		DNAME <- paste(DNAME, "and", deparse(substitute(y)))
		OK <- complete.cases(x, y)
		x <- as.factor(x[OK])
		y <- as.factor(y[OK])
		if ((nlevels(x) < 2) || (nlevels(y) < 2)) 
			stop("x and y must have at least 2 levels")
		x <- table(x, y)
	}
	if (any(x < 0) || any(is.na(x))) 
		stop("all entries of x must be nonnegative and finite")
	if ((n <- sum(x)) == 0) 
		stop("at least one entry of x must be positive")

	if (!is.matrix(x))
		stop("Could not make a 2-dimensional matrix")
		
		
	#Test of Independence
	nrows<-nrow(x)
	ncols<-ncol(x)

	sr <- apply(x,1,sum)
	sc <- apply(x,2,sum)
	E <- outer(sr,sc, "*")/n

	# no monte-carlo
	# calculate G
	g <- 0
	for (i in 1:nrows){
		for (j in 1:ncols){
			if (x[i,j] != 0) g <- g + x[i,j] * log(x[i,j]/E[i,j])
		}
	}
	q <- 1
	if (conservative){ # Do Williams correction
		row.tot <- col.tot <- 0    
		for (i in 1:nrows){ row.tot <- row.tot + 1/(sum(x[i,])) }
		for (j in 1:ncols){ col.tot <- col.tot + 1/(sum(x[,j])) }
		q <- 1+ ((n*row.tot-1)*(n*col.tot-1))/(6*n*(ncols-1)*(nrows-1))
	}
	STATISTIC <- G <- 2 * g / q
	PARAMETER <- (nrow(x)-1)*(ncol(x)-1)
	PVAL <- 1-pchisq(STATISTIC,df=PARAMETER)
	if(!conservative)
		METHOD <- "Log likelihood ratio (G-test) test of independence without correction"
	else
		METHOD <- "Log likelihood ratio (G-test) test of independence with Williams' correction"

	names(STATISTIC) <- "Log likelihood ratio statistic (G)"
	names(PARAMETER) <- "X-squared df"
	names(PVAL) <- "p.value"
	structure(list(statistic=STATISTIC,parameter=PARAMETER,p.value=PVAL,
		method=METHOD,data.name=DNAME, observed=x, expected=E),
		class="htest")
}





print.contin.table<-function(x,digits=3,prop.r=TRUE,prop.c=TRUE,prop.t=TRUE,
						expected.n=FALSE,residuals=FALSE,std.residuals=FALSE,adj.residuals=FALSE,no.tables=FALSE,...){
	tab<-x
	for(index in 1:length(tab)){
		if(class(tab[[index]])!="single.table"){
			print(tab[[index]],...)
			next
		}else if(no.tables)
			next

		ColTotal <- "Column Total"
    	RowTotal <- "Row Total"
		t<-tab[[index]]$table
		RS<-tab[[index]]$row.sums
		CS<-tab[[index]]$col.sums
		CPR<-tab[[index]]$row.prop
		CPC<-tab[[index]]$col.prop
		CPT<-tab[[index]]$total.prop
		GT<-tab[[index]]$total
		expected<-tab[[index]]$expected
		ASR <- (t - expected)/sqrt(expected * 
            ((1 - RS/GT) %*% t(1 - CS/GT)))
		StdR<-(t - expected)/sqrt(expected)
		RTitleMargin<-10
    	CWidth <- max(digits + 4, c(nchar(t), nchar(dimnames(t)[[2]]), 
    	    nchar(RS), nchar(CS), nchar(RowTotal)))
		RWidth <- max(c(nchar(dimnames(t)[[1]]), nchar(ColTotal)))
		RWidth <- RWidth+RTitleMargin
    	RowSep <- paste(rep("-", CWidth + 2), collapse = "")
    	RowSep1 <- paste(rep("-", RWidth + 1), collapse = "")
    	SpaceSep1 <- paste(rep(" ", RWidth-RTitleMargin), collapse = "")
    	SpaceSep2 <- paste(rep(" ", CWidth), collapse = "")
		FirstCol <- formatC(paste(dimnames(t)[[1]]," Count   "), width = RWidth, format = "s")
    	ColTotal <- formatC(ColTotal, width = RWidth, format = "s")
    	RowTotal <- formatC(RowTotal, width = CWidth, format = "s")
		RowData<-names(dimnames(t))[1]
		ColData<-names(dimnames(t))[2]
		if(names(tab)[index]!="No Strata")
			cat(paste("\n",SpaceSep1,"           | -- ",sep="",collapse=""),"Stratum = ",
				names(tab)[index]," --", "\n",sep="")
    	out<-capture.output(cat(paste(SpaceSep1,"          ",sep="",collapse=""), "|", ColData, "\n"))
   		out<-c(out,capture.output(cat(cat(formatC(substr(RowData,0,RWidth-1), width = RWidth, format = "s"), 
    					sep = " | ", collapse = ""), cat(formatC(dimnames(t)[[2]], 
    					width = CWidth - 1, format = "s"), sep = "  | ", 
    					collapse = ""), cat(RowTotal, sep = " | ", collapse = "\n"), 
    					sep = "", collapse = "")))
    	out<-c(out,capture.output(cat(RowSep1, rep(RowSep, ncol(t) + 1), sep = "|", collapse = "\n")))
    	for (i in 1:nrow(t)) {
			out<-c(out,capture.output(cat(cat(FirstCol[i], sep = " | ", collapse = ""), 
    					cat(formatC(c(t[i, ], RS[i]), width = CWidth - 
    					1, format = "d"), sep = "  | ", collapse = "\n"), 
    					sep = "", collapse = "")))
    	    if (prop.r) 
				out<-c(out,capture.output(cat(cat(paste(SpaceSep1," Row %   "), sep = " | ", collapse = ""), 
    	        		cat(formatC(c(CPR[i, ] * 100, 100 * RS[i]/GT), 
    	            	width = CWidth - 1, digits = digits, format = "f"), 
        	        	sep = "% | ", collapse = "\n"), sep = "", 
       	         		collapse = "")))
	       	if (prop.c) 
				out<-c(out,capture.output(cat(cat(paste(SpaceSep1," Column %"), sep = " | ", collapse = ""), 
							cat(formatC(CPC[i, ] * 100, width = CWidth - 
							1, digits = digits, format = "f"), sep = "% | ", 
							collapse = ""), cat(SpaceSep2, sep = " | ", 
							collapse = "\n"), sep = "", collapse = "")))
	        if (prop.t) 
				out<-c(out,capture.output(cat(cat(paste(SpaceSep1," Total % "), sep = " | ", collapse = ""), 
							cat(formatC(CPT[i, ] * 100, width = CWidth - 
							1, digits = digits, format = "f"), sep = "% | ", 
							collapse = ""), cat(SpaceSep2, sep = " | ", 
							collapse = "\n"), sep = "", collapse = "")))
	        if (expected.n) 
				out<-c(out,capture.output(cat(cat(paste(SpaceSep1," Expected"), sep = " | ", collapse = ""), 
							cat(formatC(expected[i, ], digits = digits, 
							format = "f", width = CWidth - 1), sep = "  | ", 
							collapse = ""), cat(SpaceSep2, sep = " | ", 
							collapse = "\n"), sep = "", collapse = "")))
	        if (residuals) 
				out<-c(out,capture.output(cat(cat(paste(SpaceSep1," Residual"), sep = " | ", collapse = ""), 
							cat(formatC(t[i,]-expected[i, ], digits = digits, 
							format = "f", width = CWidth - 1), sep = "  | ", 
							collapse = ""), cat(SpaceSep2, sep = " | ", 
							collapse = "\n"), sep = "", collapse = "")))
	        if (adj.residuals)
				out<-c(out,capture.output(cat(cat(paste(SpaceSep1,"Adj Resid"), sep = " | ", collapse = ""), 
							cat(formatC(ASR[i,], digits = digits, 
							format = "f", width = CWidth - 1), sep = "  | ", 
							collapse = ""), cat(SpaceSep2, sep = " | ", 
							collapse = "\n"), sep = "", collapse = "")))
	        if (std.residuals)
				out<-c(out,capture.output(cat(cat(paste(SpaceSep1,"Std Resid"), sep = " | ", collapse = ""), 
							cat(formatC(StdR[i,], digits = digits, 
							format = "f", width = CWidth - 1), sep = "  | ", 
							collapse = ""), cat(SpaceSep2, sep = " | ", 
							collapse = "\n"), sep = "", collapse = "")))
	        out<-c(out,capture.output(cat(RowSep1, rep(RowSep, ncol(t) + 1), sep = "|", 
	            		collapse = "\n")))
		}
		out<-c(out,capture.output(cat(cat(ColTotal, sep = " | ", collapse = ""), cat(formatC(c(CS, 
	    				GT), width = CWidth - 1, format = "d"), sep = "  | ", 
	    				collapse = "\n"), sep = "", collapse = "")))
		if(prop.c)
	        out<-c(out,capture.output(cat(cat(paste(SpaceSep1," Column %"), sep = " | ", collapse = ""), 
						cat(formatC(CS/GT * 100, width = CWidth - 
						1, digits = digits, format = "f"), sep = "% | ", 
						collapse = ""), cat(SpaceSep2, sep = " | ", 
						collapse = "\n"), sep = "", collapse = "")))
		ConsoleWidth<-getOption("width")
		printOut<-function(StartingCharPos){
			FirstWidth<-nchar(FirstCol[1])+2
			ColWidth<-CWidth+3
			NCol<-length(lapply(strsplit(out[6],"|",fixed=TRUE),nchar)[[1]])-2
			WorkingNCol<-(ConsoleWidth-FirstWidth-1)%/%ColWidth
			for(i in 1:length(out)){
				cat(substr(out[i],0,FirstWidth),substr(out[i],
					StartingCharPos+FirstWidth,StartingCharPos+FirstWidth+WorkingNCol*ColWidth),
					"\n",sep="")
			}
			if((StartingCharPos+FirstWidth+WorkingNCol*ColWidth)<=(nchar(out[6])-1)){
				cat("\n",sep="")
				printOut(StartingCharPos+WorkingNCol*ColWidth)
			}
		}
		printOut(1)
	}
}

extract.counts<-function(tables){
	result<-list()
	for(j in 1:length(tables)){
		table<-tables[[j]]
		table.list<-lapply(table,function(x) if(class(x)=="single.table") x$table)
		table.list<-table.list[sapply(table.list,function(x) !is.null(x))]
		len<-length(table.list)
		rownam<-c()
		colnam<-c()
		for(i in 1:len){
			rownam<-union(rownam,rownames(table.list[[i]]))
			colnam<-union(colnam,colnames(table.list[[i]]))
		}
		ary<-array(0,c(length(rownam),length(colnam),len),list(rownam,colnam,names(table.list)))
		for(i in 1:len){
			ary[rownames(table.list[[i]]),colnames(table.list[[i]]),i]<-table.list[[i]]
		}
		result[[names(tables)[j]]]<-ary
	}
	result
}


add.cross.strata.test<-function(tables,name,htests,types=c("asymptotic","monte.carlo","exact")){
	types<-match.arg(types,c("asymptotic","monte.carlo","exact"),TRUE)
	if(is.function(htests))
		htests<-list(htests)
	if(length(htests)!=length(types))
		stop("type and tests must be the same length")
	if(class(tables)!="contingency.tables")
		stop("tables is not a contingency.tables object")

	count.tables<-extract.counts(tables)
	for(i in 1:length(tables)){
		if(class(tables[[i]])!="contin.table")
			next
		tests<-list(stratum = "Cross Strata",
							asymptotic=if("asymptotic" %in% types) 
											htests[[which(types=="asymptotic")]](count.tables[[names(tables)[i]]]) 
										else NA,
							monte.carlo=if("monte.carlo" %in% types) 
											htests[[which(types=="monte.carlo")]](count.tables[[names(tables)[i]]]) 
										else NA,
							exact=if("exact" %in% types) 
											htests[[which(types=="exact")]](count.tables[[names(tables)[i]]]) 
										else NA)
		tests<-tests[!is.na(tests)]
		test.l<-list(tests)
		if(is.null(tables[[i]]$cross.strata.tests)){
				tables[[i]]$cross.strata.tests<-list()
				class(tables[[i]]$cross.strata.tests)<-"contin.tests"
		}
		tables[[i]]$cross.strata.tests[name]<-list(test.l)
	}
	tables
}

add.mantel.haenszel<-function(tables,conservative=FALSE){
	add.cross.strata.test(tables,"Mantel-Haenszel",list(function(x) mantelhaen.test(x,correct=conservative)),"asymptotic")
}


add.test<-function(tables,name,htests,types=c("asymptotic","monte.carlo","exact")){
	
	types<-match.arg(types,c("asymptotic","monte.carlo","exact"),TRUE)
	if(is.function(htests))
		htests<-list(htests)
	if(length(htests)!=length(types))
		stop("type and tests must be the same length")
	if(class(tables)!="contingency.tables")
		stop("tables is not a contingency.tables object")
	for(i in 1:length(tables)){
		if(class(tables[[i]])!="contin.table")
			next
		tests<-list()
		for(j in 1:length(tables[[i]])){
			tab<-tables[[i]]
			if(class(tab[[j]])!="single.table")
				next
			tests[[j]]<-list(stratum = names(tab)[j],
								asymptotic=if("asymptotic" %in% types) 
												try(htests[[which(types=="asymptotic")]](tab[[j]]$table))
											else NA,
								monte.carlo=if("monte.carlo" %in% types) 
												try(htests[[which(types=="monte.carlo")]](tab[[j]]$table)) 
											else NA,
								exact=if("exact" %in% types) 
												try(htests[[which(types=="exact")]](tab[[j]]$table))
											else NA)
			tests[[j]]<-tests[[j]][!is.na(tests[[j]])]
			invalid<-sapply(tests[[j]],function(x) class(x)=="try-error")
			htestNA<-structure(list(statistic = NA, parameter = NA, 
        		p.value = NA, method = "", data.name = ""), class = "htest")
        	for(index in 1:length(tests[[j]]))
        		if(invalid[index])
					tests[[j]][[index]]<-htestNA
			
		}
		test.l<-list(tests)
		if(is.null(tables[[i]]$tests)){
				tables[[i]]$tests<-list()
				class(tables[[i]]$tests)<-"contin.tests"
		}
		tables[[i]]$tests[name]<-test.l
	}
	tables	
}

add.chi.squared<-function(tables, simulate.p.value = FALSE, B = 10000){


	if(simulate.p.value){
		chi.func<-list(function(x) chisq.test(x, correct=FALSE),
						function(x) chisq.test(x,simulate.p.value=TRUE,B=B))
		types<-c("asymptotic","monte.carlo")
	}else{
		chi.func<-list(function(x) chisq.test(x, correct=FALSE))
		types<-c("asymptotic")
	}
	add.test(tables,"Chi Squared",chi.func,types)
}

add.likelihood.ratio<-function(tables, conservative = FALSE, simulate.p.value = FALSE, B = 10000){

	if(simulate.p.value)
		warning("monte carlo not yet implemented")
	add.test(tables,paste(c("Likelihood",if(conservative) "(conservative)" else ""),collapse="")
					,function(x) likelihood.test(x,,conservative),"asymptotic")
}

add.fishers.exact<-function(tables,  simulate.p.value = FALSE, B = 10000){
	add.test(tables,"Fishers Exact"
					,function(x) {
	tst <- fisher.test(x,simulate.p.value=simulate.p.value,B=B)
	tst$statistic <- NA
	tst$parameter <- NA
	tst
	},"exact")	
}

add.correlation<-function(tables,method=c("spearman","kendall")){
	method<-match.arg(method)
	corr<-function(x) {
		colnames(x) <- 1:ncol(x)
		rownames(x) <- 1:nrow(x)
		tmp<-table.to.data(x)
		test<-cor.test(as.numeric(as.character(tmp[[1]])),as.numeric(as.character(tmp[[2]])),exact=FALSE,method=method)
		if(is.null(test$parameter))
			test$parameter<-NA
		if(is.null(test$conf.int))
			test$conf.int<-c(NA,NA)
		test
	}
	add.test(tables,paste(c(method,"'s Correlation"),collapse=""),
					corr,
					"asymptotic")
}

add.kruskal<-function(tables,nominal=c("both","rows","cols")){
	nominal<-match.arg(nominal)
	krusk.rows<-function(x) {
		tmp<-table.to.data(x)
		test<-kruskal.test(tmp[[2]],tmp[[1]])
		test
	}
	krusk.cols<-function(x) {
		tmp<-table.to.data(x)
		test<-kruskal.test(tmp[[1]],tmp[[2]])
		test
	}
	if(nominal %in% c("rows","both"))
		tables<-add.test(tables,"Kruskal-W (nominal rows)",
					krusk.rows,
					"asymptotic")
	if(nominal %in% c("cols","both"))
		tables<-add.test(tables,"Kruskal-W (nominal cols)",
					krusk.cols,
					"asymptotic")
	tables
}



contin.tests.to.table<-function(tests,test.digits=3,...){
	if(class(tests)!="contin.tests")
		stop("tests are not of class contin.tests")
	cat("\n\n\n",sep="")
	any.mc<-any(sapply(tests , function(test) sapply(test, function(x) !is.null(x$monte.carlo))))
	any.exact<-any(sapply(tests , function(test) sapply(test, function(x) !is.null(x$exact))))
	any.a<-any(sapply(tests , function(test) sapply(test, function(x) !is.null(x$asymptotic))))
	any.aest<-any(sapply(tests , function(test) sapply(test,
					function(x) sapply(x[c("asymptotic")], 
					function(y) {!is.null(y) && !is.null(y$estimate)&& !is.na(y$estimate)}))))
	n.col<-0
	if(any.a){
		n.col<-3
		est.start<-NULL
		if(any.aest){
			aest.div<-n.col+1
			aest.start<-n.col+2
			n.col<-n.col+5
		}
	}
	if(any.mc){
		mc.div<-n.col+1
		mc.start<-n.col+2	
		n.col<-n.col+4
	}
	if(any.exact){
		exact.div<-n.col+1
		exact.start<-n.col+2
		n.col<-n.col+4
	}

	n.tests<-length(tests)*max(sapply(tests,function(test) length(test)))
	table.form<-as.table(matrix(NA,ncol=n.col,nrow=n.tests+1))
	rownames(table.form)[1]<-"Test"
	for(test.index in 1:length(tests)){
		for(stratum.index in 1:length(tests[[test.index]])){
			row.index<-test.index+length(tests)*(stratum.index-1)+1
			test.list<-tests[[test.index]][[stratum.index]]
			if(test.index==1 && test.list$stratum!="No Strata" && test.list$stratum!="Cross Strata"){
				rownames(table.form)[row.index]<-paste(
									c("Stratum =",test.list$stratum,"  ",
									names(tests)[[test.index]]),collapse=" ")
			}else
				rownames(table.form)[row.index]<-names(tests)[[test.index]]
			if(!is.null(test.list$asymptotic)){
				table.form[row.index,1]<-round(test.list$asymptotic$statistic,test.digits)
				table.form[row.index,2]<-round(test.list$asymptotic$parameter,test.digits)
				table.form[row.index,3]<-round(test.list$asymptotic$p.value,test.digits)
			
				if(is.finite(table.form[row.index,3]) && table.form[row.index,3]<10^(-test.digits))
					table.form[row.index,3]<-paste(c("<0.",rep("0",test.digits-1),"1"),collapse="")
			}
			if(!is.null(test.list$monte.carlo)){
				table.form[row.index,mc.start]<-round(test.list$monte.carlo$statistic,test.digits)
				table.form[row.index,mc.start+1]<-round(test.list$monte.carlo$parameter,test.digits)
				table.form[row.index,mc.start+2]<-round(test.list$monte.carlo$p.value,test.digits)		
				if(is.finite(table.form[row.index,mc.start+2]) && 
						table.form[row.index,mc.start+2]<(10^-test.digits))
					table.form[row.index,mc.start+2]<-paste(c("<0.",rep("0",test.digits-1),"1"),
																collapse="")	
			}
			if(!is.null(test.list$exact)){
				table.form[row.index,exact.start]<-round(test.list$exact$statistic,test.digits)
				table.form[row.index,exact.start+1]<-round(test.list$exact$parameter,test.digits)
				table.form[row.index,exact.start+2]<-round(test.list$exact$p.value,test.digits)		
				if(is.finite(table.form[row.index,exact.start+2]) && 
						table.form[row.index,exact.start+2]<(10^-test.digits))
					table.form[row.index,exact.start+2]<-paste(c("<0.",rep("0",test.digits-1),"1"),
																collapse="")	
			}
			if(!is.null(test.list$asymptotic$estimate) && !is.null(names(test.list$asymptotic$estimate))){
				test<- test.list$asymptotic
				lim<-names(test$conf.int)
				table.form[row.index,aest.start]<-names(test$estimate)
				table.form[row.index,aest.start+1]<-round(test$estimate,test.digits)
				if(!is.null(test$conf.int)){
					if(!is.na(test$conf.int[1]))
						table.form[row.index,aest.start+2]<-paste(c(round(test$conf.int[1],test.digits)," (",
														as.numeric(lim[1])*100,")"),collapse="")
					if(!is.na(test$conf.int[1]))	
						table.form[row.index,aest.start+3]<-paste(c(round(test$conf.int[2],test.digits)," (",
														as.numeric(lim[2])*100,")"),collapse="")
				}
			}
		}
	}
	if(any.a){
		colnames(table.form)[1:3]<-c("Large Sample","","")
		table.form[1,1:3]<-c("Statistic","DF","p-value")
	}
	if(any.mc){
		colnames(table.form)[mc.start:(mc.start+2)]<-c("Monte Carlo","","")
		colnames(table.form)[mc.div]<-"|"
		table.form[1,mc.start:(mc.start+2)]<-c("Statistic","DF","p-value")
		table.form[,mc.div]<-"|"
	}
	if(any.exact){
		colnames(table.form)[exact.start:(exact.start+2)]<-c("Exact","","")
		colnames(table.form)[exact.div]<-"|"
		table.form[1,exact.start:(exact.start+2)]<-c("Statistic","DF","p-value")
		table.form[,exact.div]<-"|"
	}
	if(any.aest){
		table.form[1,aest.start:(aest.start+3)]<-c("Effect Size","est.","Lower (%)","Upper (%)")
		colnames(table.form)[c(aest.div,aest.start:(aest.start+3))]<-""
		table.form[,aest.div]<-"|"
	}
	
	rownames(table.form)<-formatC(rownames(table.form),width=max(nchar(rownames(table.form))))
	
	table.form
}

print.contin.tests<-function(x,test.digits=3,...){
	tests<-x
	class(tests)<-"contin.tests"
	table.form<-contin.tests.to.table(tests,test.digits,...)
	div<-which(colnames(table.form)=="|")
	div<-c(div,dim(table.form)[2]+1)
	if(dim(table.form)[2]<=7){
		print(table.form)
	}else{
		start.ind<-1
		if(div[1]==1)
			div<-div[-1]
		for(ind in div){
			tmp<-table.form[,start.ind:(ind-1)]
			print(as.table(tmp))
			cat(rep("-",nchar(rownames(tmp)[1])),"\n",sep="")
			start.ind<-ind+1
		}
	}
cat("\n\n\n")
}


contingency.tables<-function (row.vars, col.vars,stratum.var, data=NULL, missing.include=FALSE ) {
	arguments <- as.list(match.call()[-1])
	if(missing(row.vars) || missing(col.vars))
		stop("Please specify the row variables (row.vars), and column variables (col.vars)")
	row.vars <- eval(substitute(row.vars),data, parent.frame())
	col.vars <- eval(substitute(col.vars),data, parent.frame())
	if(length(dim(row.vars))<1.5){
		row.vars <- as.data.frame(row.vars)
		names(row.vars) <- as.character(arguments$row.vars)
	}
	if(length(dim(col.vars))<1.5){
		col.vars<-as.data.frame(col.vars)
		names(col.vars) <- as.character(arguments$col.vars)
	}
	
	#print(row.vars)
	if(!missing(stratum.var))
		stratum.var<-eval(substitute(stratum.var),data, parent.frame())
	else
		stratum.var<-NULL
    vector.x <- FALSE
	num.row.vars<-dim(row.vars)[2]
	num.col.vars<-dim(col.vars)[2]


	##takes a dataframe whose first two column are the row and column vectors for the table.
	##optionally a third column indicates stratum
	single.table<-function(dat,dnn){
		x<-dat[[1]]
		y<-dat[[2]]
		if(is.null(stratum.var))
			stratum.var<-rep("No Strata",length(x))
        if (length(x) != length(y)) 
            stop("all row.vars and col.vars must have the same length")
        if (missing.include) {
            x <- factor(x, exclude = c())
            y <- factor(y, exclude = c())
			strata <- factor(stratum.var, exclude = c())
        }
        else {
            x <- factor(x)
            y <- factor(y)
			strata <- factor(stratum.var)
        }
		lev<-levels(strata)
		table.list<-list()
		for(level in lev){
			temp.x<-x[strata==level]
			temp.y<-y[strata==level]
        	t <- table(temp.x, temp.y,dnn=dnn)
    		RS <- rowSums(t)
    		CS <- colSums(t)
			t<-t[RS>0,CS>0,drop=FALSE]
    		RS <- rowSums(t)
    		CS <- colSums(t)
    		CPR <- prop.table(t, 1)
    		CPC <- prop.table(t, 2)
    		CPT <- prop.table(t)
			CST <- try(suppressWarnings(chisq.test(t, correct = FALSE)),silent=TRUE)
			if(class(CST)!="htest")
				CST<-list(expected=t*NA)
    		GT <- sum(t)
    		if (length(dim(x) == 2)) 
        		TotalN <- GT
    		else TotalN <- length(temp.x)
			table.list[[level]]<-list(table=t,row.sums=RS,col.sums=CS,
										total=GT,row.prop=CPR,col.prop=CPC,total.prop=CPT,
										expected=CST$expected)
			class(table.list[[level]])<-"single.table"
		}
		class(table.list)<-"contin.table"
		table.list
	}
	result<-list()
	count<-1;
	for(i in 1:num.row.vars){
		for(j in 1:num.col.vars){
			result[[paste(names(row.vars)[i],"by",names(col.vars)[j])]]<-single.table(
						data.frame(row.vars[,i],col.vars[,j]),
								dnn=c(names(row.vars)[i],names(col.vars)[j]))
			count<-count+1
		}
	}
	class(result)<-"contingency.tables"
	result
}

print.contingency.tables<-function(x,digits=3,prop.r=TRUE,prop.c=TRUE,prop.t=TRUE,expected.n=FALSE,no.tables=FALSE,...){
	tables<-x
	ConsoleWidth<-getOption("width")
	MainSep<-paste(rep("=",ConsoleWidth),collapse="")
	SubSep<-paste(c(rep(" ",15),rep("=",ConsoleWidth-30),rep(" ",15)),collapse="")
	cat(MainSep,"\n",sep="")
	cat("\n",sep="")
	for(i in 1:length(tables)){
		cat(SubSep,"\n",sep="")
		name<-paste(c("Table: ",names(tables)[i]),collapse="")
		len<-nchar(name)+2
		space<-floor((ConsoleWidth-len)/2)
		cat(paste(rep(" ",max(0,space-10)),collapse=""),paste(rep("=",10),collapse="")," ",name,
			" ",paste(rep("=",10),collapse=""),paste(rep(" ",max(0,space-10)),collapse=""),"\n",sep="")
		print(tables[[i]],digits,prop.r,prop.c,prop.t,expected.n,no.tables=no.tables,...)
		cat("\n",sep="")
	}
	cat("\n",MainSep,"\n",sep="")
}




