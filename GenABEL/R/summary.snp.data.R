"summary.snp.data" <-
		function(object, ... ) {
	if (!is(object,"snp.data")) stop("wrong object class, should be snp.data")
	
	
#	if (class(object) == "snp.data") {
#		gtNrow <- dim(object)[1]
#		gtNcol <- dim(object)[2]
#		genodata <- object@gtps
#	} else if (class(object) == "matrix") {
#		gtNrow <- dim(object)[1]
#		gtNcol <- dim(object)[2]
#		genodata <- object
#		storage.mode(genodata) <- "double"
#	} else if (class(object)=="databel") {
#		gtNrow <- dim(object)[1]
#		gtNcol <- dim(object)[2]
#		genodata <- object@data
#	} else {
#		stop(paste("object's class not recognised ('",class(object),"')",sep=""))
#	}
#	print(class(object))
	
#	print(class(object))
#	print(class(object@gtps))
	if (class(object@gtps) != "snp.mx" & 
			class(object@gtps) != "matrix" &
			class(object@gtps) != "databel") 
		stop(paste("object@gtps's class not recognised ('",
						class(object@gtps),"')",sep=""))
	if (class(object@gtps) == "databel") if (!require(DatABEL)) 
			stop ("this function requires DatABEL package to be installed");
	
	gtNrow <- dim(object)[1]
	gtNcol <- dim(object)[2]
	res <- .Call("iteratorGA", object@gtps, 
			as.integer(gtNrow), as.integer(gtNcol),
			as.character("snp_summary_exhwe"),
			"R", # output type 
			as.integer(2), #margin 
			as.integer(1), # STEP 
			as.integer(1),999)
	
	
	#res <- .C("snp_summary_exhwe",as.raw(object@gtps),as.integer(object@nids),
	#		as.integer(object@nsnps), out = double(object@nsnps*9), PACKAGE="GenABEL")$out
	
	
	#dim(res) <- c(object@nsnps,9)
	#res <- as.data.frame(res)
	res[,9] <- pchisq(res[,9],1,lower.tail=F)
# X-chromosome
	if (any(chromosome(object) == "X")) {
		vec <- (chromosome(object) == "X")
		if (any(male(object) == 0)) {
			oX <- object[which(male(object) != 1),which(chromosome(object) == "X")]			
			#resX <- .C("snp_summary_exhwe",as.raw(oX@gtps),as.integer(oX@nids),as.integer(oX@nsnps), out = double(oX@nsnps*9), PACKAGE="GenABEL")$out
			gtNrow <- dim(oX)[1]
			gtNcol <- dim(oX)[2]
			resX <- .Call("iteratorGA", oX@gtps, 
					as.integer(gtNrow), as.integer(gtNcol),
					as.character("snp_summary_exhwe"),
					"R", # output type 
					as.integer(2), #margin 
					as.integer(1), # STEP 
					as.integer(1),999)			
			res[vec,7] <- resX[(nsnps(oX)*6+1):(nsnps(oX)*7)]
			res[vec,8] <- resX[(nsnps(oX)*7+1):(nsnps(oX)*8)]
			res[vec,9] <- pchisq(resX[(nsnps(oX)*8+1):(nsnps(oX)*9)],1,lower.tail=F)
			rm(oX,vec,resX);gc(verbose=FALSE)
		} else {
			res[vec,7] <- rep(1,sum(vec))
			res[vec,8] <- rep(1,sum(vec))
			res[vec,9] <- rep(1,sum(vec))
		}
	}
# Y- chromosome
	if (any(chromosome(object) == "Y")) {
		vec <- (chromosome(object) == "Y")
		res[vec,7] <- rep(1,sum(vec))
		res[vec,8] <- rep(1,sum(vec))
		res[vec,9] <- rep(1,sum(vec))
	}
# mtDNA
	if (any(chromosome(object) == "mt")) {
		vec <- (chromosome(object) == "mt")
		res[vec,7] <- rep(1,sum(vec))
		res[vec,8] <- rep(1,sum(vec))
		res[vec,9] <- rep(1,sum(vec))
	}
	
# report
	rownames(res) <- snpnames(object)
	colnames(res) <- c("NoMeasured","CallRate","Q.2",
			"P.11","P.12","P.22","Pexact","Fmax","Plrt")
	#res$Plrt[res$Plrt==0] <- 9.99e-17;
	res <- cbind(annotation(object),res)
	res
}
