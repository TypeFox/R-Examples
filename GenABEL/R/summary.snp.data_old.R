"summary.snp.data_old" <-
	function(object, ... ) {
		if (!is(object,"snp.data")) stop("wrong object class, should be snp.data")
		res <- .C("snp_summary_exhwe",as.raw(object@gtps),as.integer(object@nids),as.integer(object@nsnps), out = double(object@nsnps*9), PACKAGE="GenABEL")$out
		dim(res) <- c(object@nsnps,9)
		res <- as.data.frame(res)
		res[,9] <- pchisq(res[,9],1,lower.tail=F)
# X-chromosome
		if (any(object@chromosome == "X")) {
		  vec <- (object@chromosome == "X")
		  if (any(object@male == 0)) {
 		    oX <- object[(object@male != 1),(object@chromosome == "X")]
		    resX <- .C("snp_summary_exhwe",as.raw(oX@gtps),as.integer(oX@nids),as.integer(oX@nsnps), out = double(oX@nsnps*9), PACKAGE="GenABEL")$out
		    res[vec,7] <- resX[(oX@nsnps*6+1):(oX@nsnps*7)]
		    res[vec,8] <- resX[(oX@nsnps*7+1):(oX@nsnps*8)]
		    res[vec,9] <- pchisq(resX[(oX@nsnps*8+1):(oX@nsnps*9)],1,lower.tail=F)
		    rm(oX,vec,resX);gc(verbose=FALSE)
		  } else {
		    res[vec,7] <- rep(1,sum(vec))
		    res[vec,8] <- rep(1,sum(vec))
		    res[vec,9] <- rep(1,sum(vec))
		  }
		}
# Y- chromosome
		if (any(object@chromosome == "Y")) {
		  vec <- (object@chromosome == "Y")
		  res[vec,7] <- rep(1,sum(vec))
		  res[vec,8] <- rep(1,sum(vec))
		  res[vec,9] <- rep(1,sum(vec))
		}
# mtDNA
		if (any(object@chromosome == "mt")) {
		  vec <- (object@chromosome == "mt")
		  res[vec,7] <- rep(1,sum(vec))
		  res[vec,8] <- rep(1,sum(vec))
		  res[vec,9] <- rep(1,sum(vec))
		}

# report
		rownames(res) <- snpnames(object)
		colnames(res) <- c("NoMeasured","CallRate","Q.2",
				"P.11","P.12","P.22","Pexact","Fmax","Plrt")
		res$Plrt[res$Plrt==0] <- 9.99e-17;
		res <- cbind(annotation(object),res)
		res
	}
