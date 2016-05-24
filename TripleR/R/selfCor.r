# find the correct effect based on its attribute
findEff <- function(df, type) {
	for (i in 1:ncol(df)) {
		if (!is.null(attr(df[,i], "type"))) {
			if (attr(df[,i], "type") == type) {return(i)}	
		} 
	}
	return(NA)
}

noVar <- function(x) {
	if (length(x[!is.na(x)]) <= 2) {return(TRUE)}
	if (var(x, na.rm=TRUE) <= 0) {return(TRUE)}
	return(FALSE)
}

# x is an univariate RR-object
#' @export
selfCor <- function(x, digits=3, measure=NA) {
		
	print(attr(x, "group.id"))	
		
	## print (partial) correlations with self ratings:
	if (attr(x, "self") == TRUE) {
		
		if (is.na(measure)) {
			measure <- localOptions$style
		} else {
			measure <- match.arg(measure, c("behavior", "perception", "metaperception"))
		}

		if (noVar(x$effects[,findEff(x$effects, "self")])) {
			print("Warning: self ratings have zero variance or to few data points; skipping partial correlations with self ratings!")
			return(NULL);
		}
		
		
		if (length(x$groups) <= 1) {	
			
			cat("\n\nCorrelations with self ratings:\n")
			
			if (noVar(x$effects[,findEff(x$effects, "actor")])==FALSE) {
				c.a <- cor.test(x$effects[,findEff(x$effects, "actor")], x$effects[,findEff(x$effects, "self")], use="p")
				res.a <- data.frame(r=c.a$estimate, t=c.a$statistic, df=c.a$parameter, p=c.a$p.value)
			} else {
				res.a <- data.frame(r=NA, t=NA, df=NA, p=NA)
			}
			if (noVar(x$effects[,findEff(x$effects, "partner")])==FALSE) {
				c.p <- cor.test(x$effects[,findEff(x$effects, "partner")], x$effects[,findEff(x$effects, "self")], use="p")
				res.p <- data.frame(r=c.p$estimate, t=c.p$statistic, df=c.p$parameter, p=c.p$p.value)
			} else {
				res.p <- data.frame(r=NA, t=NA, df=NA, p=NA)
			}
			
			
			res0 <- rbind(res.a, res.p)
			res <- f2(res0)
			
		} else {
			cat("\n\nPartial correlations with self ratings (controlled for group membership):\n")
			
			c.a <- parCor(x$effects[,findEff(x$effects, "actor")], x$effects[,findEff(x$effects, "self")], x$effects[,2])
			c.p <- parCor(x$effects[,findEff(x$effects, "partner")], x$effects[,findEff(x$effects, "self")], x$effects[,2])
			res0 <- data.frame(r=c(c.a$par.cor, c.p$par.cor), t=c(c.a$t.value, c.p$t.value), df=c(c.a$df, c.p$df), p=c(c.a$p, c.p$p))
			res <- data.frame(r=f2(c(c.a$par.cor, c.p$par.cor), digits), t=f2(c(c.a$t.value, c.p$t.value)), df=c(c.a$df, c.p$df), p=f2(c(c.a$p, c.p$p), digits))
			
			if (noVar(x$effects[,findEff(x$effects, "actor")])==FALSE) {
				c.a <- parCor(x$effects[,findEff(x$effects, "actor")], x$effects[,findEff(x$effects, "self")], x$effects[,2])
				res.a <- data.frame(r=c.a$par.cor, t=c.a$t.value, df=c.a$df, p=c.a$p)
			} else {
				res.a <- data.frame(r=NA, t=NA, df=NA, p=NA)
			}
			if (noVar(x$effects[,findEff(x$effects, "partner")])==FALSE) {
				c.p <- parCor(x$effects[,findEff(x$effects, "partner")], x$effects[,findEff(x$effects, "self")], x$effects[,2])
				res.p <- data.frame(r=c.p$par.cor, t=c.p$t.value, df=c.p$df, p=c.p$p)
			} else {
				res.p <- data.frame(r=NA, t=NA, df=NA, p=NA)
			}
			
			
			res0 <- rbind(res.a, res.p)
			res <- f2(res0)
		}
				
		rownames(res) <- paste("self rating with",role[[measure]][1:2],"effect")

		if (measure=="perception") {
			rownames(res)[1] <- rownames(res0)[1] <-paste(rownames(res)[1], "(assumed similarity)")
			rownames(res)[2] <- rownames(res0)[2] <-paste(rownames(res)[2], "(self-other agreement)")
		}
		
		if (x$varComp[1,2] <= 0) {res[1,] <- NA}
		if (x$varComp[2,2] <= 0) {res[2,] <- NA}
		
		print(res, quote=FALSE)
		cat("\n\n")
		
		return(invisible(res0))
	}

}
