# Workaround: When @export does not recognize that this is a S3-method, you need the extra @method statement
#' @method print RRmulti
#' @export 
print.RRmulti <- function(x, ...) {
	print.RR(x, ...)
}

#' @method print RRuni
#' @export
print.RRuni <- function(x, ...) {
	print.RR(x, ...)
}

#' @method print RRbi
#' @export
print.RRbi <- function(x, ...) {
	print.RR(x, ...)
}


# x muss hier direkt auf das univariate-Objekt verweisen
print.uni <- function(x, ..., measure=NA, digits=3, r.names=NULL, minVar=0) {
	
	# print descriptivers for multi group
	if (length(x$groups) > 1) {
		groupsizes <- laply(x$groups, function(y) return(attr(y, "group.size")))
		av.groupsize <- round(mean(groupsizes), 2)
		
		 cat(paste("Group descriptives: n = ",length(x$groups),"; average group size = ",av.groupsize, "; range: ", min(groupsizes), "-", max(groupsizes), "\n"))
	} else {
		cat("Round robin analysis for a single group; using the formula of Lashley & Bond (1997).\n\n")
	}
	
	uni <- round(x$varComp[,2:ncol(x$varComp)], digits)	
	
	# remove SEVAR columns form output
	uni <- uni[, !grepl("SEVAR", colnames(uni))]
	
	if (checkVar(uni[1, 2], minVar)) {uni[5, 2:5] <- NA}
	if (checkVar(uni[2, 2], minVar)) {uni[5, 2:5] <- NA}
	if (checkVar(uni[3, 2], minVar)) {uni[6, 2:5] <- NA}
		
	if (is.na(measure)) {
		measure <- localOptions$style
	} else {
		measure <- match.arg(measure, c("behavior", "perception", "metaperception"))
	}
	
	if (!is.null(r.names)) {rownames(uni) <- r.names} else {
		if (measure == "behavior") rownames(uni) <- unilabels_b
		if (measure == "perception") rownames(uni) <- unilabels_p
		if (measure == "metaperception") {
			warning("Warning: the current RR-object only consists of a single variable. Labels for metaperception are only provided when two variables are specified.", call.=FALSE)
			rownames(uni) <- unilabels_b
		}
	}
	
	print(uni)
	
	
	# Actor effect reliability
	if (!is.null(x$effects[,grep(localOptions$suffixes[1], colnames(x$effects), fixed=TRUE)])) cat(paste(role[[measure]][1], "effect reliability:",f2(attr(x$effects[,grep(localOptions$suffixes[1], colnames(x$effects), fixed=TRUE)], "reliability"), 3), "\n"))
	
	# Partner effect reliability
	if (!is.null(x$effects[,grep(localOptions$suffixes[2], colnames(x$effects), fixed=TRUE)])) cat(paste(role[[measure]][2], "effect reliability:",f2(attr(x$effects[,grep(localOptions$suffixes[2], colnames(x$effects), fixed=TRUE)], "reliability"), 3), "\n"))
	
	# Relationship effect reliability
	if (!is.null(attr(x$effectsRel$relationship, "reliability"))) cat(paste(role[[measure]][3], "effect reliability:",f2(attr(x$effectsRel$relationship, "reliability"), 3), "\n"))
	
	selfCor(x, measure=measure)
}



# Here the default print method for RR-objects gets overwritten, so that 
# the information in the RR-class is displayed in a convenient way
print.RR <- function(x, ..., measure1=NA, measure2=NA, digits=3, measure=NULL) {
	
	if (is.na(measure1)) {
		measure1 <- localOptions$style
	} else {
		measure1 <- match.arg(measure1, c("behavior", "perception", "metaperception"))
	}
	if (is.na(measure2)) {
		measure2 <- measure1
	} else {
		measure2 <- match.arg(measure2, c("behavior", "perception", "metaperception"))
	}
	
	cat("Round-Robin object ('RR'), calculated by TripleR\n")
	cat("------------------------------------------------\n")
	cat(x$anal.type)
	cat("\n\n")
	
	
	if (!is.null(measure)) {measure1 <- measure}
	
	# bivariate case
	if (length(x$univariate) == 2) {
		
		uni <- lapply(x$univariate, function(x) return(x))
		bi <- round(x$bivariate[,2:ncol(x$bivariate)], digits)
		
		# remove SEVAR columns form output
		bi <- bi[, !grepl("SEVAR", colnames(bi))]
				
		# Erase bivariate correlations for variance components < minVar
		if (checkVar(uni[[1]]$varComp[1, 2], x$minVar)) {bi[c(1,3), 2:5] <- NA}
		if (checkVar(uni[[1]]$varComp[2, 2], x$minVar)) {bi[c(2,4), 2:5] <- NA}
		if (checkVar(uni[[1]]$varComp[3, 2], x$minVar)) {bi[c(5,6), 2:5] <- NA}
		if (checkVar(uni[[2]]$varComp[1, 2], x$minVar)) {bi[c(1,4), 2:5] <- NA}
		if (checkVar(uni[[2]]$varComp[2, 2], x$minVar)) {bi[c(2,3), 2:5] <- NA}
		if (checkVar(uni[[2]]$varComp[3, 2], x$minVar)) {bi[c(5,6), 2:5] <- NA}
		                                     
		r.names1 <- r.names2 <- NULL
		if (measure1 == "behavior" & measure2 == "behavior") {
			rownames(bi) <- bilabels_bb
		} else
   		if (measure1 == "behavior" & measure2 == "perception") {
				rownames(bi) <- bilabels_bp
		} else
		if (measure1 == "perception" & measure2 == "behavior") {
				rownames(bi) <- bilabels_pb
		} else
		if (measure1 == "perception" & measure2 == "perception") {
				rownames(bi) <- bilabels_pp
		} else
		if (measure1 == "perception" & measure2 == "metaperception") {
			r.names1 <- unilabels_b_meta1
			r.names2 <- unilabels_b_meta2
			rownames(bi) <- bilabels_meta
		} else {
			stop("This combination of measurement labels does not fit.")
		}
		cat(paste("Univariate analyses for:", attr(uni[[1]], "varname"), "\n"))
		cat("---------\n")		
		print.uni(uni[[1]], measure=measure1, r.names=r.names1, minVar=x$minVar, digits=digits)
		cat("\n\n")
		
		cat(paste("Univariate analyses for:", attr(uni[[2]], "varname"), "\n"))
		cat("---------\n")
		print.uni(uni[[2]], measure=measure2, r.names=r.names2, minVar=x$minVar, digits=digits)
		cat("\n\n")
		cat(paste0("Bivariate analyses:", "\n"))
		cat("---------\n")
		
		print(bi)
		
		if (length(uni[[1]]$groups) != length(uni[[2]]$groups)) {
			warning(paste("Note: Univariate analyses of both variables are based on different numbers of groups. Bivariate analyses therefore are based on the common groups of both variables (n=",min(length(uni[[1]]$groups), length(uni[[2]]$groups)),")", sep=""), call.=FALSE)
		}
		
	} else
	
	# univariate case
	{
		cat(paste("Univariate analyses for:", attr(x, "varname"), "\n"))
		cat("---------\n")
		print.uni(x, measure=measure1, minVar=x$minVar, digits=digits)
	}
}
