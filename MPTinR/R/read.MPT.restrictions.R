

#rules for restriction files:
# 1. inequalities first
# 2. If a variable appears in an inequality restriction, it can not be on a LHS of any restriction
# 3. If a variable appears on RHS of an equality restriction, it can not appear on LHS of an equality restriction,
# 4. Parameters can only be restricted to certain numbers.

.read.MPT.restrictions <- function(filename) {
	#min.restriction <- c(0, 0.001) 
	#max.restriction <- 0.99999
	if (!is.list(filename)) {
		whole <- readLines(filename)
		model <- vector("list", length(whole))
		whole <- gsub("#.*", "", whole)
		c2 <- 1
		for (c1 in 1:length(whole)) {
			if (!(grepl("^[[:space:]]*$", whole[c1]))) {
				model[[c2]] <- whole[c1]
				fin <- c2
				c2 <- c2 + 1
			}
		}
		if(!exists("fin")) stop("Restrictions file seems to be empty!")
		tmp.restrictions <- model[1:fin]
	} else {
		if (length(filename) == 0) stop("Restrictions list seems to be empty!")
		tmp.restrictions <- filename
	}
	
	restrictions <- list()
	
   
	if (sum(grepl("[>/\\+\\*\\!]", unlist(tmp.restrictions)))) stop("Error getting Restrictions: Non supported operators (>, +, *, /, !) found in restriction file.")
	restrictions.split <- sapply(tmp.restrictions,strsplit, split = "")
	#if (length(tmp.restrictions) > 1) {
		c.inequality <- vapply(lapply(restrictions.split, grepl, pattern = "<"), sum, 0)
	#} else c.inequality <- vapply(lapply(list(restrictions.split), grepl, pattern = "<"), sum, 0)
	#if (length(tmp.restrictions) > 1) {
		n.equality <- vapply(lapply(restrictions.split, grepl, pattern = "="), sum, 0)
	#} else n.equality <- vapply(lapply(list(restrictions.split), grepl, pattern = "="), sum, 0)
	if (any(pmax(c.inequality, n.equality) < (c.inequality + n.equality))) stop("Error getting Restrictions: A line contains = and <!")
	prevEquality <- FALSE
	c.x.all <- 1
	alsoUsed <- NULL
	for (c.restr in 1:length(tmp.restrictions)) {
		usedParamsList <- sapply(restrictions, function(x) return(x[1]))
		usedParams <- c(sapply(restrictions, function(x) return(x[1])), alsoUsed)
        #browser()
		if (c.inequality[c.restr] > 1) {
			if (prevEquality) stop("Error getting Restrictions: Inequality after Equality.")
			tmpStr <- strsplit(tmp.restrictions[[c.restr]], "[\t <]")[[1]]
			tmpStr <- grep(".", tmpStr, value=TRUE)
			refParam <- tmpStr[length(tmpStr)]
			alsoUsed <- c(alsoUsed, refParam)
			if (refParam %in% usedParamsList) refParam <- restrictions[[match(refParam, usedParamsList)]][2]
			tmpStr.2 <- tmpStr
			tmpStr <- tmpStr[-length(tmpStr)]
			if (sum(tmpStr %in% usedParams)) stop("Error getting Restrictions: In inequality restrictions, any variable can only appear once not as the rightmost variable (and only the first time it appears).")
			tmpStr <- rev(tmpStr)
			tmpRestr <- vector("list", length(tmpStr))
			dummy <- ""
			for (c.ineq in 1:length(tmpStr)) {
				dummy.x <- paste("hank", c.x.all, "y", c.ineq, sep = "")
				dummy <- paste(dummy, dummy.x, sep = "")
				tmpRestr[[c.ineq]] <- c(tmpStr[c.ineq], paste(refParam, " * ", dummy, sep = ""), "<", tmpStr.2[length(tmpStr.2)-(c.ineq-1)])
				dummy <- paste(dummy, " * ", sep = "")
			}
		}
		else {
			if (c.inequality[c.restr] == 1) {
				if (prevEquality) stop("Error getting Restrictions: Inequality after Equality.")
				tmpStr <- strsplit(tmp.restrictions[[c.restr]], "[\t <]")[[1]]
				tmpStr <- grep(".", tmpStr, value=TRUE)
				if (tmpStr[1] %in% usedParams) stop("Error getting Restrictions: In inequality restrictions, any variable can only appear once not as the rightmost variable (and only the first time it appears).")
				if (tmpStr[2] %in% usedParamsList) {
					tmpRestr <- list(c(tmpStr[1], paste(restrictions[[match(tmpStr[2], usedParams)]][2], " * x", c.x.all, sep = ""), "<", restrictions[[match(tmpStr[2], usedParams)]][2]))
				}
				else {
					tmpRestr <- list(c(tmpStr[1], paste(tmpStr[2], " * hank", c.x.all, sep = ""), "<", tmpStr[2]))
				}
			}
			else {
				if (!(sum(grepl("=", tmp.restrictions[[c.restr]])))) stop("Error getting Restrictions: No recognized operator (< or =) in at least one of the restrictions.")
				tmpStr <- strsplit(tmp.restrictions[[c.restr]], "[\t =]")[[1]]
				tmpStr <- grep(".", tmpStr, value=TRUE)
				pos.ref <- length(tmpStr)
				tmpRestr <- vector("list", pos.ref - 1)
				if (!(grepl("[[:alpha:]]", tmpStr[pos.ref]))) {
					tmpStr[pos.ref] <- as.numeric(tmpStr[pos.ref])
# 					if (as.numeric(tmpStr[pos.ref]) > min.restriction[1] & as.numeric(tmpStr[pos.ref]) < min.restriction[2]) {
# 						warning(paste("Restriction starting with ", tmpStr[1], ": Constant is in the not allowed interval from ", min.restriction[1], " to ", min.restriction[2], " and therefore changed to ", 0, ".", sep = ""))
# 						tmpStr[pos.ref] <- 0
# 					}
# 					if (as.numeric(tmpStr[pos.ref]) > max.restriction) {
# 						warning(paste("Restriction starting with ", tmpStr[1], ": Constant is bigger than ", max.restriction, " and therefore changed to ", max.restriction, ".", sep = ""))
# 						tmpStr[pos.ref] <- max.restriction
# 					}
					if (as.numeric(tmpStr[pos.ref]) >= 1 | as.numeric(tmpStr[pos.ref]) <= 0) {
					  message(paste("Restriction starting with ", tmpStr[1], ": Constant is either equal to 0 or 1 or outside the interval from 0 to 1. This may lead to problems.", sep = ""))
					}
				}
				for (ci in 1:(pos.ref-1)) {
						if (tmpStr[ci] %in% usedParams) stop("Error getting Restrictions: In equality restrictions, no variable can appear on RHS of restrtcitons and then on LHS of restriction.")
				}
				alsoUsed <- c(alsoUsed, tmpStr[pos.ref])
				for (ci in 1:(pos.ref-1)) {
					if (tmpStr[pos.ref] %in% usedParamsList) {
                        if (restrictions[[match(tmpStr[pos.ref], usedParams)]][3] == "<") tmpRestr[[ci]] <- c(tmpStr[ci], restrictions[[match(tmpStr[pos.ref], usedParams)]][2], "<", restrictions[[match(tmpStr[pos.ref], usedParams)]][4])
						else tmpRestr[[ci]] <- c(tmpStr[ci], restrictions[[match(tmpStr[pos.ref], usedParams)]][2], "=")
					}
					else {
						tmpRestr[[ci]] <- c(tmpStr[ci], tmpStr[pos.ref], "=")
					}
				}
				prevEquality <- TRUE
			}
		}
		c.x.all <- c.x.all + 1
		restrictions <- c(restrictions, tmpRestr)
	}
	#recover()
	return(restrictions)
}

