# corrects values < -1 to -1 and values > 1 to 1
clamp <- function(...) {
	x <- c(...)
	x[x < -1] <- -1
	x[x >  1] <-  1
	return(x)
}



# helper-functions
posOrNA <- function(x) {
	return(ifelse(x>=0, x, NA))
}



# simple wrapper: formats a number in f.2 format
f2 <- function(x, digits=3, prepoint=0) {
	if (length(dim(x)) == 2) {
		apply(x, 2, function(x2) {gsub("0.", ".", sprintf(paste("%",prepoint,".",digits,"f",sep=""), x2) , fixed=TRUE)})
	} else {
		gsub("0.", ".", sprintf(paste("%",prepoint,".",digits,"f",sep=""), x) , fixed=TRUE)
	}
}


checkVar <- function(x, minVar=0) {
	if (is.null(minVar)) return(FALSE)
	if (is.na(minVar)) return(FALSE)
	if (is.null(x)) return(TRUE)
	if (is.nan(x)) return(TRUE)
	if (is.na(x)) return(TRUE)
	if (x < minVar) return(TRUE)
	return(FALSE)
}




# @param formula The RR formula - see ?RR
# @param data A long format data frame
# @example
# data(multiLikingLong)
# multi2 <- exportForREML(liking_a ~ perceiver.id*target.id|group.id, data=multiLikingLong)
exportForREML <- function(formula, data) {
		# Sort data according to the matrix formulation of the SRM
		f1 <- formula
		lhs <- strsplit(gsub(" ", "", as.character(f1)[2], fixed=TRUE), "+", fixed=TRUE)[[1]]
		rhs <- strsplit(gsub(" ", "", as.character(f1)[3], fixed=TRUE),"\\*|\\|", perl=TRUE)[[1]]

		actor.id <- rhs[1]
		partner.id <- rhs[2]
		if (length(rhs)>=3) {group.id <- rhs[3]} else {
			data[, ".group"] <- 1
			group.id=".group"
		}
		
		dat1 <- subset(data, data[, actor.id] < data[, partner.id])
		dat1 <- dat1[order(dat1[, group.id], dat1[, actor.id], dat1[, partner.id]), ]
		dat2 <- subset(data, data[, actor.id] > data[, partner.id])
		dat2 <- dat2[order(dat2[, group.id], dat2[, partner.id], dat2[, actor.id]), ]
		dat3   <- rbind(dat1, dat2)
		return(dat3)
}


clearLongData <- function(formule, data, minData=1) {
	ll1 <- long2matrix(formule, data, reduce=TRUE, minData=minData)
	
	lhs <- strsplit(gsub(" ","",as.character(formule)[2], fixed=TRUE), "+", fixed=TRUE)[[1]]
	rhs <- strsplit(gsub(" ","",as.character(formule)[3], fixed=TRUE),"\\*|\\|", perl=TRUE)[[1]]
	
	var.id <- lhs
	actor.id <- rhs[1]
	partner.id <- rhs[2]
	if (length(rhs)>=3) {group.id <- rhs[3]} else {group.id="group.id"}
	
	ll2 <- ldply(ll1, function(x) {
		matrix2long(x, new.ids=FALSE, var.id=var.id)
	})
	colnames(ll2)[1:3] <- c(group.id, actor.id, partner.id)

	return(ll2)
}

