`randomNames` <- 
function(
	n, 
	gender, 
	ethnicity, 
	which.names="both", 
	name.order="last.first", 
	name.sep=", ") {

	.N <- V1 <- NULL ## To prevent R CMD check warnings
                        
	first_names <- function(tmp.gender, tmp.ethnicity, tmp.number) {
		tmp.gender <- tmp.gender[1]; tmp.ethnicity <- tmp.ethnicity[1]
		tmp <- randomNames::randomNamesData[[paste("first_names_e", tmp.ethnicity, "_g", tmp.gender, sep="")]]
		suppressWarnings(sample(rownames(tmp), tmp.number, replace=TRUE, prob=tmp))
	} 

	last_names <- function(tmp.ethnicity, tmp.number) {
		tmp.ethnicity <- tmp.ethnicity[1]
		tmp <- randomNames::randomNamesData[[paste("last_names_e", tmp.ethnicity, sep="")]]
		suppressWarnings(sample(rownames(tmp), tmp.number, replace=TRUE, prob=tmp))
	}

	if (missing(n) & missing(gender) & missing(ethnicity)) n <- 1
	if (missing(n)) n <- NA
	if (missing(gender)) gender <- NA
	if (missing(ethnicity)) ethnicity <- NA

	tmp.length <- max(n, length(gender), length(ethnicity), na.rm=TRUE)
	if (is.na(tmp.length)) tmp.length <- 1

	if (missing(gender)) gender <- round(runif(tmp.length))
	if (missing(ethnicity)) ethnicity <- round(runif(tmp.length,min=1,max=5))

	gender <- rep(gender, length=tmp.length)
	ethnicity <- rep(ethnicity, length=tmp.length)

	if (!all(gender %in% c(0,1))) {
		tmp <- rep(NA, length(gender))
		tmp[grep("Male|M", gender, ignore.case=TRUE)] <- 0
		tmp[grep("Female|F", gender, ignore.case=TRUE)] <- 1
		gender <- tmp
	}

	if (!all(ethnicity %in% 1:5)) {
		tmp <- rep(NA, length(ethnicity))
		tmp[grep("Indian|Alaska|Native American", ethnicity, ignore.case=TRUE)] <- 1
		tmp[grep("Hawaii|Pacific|Asian", ethnicity, ignore.case=TRUE)] <- 2
		tmp[grep("Black|African", ethnicity, ignore.case=TRUE)] <- 3
		tmp[grep("Latino|Hispanic", ethnicity, ignore.case=TRUE)] <- 4
		tmp[grep("White|Caucasian", ethnicity, ignore.case=TRUE)] <- 5
		ethnicity <- tmp
	}

	gender[is.na(gender)] <- round(runif(length(gender[is.na(gender)]))) 
	ethnicity[is.na(ethnicity)] <- round(runif(length(ethnicity[is.na(ethnicity)]),min=1,max=5)) 

	tmp.dt <- data.table(seq.int(tmp.length), gender, ethnicity)
	setkeyv(tmp.dt, c("gender", "ethnicity"))

	if (which.names=="first" | which.names=="both") {
		tmp.dt$tmp_first <- tmp.dt[,first_names(gender, ethnicity, .N), by=list(gender, ethnicity)]$V1
	}

	if (which.names=="last" | which.names=="both") {
		tmp.dt$tmp_last <- tmp.dt[,last_names(ethnicity, .N), by=list(ethnicity)]$V1
	}

	setkey(tmp.dt, V1)
	if (which.names=="first") return(tmp.dt$tmp_first)
	if (which.names=="last") return(tmp.dt$tmp_last)
	if (which.names=="both" & name.order=="last.first") return(paste(tmp.dt$tmp_last, tmp.dt$tmp_first, sep=name.sep))
	if (which.names=="both" & name.order=="first.last") return(paste(tmp.dt$tmp_first, tmp.dt$tmp_last, sep=name.sep))

} ## END randomNames Function

