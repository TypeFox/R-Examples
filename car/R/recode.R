# recode function (J. Fox)
# last modified 2014-08-04 by J. Fox

recode <- function(var, recodes, as.factor.result, as.numeric.result=TRUE, levels){
  lo <- -Inf
  hi <- Inf
	recodes <- gsub("\n|\t", " ", recodes)
	recode.list <- rev(strsplit(recodes, ";")[[1]])
	is.fac <- is.factor(var)
	if (missing(as.factor.result)) as.factor.result <- is.fac
	if (is.fac) var <- as.character(var)
	result <- var
	for (term in recode.list){
		if (0 < length(grep(":", term))) {
			range <- strsplit(strsplit(term, "=")[[1]][1],":")
			low <- try(eval(parse(text=range[[1]][1])), silent=TRUE)
            if (class(low) == "try-error"){
                stop("\n  in recode term: ", term, 
                    "\n  message: ", low)
            }
			high <- try(eval(parse(text=range[[1]][2])), silent=TRUE)
		    if (class(high) == "try-error"){
		        stop("\n  in recode term: ", term, 
		            "\n  message: ", high)
		    }
			target <- try(eval(parse(text=strsplit(term, "=")[[1]][2])), silent=TRUE)
		    if (class(target) == "try-error"){
		        stop("\n  in recode term: ", term, 
		            "\n  message: ", target)
		    }
			result[(var >= low) & (var <= high)] <- target
		}
		else if (0 < length(grep("^else=", squeezeBlanks(term)))) {
			target <- try(eval(parse(text=strsplit(term, "=")[[1]][2])), silent=TRUE)
		    if (class(target) == "try-error"){
		        stop("\n  in recode term: ", term, 
		            "\n  message: ", target)
		    }
			result[1:length(var)] <- target
		}
		else {
			set <- try(eval(parse(text=strsplit(term, "=")[[1]][1])), silent=TRUE)
		    if (class(set) == "try-error"){
		        stop("\n  in recode term: ", term, 
		            "\n  message: ", set)
		    }
			target <- try(eval(parse(text=strsplit(term, "=")[[1]][2])), silent=TRUE)
		    if (class(target) == "try-error"){
		        stop("\n  in recode term: ", term, 
		            "\n  message: ", target)
		    }
			for (val in set){
				if (is.na(val)) result[is.na(var)] <- target
				else result[var == val] <- target
			}
		}
	}
	if (as.factor.result) {
		result <- if (!missing(levels)) factor(result, levels=levels) 
			else as.factor(result)
	}
	else if (as.numeric.result && (!is.numeric(result))) {
		result.valid <- na.omit(result)
		opt <- options("warn"=-1)
		result.valid <- as.numeric(result.valid)
		options(opt)
		if (!any(is.na(result.valid))) result <- as.numeric(result)
	}
	result
}

Recode <- function (...) car::recode(...)

