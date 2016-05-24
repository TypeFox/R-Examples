## Variable extraction from large formulas (more efficient than terms.formula)
## '...' can notably contain "backtick=TRUE"
## Author : Sylvain Mareschal <maressyl@gmail.com>
var.formula <- function(formula, asCall=FALSE, ...) {
	# Checks
	if(!is(formula, "formula")) stop("'formula' must be a formula object")
	
	# Primary call
	if(length(formula) == 3)        { k <- formula[[3]]
	} else if(length(formula) == 2) { k <- formula[[2]]
	} else                          { stop("'formula' length is supposed to be 2 or 3")
	}
	
	# (Tail) Recursive name extraction
	out <- list()
	while(is.call(k)) {
		out <- c(k[[3]], out)
		k <- k[[2]]
	}
	out <- c(k, out)
	
	# Return mode
	if(isTRUE(asCall)) { return(out)
	} else             { return(sapply(out, deparse, ...))
	}
}
