setClass("anoint.fit",
	representation(
		K = "numeric",
		responsiveness = "list",
		tests = "list",
		pvalues = "list",
		fits = "list"
	)
)


anointfit.table <- 	function(x,...){
		test.table <- cbind(unlist(x@tests))
		row.names(test.table) <- c(
			"OBO",
			"OBO (adj.)",
			"UIM",
			"PIM (exact)",
			"PIM (approx)",
			"PIM/OBO (adj.)",
			"PIM/UIM"
		)
		colnames(test.table) <- "Global null rejected"

test.table
}
	
	
setMethod("print","anoint.fit",function(x,...) print(anointfit.table(x)))

setMethod("show","anoint.fit",function(object) print(anointfit.table(object)))

anointfit.summary <- function(object,...){
	
	test.table <- anointfit.table(object)
	
	pvalues <- cbind(unlist(object@pvalues))
	
	row.names(pvalues) <- c(
		"OBO (max)",
		"UIM",
		"PIM (exact)",
		"PIM (approx)"
	)
	
	colnames(pvalues) <- "Global LRT (p-value)"

	print(test.table)
	cat("\n")
	print(pvalues)
	
	invisible(list(tests=object@tests,pvalues=object@pvalues))	
}

setMethod("summary","anoint.fit",anointfit.summary)

anointfit.fit <- function(object,type=c("obo","uim","pim.exact","pim.approx")){
	object@fits[type]
}

setGeneric("fits",function(object,...)standardGeneric("fits"))

setMethod("fits","anoint.fit",anointfit.fit)
