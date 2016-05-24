setClass("formula.anoint",
	representation(
		formula = "formula",
		uni = "list",
		prognostic = "formula",
		prognostic.trt = "formula",
		trt = "character",
		family = "character"
))


anoint.formula <- function(formula=y~(a+b)*trt,family="binomial"){
	
	# OBO LIST OF FORMULAS
	
formulas.onebyone <- function(formula=y~(a+b)*trt){
	
	# GET ONE-BY-ONE FORMULA
	factors <- as.character(formula)
	response <- factors[2]
	factors <- factors[3]
	
	factors <- strsplit(factors,split="\\*")[[1]]
	trt <- factors[2]
	factors <- factors[1]
	
	factors <- sub("^\\((.*)\\) ?","\\1",factors)
	factors <- strsplit(factors,split=" ")[[1]]
	factors <- factors[grep("[a-zA-z]",factors)]
	
	formulas <- sapply(paste(response,"~",factors,"*",trt),as.formula)
	
 formulas
}
	uni <- formulas.onebyone(formula)
	
	trt <- sub("(.*)(\\* ?)(.*)","\\3",as.character(formula)[3])
	remove.trt <- paste("~",sub("(.*)(\\*.*)","\\1",as.character(formula)[3]),
						sep="",collapse="")
	add.trt <- paste("~.",trt,sep="+",collapse="")

	prognostic <- update(formula,remove.trt)
	prognostic.trt <- update(prognostic,add.trt)
	
	new("formula.anoint",
		uni = uni,
		formula = formula,
		prognostic = prognostic,
		prognostic.trt = prognostic.trt,
		trt = trt,
		family = family
	)
	
}

setMethod("print","formula.anoint",
		function(x,...) print(x@formula,showEnv=FALSE)
)

setMethod("show","formula.anoint",
		function(object) print(object@formula,showEnv=FALSE)
)

setMethod("update","formula.anoint",
		function(object,formula,...){
			update.formula(object@formula,formula)
		}
)
