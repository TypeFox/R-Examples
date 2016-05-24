.isREMLFit <-
isREML <- 
function(x) {
	#if (inherits(x, "mer")) return (x@dims[["REML"]] != 0)
	if (inherits(x, "merMod")) return (lme4::isREML(x))
	if (inherits(x, c("lme", "gls", "gam")) && !is.null(x$method))
		return(x$method %in% c("lme.REML", "REML"))
	#if (inherits(x, c("lmer", "glmer")))
		#return(x@status["REML"] != 0)
	return(FALSE)
}

isGEE <- 
function(object) 
inherits(object, c("geeglm", "geese", "gee", "geem",
	"yagsResult"))
