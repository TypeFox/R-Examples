# CLASS FOR IPD_LME META-ANALYTIC MODEL

setClass("ipdlme",
	representation(
		fixef = "matrix",
		ranef = "matrix",
		vcov.fixef = "matrix",
		vcov.ranef = "list",
		sigma2 = "numeric",
		VarCorr = "matrix",
		convergence.trace = "list",
		converged = "logical",
		n.iter = "numeric",
		max.iter = "numeric",
		tol = "numeric",
		df = "numeric"
	)
)

# METHODS

print_ipdlme <- function(x,...)
{
  	cat("Fixed effects:\n")
	print(x@fixef)
	
	cat("\nRandom effects:\n")
	print(x@ranef)
	
	cat("\nResidual variances:\t",x@sigma2)
	
	cat("\nBetween-study variance:\n")
	
	print(x@VarCorr)
}

setMethod("print", "ipdlme", print_ipdlme)
setMethod("show", "ipdlme", function(object) print_ipdlme(object))

setMethod("summary", signature(object = "ipdlme"),
	  function(object, ...)
      {

          fcoef <- object@fixef
          vcov <- object@vcov.fixef
          corF <- sqrt(diag(vcov))

           coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF) #, DF = DF)
           coefs <- coefs[, 1:2, drop = FALSE]
           stat <- coefs[,1]/coefs[,2]
           pval <- 2*pt(abs(stat), df = object@df, lower = FALSE)
           coefs <- cbind(coefs, "t value" = stat, "Pr(>|t|)" = pval)
  
		 print(coefs)  
         
         cat("\nRandom effects:\n")
		print(object@ranef)
	
		cat("\nResidual variance:\t",object@sigma2,"\n")
	
		cat("\nBetween-study variance:\n")
	
		print(object@VarCorr)


         invisible(coefs)
})
      
setMethod("confint", signature(object = "ipdlme"),
              function(object, parm, level = 0.95, ...){
              		
              		# parm ignored
               		
               		fcoef <- object@fixef
          			se <- sqrt(diag(object@vcov.fixef))
          			q <- qt(1-(1-level)/2,df=object@df)
          			lower <- fcoef-q*se
          			upper <- fcoef+q*se
          			
          			result <- cbind(fcoef, lower, upper)
          			colnames(result) <- c("Estimate","LowerCI","UpperCI")
          		result
             })
              
# NEW GENERICS
setMethod("fixef","ipdlme",function(object) object@fixef)
setMethod("ranef","ipdlme",function(object) object@ranef)
setMethod("coef","ipdlme",function(object) object@fixef)
setMethod("vcov","ipdlme",function(object) object@vcov.fixef)

# NEW GENERICS
setGeneric("Var", function(object) standardGeneric("Var"))
setGeneric("sigma2", function(object) standardGeneric("sigma2"))
setGeneric("VcovFixef", function(object,...) standardGeneric("VcovFixef"))
setGeneric("VcovRanef", function(object,...) standardGeneric("VcovRanef"))
setGeneric("convergence", function(object) standardGeneric("convergence"))
setGeneric("n.iter", function(object) standardGeneric("n.iter"))
setGeneric("MaxIter", function(object,...) standardGeneric("MaxIter"))
setGeneric("tol", function(object) standardGeneric("tol"))
setGeneric("converged", function(object) standardGeneric("converged"))


setMethod("Var","ipdlme",function(object) object@VarCorr)
setMethod("sigma2","ipdlme",function(object) object@sigma2)
setMethod("VcovFixef","ipdlme",function(object) object@vcov.fixef)
setMethod("VcovRanef","ipdlme",function(object) object@vcov.ranef)
setMethod("convergence","ipdlme",function(object) object@convergence.trace)
setMethod("n.iter","ipdlme",function(object) object@n.iter)
setMethod("MaxIter","ipdlme",function(object) object@max.iter)
setMethod("tol","ipdlme",function(object) object@tol)
setMethod("converged","ipdlme",function(object) object@converged)

