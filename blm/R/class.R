setClass("lexpit",
	representation(
		coef.linear = "numeric",
		coef.expit = "numeric",
		vcov.linear = "matrix",
		vcov.expit = "matrix",
		formula.linear = "formula",
		formula.expit = "formula",
		df.residual = "numeric",
		p = "numeric",
		q = "numeric",
		data = "data.frame",
		which.kept = "numeric",
		y = "numeric",
		weights = "numeric",
		strata = "factor",
		converged = "logical",
		par.init = "numeric",
		loglik = "numeric",
		loglik.null = "numeric",
		barrier.value = "numeric",
		control.lexpit = "list"
	)
)

setMethod("print","lexpit",
			function(x,...){
				cat("Linear model","\n")
				print(x@formula.linear)
				print(x@coef.linear)
				cat("\n")
				cat("Expit model","\n")
				print(x@formula.expit)
				print(x@coef.expit)
			}
)


setMethod("show","lexpit",
			function(object) print(object))

aic <- function(object){
	if(class(object)=="lexpit")
	-2*object@loglik+2*(object@p+object@q)
	else
	-2*object@loglik+2*length(coef(object))
}

setMethod("confint","lexpit",
		function(object, parm, level = 0.95, ...){
			
			z <- qnorm(1-(1-level)/2)
			beta.se <- sqrt(diag(object@vcov.linear))
			gamma.se <- sqrt(diag(object@vcov.expit))
			l <- c(object@coef.linear-z*beta.se,object@coef.expit-z*gamma.se)
			u <- c(object@coef.linear+z*beta.se,object@coef.expit+z*gamma.se)
			
			mat <- cbind(c(object@coef.linear,object@coef.expit),l,u)
			
			if(!missing(parm)){
				mat <- rbind(mat[pmatch(parm, names(coef(object))),])
			}
			colnames(mat) <- c("Est.","Lower","Upper")
		mat
		}
)


setMethod("summary","lexpit",
			function(object,...){
			
			digits = max(3, getOption("digits") - 3)
			
			#LINEAR		
			t <- object@coef.linear/sqrt(diag(object@vcov.linear))
			p <- sapply(t,function(x)2*pt(abs(x),df=object@df.residual,lower.tail=FALSE))
			
			linear.mat <- cbind(
				object@coef.linear,
				sqrt(diag(object@vcov.linear)),
				t,p
			)	
			
			
			row.names(linear.mat) <- names(object@coef.linear)
			colnames(linear.mat) <- c("Estimate","Std.Err","t value","Pr(>|t|)")
			
#			print.linear.mat <- format(linear.mat, digits=digits)
			
			#EXPIT
			t <- object@coef.expit/sqrt(diag(object@vcov.expit))
			p <- sapply(t,function(x) 2*pt(abs(x),df=object@df.residual,lower.tail=FALSE))
			
			expit.mat <- cbind(
				object@coef.expit,
				sqrt(diag(object@vcov.expit)),
				t,p
			)	
			
			print.expit.mat <- format(expit.mat, digits=digits)
			
			row.names(expit.mat) <- names(object@coef.expit)
			colnames(expit.mat) <- c("Estimate","Std.Err","t value","Pr(>|t|)")
			
	#		print.expit.mat <- format(expit.mat, digits = digits)
			
			cat("Linear effects:\n")
			printCoefmat(linear.mat)
			
			cat("\nExpit effects:\n")
			printCoefmat(expit.mat)
			
			cat("\nConverged:",object@converged,"\n")
			
invisible(list(linear=linear.mat,expit=expit.mat))
})

setMethod("coef","lexpit",
	function(object,...) c(object@coef.linear,object@coef.expit)
)

setMethod("vcov","lexpit",
	function(object,...) list(linear=object@vcov.linear,expit=object@vcov.expit)
)


###

setClass("blm",
	representation(
		coef = "numeric",
		vcov = "matrix",
		formula = "formula",
		df.residual = "numeric",
		data = "data.frame",
		which.kept = "numeric",
		y = "numeric",
		weights = "numeric",
		strata = "factor",
		converged = "logical",
		par.init = "numeric",
		loglik = "numeric",
		loglik.null = "numeric",
		barrier.value = "numeric"
	)
)


setMethod("print","blm",
			function(x,...){
				print(x@formula)
				print(x@coef)
			}
)


setMethod("show","blm",
			function(object) print(object))

setMethod("confint","blm",
		function(object, parm = NULL, level = 0.95, ...){
			
			z <- qnorm(1-(1-level)/2)
			beta.se <- sqrt(diag(object@vcov))
			l <- c(object@coef-z*beta.se)
			u <- c(object@coef+z*beta.se)
			
			mat <- cbind(object@coef,l,u)
			
		   if(!missing(parm)){
				mat <- rbind(mat[pmatch(parm, names(coef(object))),])
			}
			
			colnames(mat) <- c("Est.","Lower","Upper")
		mat
		}
)


setMethod("summary","blm",
			function(object,...){

			digits = max(3, getOption("digits") - 3)
	
			#LINEAR		
			t <- object@coef/sqrt(diag(object@vcov))
			p <- sapply(t,function(x) 2*pt(abs(x),df=object@df.residual,lower.tail=FALSE))
			
			linear.mat <- cbind(
				object@coef,
				sqrt(diag(object@vcov)),
				t,p
			)	
			
			print.linear.mat <- format(linear.mat, digits=digits)
			
			row.names(linear.mat) <- names(object@coef)
			colnames(linear.mat) <- c("Estimate","Std.Err","t value","Pr(>|t|)")
			
			printCoefmat(linear.mat)
			
			cat("\nConverged:",object@converged,"\n")
invisible(linear.mat)
})

setMethod("coef","blm",
	function(object,...) object@coef
	)

setMethod("vcov","blm",
	function(object,...) object@vcov
	)

setMethod("resid", "blm", function(object) object@y - predict(object))
setMethod("resid", "lexpit", function(object) object@y - predict(object))

setMethod("logLik","blm", function(object,..){
	LL <- object@loglik
	class(LL) <- "logLik"
	attr(LL,"df") <- object@df.residual
	attr(LL,"nobs") <- nrow(object@data)
LL	
})

setMethod("logLik","lexpit", function(object,..){
	LL <- object@loglik
	class(LL) <- "logLik"
	attr(LL,"df") <- object@df.residual
	attr(LL,"nobs") <- nrow(object@data)
LL	
})

setGeneric("leverage",function(object) standardGeneric("leverage"))

setMethod("leverage","blm",function(object) blm.leverage(object))
setMethod("leverage","lexpit",function(object) lexpit.leverage(object))

setGeneric("displacement",function(object) standardGeneric("displacement"))

setMethod("displacement","blm",function(object) C(object))
setMethod("displacement","lexpit",function(object) C(object))

setGeneric("model.formula",function(object) standardGeneric("model.formula"))

setMethod("model.formula","blm",function(object) object@formula)
setMethod("model.formula","lexpit",function(object) list(linear = object@formula.linear, expit = object@formula.expit))
