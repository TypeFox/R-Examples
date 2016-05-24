setMethod("getQuan", "Fa", function(obj) obj@n.obs)
setMethod("getCenter", "Fa", function(obj) obj@center)
setMethod("getLoadings", "Fa", function(obj) obj@loadings)
setMethod("getEigenvalues", "Fa", function(obj) obj@eigenvalues)
setMethod("getSdev", "Fa", function(obj) sqrt(obj@eigenvalues))
setMethod("getScores", "Fa", function(obj) obj@scores)
setMethod("getFa", "Fa", function(obj) {
    res = list(loadings = obj@loadings,
		uniquenesses = obj@uniquenesses,
		covariance = obj@covariance,
		correlation = obj@correlation,
		usedMatrix = obj@usedMatrix,
		scores = obj@scores,
		scoresMethod = obj@scoresMethod,
		sdev = sqrt(obj@eigenvalues))
    class(res) = "fa"
    res
})

##
## Follow the standard methods: show, print, plot
##
## ?show, we see show(object)
## For "show", the signature is object = "Fa", so we should use object as argument.
## Methods may be defined for arguments: object
## Use  showMethods("show")  for currently available ones.
setMethod("show", "Fa", function(object) myFaPrint(object))

myFaPrint = function(object, print.x = FALSE) {
	cl = object@call
    if(!is.null(cl)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }

    cat("Standard deviations:\n"); print(getSdev(object))
    cat("\nLoadings:\n");          print(getLoadings(object))

    if (print.x) {
        cat("\nRotated variables:\n"); print(getScores(object))
    }
    invisible(object)
}

## ?summary, we see summary(object)
## For "summary", the signature is object = "Fa", so we should use object as argument.
## Use  showMethods("summary")  for currently available ones.
setMethod("summary", "Fa", function(object){
	diag_S = diag(object@usedMatrix) 
	sum_rank = switch(object@method, 
                  pca = sum(diag_S),
                  pfa = sum(diag_S),
                  mle = nrow(getLoadings(object))) # mle uses the correlation matrix!
    m = object@factors; A = getLoadings(object)[,1:m]
    rowname = c("SS loadings", "Proportion Var", "Cumulative Var")
    colname = colnames(A)
    B = matrix(0, nrow = 3, ncol = m, dimnames = list(rowname, colname))
    for (i in 1:m){
       B[1,i] = sum(A[,i]^2)
       B[2,i] = B[1,i]/sum_rank
       B[3,i] = sum(B[1,1:i])/sum_rank
    }
    new("SummaryFa", faobj = object, importance = round(B,3))
})

## ?show, we see show(object)
## For "show", the signature is object = "SummaryFa", so we should use object as argument.
setMethod("show", "SummaryFa", function(object){

    cat("\nCall:\n")
    print(object@faobj@call)

    digits = max(3, getOption("digits") - 3)

    cat("Importance of components:\n")
    print(object@importance, digits = digits)
    invisible(object)
})

## ?print, we see print(x)
## For "print", the signature is x = "Fa", so we should use x as argument.
## Use  showMethods("print")  for currently available ones.
setMethod("print", "Fa", function(x, ...) myFaPrint(x, ...))

## ?predict, we see predict(object)
## For "predict", the signature is object = "Fa", so we should use object as argument.
## Use  showMethods("predict")  for currently available ones.
setMethod("predict", "Fa", function(object, newdata){
    res = getFa(object)
    if (missing(newdata)) { 
        if (!is.null(res$scores)) 
            return(res$scores)
        else stop("no scores and newdata")
    }
    ## !missing(newdata)
	## if newdata is not missing, newdata should be scaled before "predict".
    else { 
    F = switch(res$scoresMethod, 
             none = stop("no scores: try scoresMethod = 'regression' or 'Bartlett'"),
             regression = newdata %*% solve(res$usedMatrix) %*% res$loadings[],
             Bartlett = newdata %*% diag(1/res$uniquenesses) %*% res$loadings[] %*% solve(t(res$loadings[]) %*% diag(1/res$uniquenesses) %*% res$loadings[])
            ) 
    }
    
})

## ?plot, we see plot(x, y)
## For "plot", the signature is x = "Fa", y = "missing", so we should use x, y as argument.
## Use  showMethods("plot")  for currently available ones.
setMethod("plot", signature(x = "Fa", y = "missing"), 
function(x, y = "missing", which = c("factorScore", "screeplot"), choices = 1:2){
    if (which == "factorScore"){
        c.r = if (is.null(x@cov.control)) "classical" else "robust"
        main = eval(expression(paste("scores = ", x@scoresMethod, ", ", c.r, sep = "")))
        plot(x@scores[, choices], type = "n", main = main)
        text(x@scores[, choices][,1], x@scores[, choices][,2])
    }
    else # which == "screeplot"
        {stats:::screeplot.default(getFa(x), type = "lines", main = "Scree plot")} 
})
