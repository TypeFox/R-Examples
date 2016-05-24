asggm <- function(x, ...) UseMethod("asggm")

asggm.default <- function(x, iterations = 100000000, init = NULL, epsilon = 0.001, ...)
{
    x <- as.matrix(x)
    if (!is.null(init)) {
      init <- as.matrix(init)
    }

    result <- rCSL(x, iterations, init, epsilon, ...)

    class(result) <- "asggm"
    result
}

# TODO do we really want a formula interface?
# if so, we need to add the model response as a column in x instead of putting it in y
asggm.formula <- function(formula, data=list(), ...)
  {
    mf <- model.frame(formula=formula, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    y <- model.response(mf)
    
    result <- asggm.default(x, ...)
    result
  }  

  
genL = function(kNodes, spP) {
	L = matrix(0, kNodes, kNodes)
	for (i in 2:kNodes) {
		for (j in 1:(i-1)) {
			L[i,j] = (runif(1) < spP) * rnorm(1, 0, 1)
		}
	}
	# diag(L) = 1
	for (k in 1:kNodes) {
		L[k,k] = rnorm(1, 1, 0.1)
	}
	L
}

genData = function(L, nSamples) {
	kNodes = nrow(L)
	xTemp = matrix(rnorm(kNodes*nSamples, mean = 0, sd = 1), nrow = nSamples, ncol = kNodes)
	#### x = t(t(solve(L))%*%t(xTemp))
	x = xTemp%*%solve(L)
	
	x
}

rCSL = function(x, iterations = 500, init = NULL, epsilon = 1e-5, ansL = NULL){
    # TODO, compute init if it's null (init is L matrix)
    # TODO call C++ CSL
    # TODO return something
	
	LInit = init
	nSamples = nrow(x)
	kNodes = ncol(x)
	S = cov(x) * (nSamples-1)
	
	if (is.null(init)) {
		# LInit = .Call("olsInit", x, PACKAGE = "AdaptiveSparsity")
		# LInit = LInit[["OLS L"]]
		# LInit = t(chol(solve(S)))

		require(MASS)
		require(Matrix)
		Sinv = ginv(S)
		LUSinv = expand(lu(Sinv))
		LInit = as.matrix(LUSinv$L)
	}
	if (is.null(ansL)) {
		ansQ = LInit%*%t(LInit)
	} else {
		ansQ = ansL%*%t(ansL)
	}
	
	print(norm(ansQ-LInit%*%t(LInit),'f'))

	CSLOut = .Call("CSL", iterations, LInit, kNodes, nSamples, epsilon, S, ansQ, PACKAGE = "AdaptiveSparsity")
	CSLOut
	
}
