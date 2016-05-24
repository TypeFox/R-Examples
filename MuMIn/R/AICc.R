# AIC = -2 * LL + 2 * n
# LL = -1/2 * log(dev) * n + C
# AIC = log(dev) * n + C + 2 * n

.get.nobs <- function(llik, error = FALSE) {
	no <- attr(llik, "nall")
	if(is.null(no)) no <- attr(llik, "nobs")
	if (error && is.null(no)) stop("'logLik' object must have a \"nobs\" attribute")
	no
}

.aic <- function(objectlist, chat, k, REML, ICFun, ICName) {
	if(chat < 1) {
		warning("'chat' given is < 1, increased to 1")
		chat <- 1
	}
	npar.adj <- if(chat == 1) 0 else 1
	
	llCall <-  call("logLik", as.name("object"))
	if(!is.null(REML)) llCall$REML <- REML
	ll <- function(object) fixLogLik(NA, object)
	body(ll)[[2L]] <- llCall
	
	if(length(objectlist) > 1L) {
		lls <- lapply(objectlist, ll)
		val <- data.frame(df = vapply(lls, attr, 1, "df"),
						  ic = sapply(lls, ICFun, chat = chat, k = k, npar.adj = npar.adj))
		Call <- match.call(sys.function(sys.parent()), call = sys.call(sys.parent()))
		Call$chat <- Call$REML <- Call$k <- NULL
		dimnames(val) <- list(as.character(Call[-1L]), c("df", ICName))
		return(val)
	} else {
		return(ICFun(ll(objectlist[[1L]]), chat = chat, k = k, npar.adj = npar.adj))
	}	
}

`QAICc` <-
function(object, ..., chat, k = 2, REML = NULL) {
	.aic(list(object, ...), chat, k, REML, function(ll, chat, k, npar.adj) {
		no <- .get.nobs(ll, error = FALSE)
		# df is the number of parameters plus 1 for estimating c-hat
		df <- attr(ll, "df") + npar.adj
		#neg2ll <- log(deviance(object)) * n # + Constant...
		neg2ll <- -2 * c(ll)
		ret <- (neg2ll / chat) + (k * df) * (1 + ((df + 1) / (no - df - 1)))
		return (ret)
	}, "QAICc") 
}

`QAIC` <-
function(object, ..., chat, k = 2, REML = NULL) {
	.aic(list(object, ...), chat, k, REML, function(ll, chat, k, npar.adj) {
		df <- attr(ll, "df") + npar.adj
		neg2ll <- -2 * c(ll)
		#ret <- (neg2ll * no / chat) + k * df
		#ret <- -2 * ll / chat + k * df
		ret <- neg2ll / chat + k * df
		return (ret)
	}, "QAIC") 
}

`AICc` <-
function(object, ..., k = 2, REML = NULL) {
	.aic(list(object, ...), 1, k, REML, function(ll, chat, k, npar.adj) {
		no <- .get.nobs(ll, error = FALSE)
		df <- attr(ll, "df")
		ret <- (-2 * c(ll)) + (k * df) * (1 + ((df + 1) / (no - df - 1)))
		return (ret)
	}, "AICc") 
}


#QAIC = -2log Likelihood/c-hat + 2K
#QAICc = -2log Likelihood/c-hat + 2K + 2K(K + 1)/(n - K - 1)
