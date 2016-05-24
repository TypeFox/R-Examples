"mps" <- function(x, method="whittaker", all=FALSE){
	METHODS <- c("whittaker", "additive", "inverse.whittaker", "harrison", "diserud", "harrison.turnover", "williams")
	method <- pmatch(method, METHODS)
	# x <- as.matrix(x > 0)*1
	gamma <- ncol(x[,colSums(x)>0])
	alphas <- colSums(x>0)
	n <- nrow(x)
	betaW <- gamma / (sum(alphas)/n)
	betaADD <- gamma - (sum(alphas)/n)
	betaWrev <- (sum(alphas)/n) / gamma
	mpsH <- (betaW - 1) / (n - 1)
	mpsDO <- 1 - mpsH
	mpsHt <- (gamma/max(alphas) - 1) / (n - 1)
	mpsWt <- 1 - (max(alphas)/gamma)
	res <- c(betaW, betaADD, betaWrev, mpsH, mpsDO, mpsHt, mpsWt)
	if(length(res)!=7){
		res <- rep(NA, 7)
		}
	names(res) <- c("Whittaker's beta", "additive partitioning", "inverse Whittaker's beta","Harrison mps", "Diserud & Odegaard mps", "Harrison mpst", "Williams mpst")
	if(!all){
		res <- res[method]
	}
	res	
}