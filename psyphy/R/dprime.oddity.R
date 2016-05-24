`dprime.oddity` <-
function(Pc.tri) {
	if (Pc.tri < 1/3) stop("Only valid for Pc.tri > 1/3")
	root3 <- sqrt(3)
	root2.3 <- sqrt(2)/root3
	est.dp <- function(dp) {
		pr <- function(x) {
			2 *(pnorm(-x * root3 + dp * root2.3) + 
			   pnorm(-x * root3 - dp * root2.3)) * dnorm(x)
			}
		Pc.tri - integrate(pr, lower=0, upper=Inf)$value
		}
	dp.res <- uniroot(est.dp, interval = c(0,10))
	dp.res$root
}

