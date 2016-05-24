specTest <- function(x, ...)
{
  UseMethod("specTest")
}

specTest.gmm <- function(x, ...)
	{
	if (x$infWmatrix == "ident")
		{
		gb <- colMeans(x$gt)
		j <- crossprod(gb,solve(x$v,gb))*x$n
		}
	else if ( (x$infVcov!="TrueFixed") & !is.null(x$weightsMatrix) )
		{
		gb <- colMeans(x$gt)
		j <- crossprod(gb,solve(x$v,gb))*x$n
		}
	else
		j <- x$objective*x$n

	J_test <- noquote(paste("J-Test: degrees of freedom is ",x$df,sep=""))
	j <- noquote(cbind(j, ifelse(x$df>0,pchisq(j,x$df, lower.tail = FALSE),"*******")))
	dimnames(j) <- list("Test E(g)=0:  ", c("J-test", "P-value"))
	ans<-list(ntest=J_test, test = j)
	class(ans) <- "specTest"
	ans
	}

print.specTest <- function(x, digits=5, ...)
	{
	cat("\n","## ",x$ntest," ##","\n\n")
	print.default(format(x$test, digits=digits),
                      print.gap = 2, quote = FALSE)
	cat("\n")
	invisible(x)
	}

specTest.gel <- function(x, ...)
	{
	n <- nrow(x$gt)
	khat <- x$khat
	gbar <- colMeans(x$gt)
	LR_test <- 2*x$objective*n*x$k2/(x$bwVal*x$k1^2)
	LM_test <- n*crossprod(x$lambda, crossprod(khat, x$lambda))/(x$bwVal^2)
	J_test <- n*crossprod(gbar, solve(khat, gbar))/(x$k1^2)
	test <- c(LR_test, LM_test, J_test)
	df <- (ncol(x$gt) - length(x$coef))
	ntest <- noquote(paste("Over-identifying restrictions tests: degrees of freedom is ", df, sep = ""))
	vptest <- pchisq(test,df,lower.tail = FALSE)
	test <- cbind(test,vptest)
	dimnames(test) <- list(c("LR test", "LM test", "J test"), c("statistics", "p-value"))	
	ans <- list(test = test, ntest = ntest)
	class(ans) <- "specTest"
	ans
	}

