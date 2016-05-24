##' also called from ../data/qDiptab.R :
if(getRversion() < "2.13.0") {
   rdRDS <- function(..., package = "diptest")
    .readRDS(system.file(..., package=package))
} else
   rdRDS <- function(..., package = "diptest")
    readRDS(system.file(..., package=package, mustWork=TRUE))

dip.test <- function(x, simulate.p.value = FALSE, B = 2000)
{
    DNAME <- deparse(substitute(x))
    x <- sort(x[complete.cases(x)])
    stopifnot(is.numeric(x))
    n <- length(x) # *is* integer
    D <- dip(x)
    method <- "Hartigans' dip test for unimodality / multimodality"
    if(n <= 3) {
	P <- 1
    } else if(simulate.p.value) {
	method <- paste(method,
			"with simulated p-value\n\t (based on", B, "replicates)")
	P <- mean(D <= replicate(B, dip(runif(n))))
    } else {
	## Long "codetools-compatible" way of	data(qDiptab) :
	qDiptab <- rdRDS("extraData", "qDiptab.rds")
	dn <- dimnames(qDiptab)
	max.n <- max(nn <- as.integer(dn[["n"]]))
	P.s <- as.numeric(dn[["Pr"]])

        if(n > max.n) { ## extrapolate, or rather just use the last n as == "asymptotic"
            message("n = ",n," > max_n{n in table} = ",max.n,
                    " -- using that as asymptotic value.")
            n1 <- n0 <- max.n
            i2 <- i.n <- length(nn)
            f.n <- 0
        } else {
            n0 <- nn[i.n <- findInterval(n, nn)]
            n1 <- nn[(i2 <- i.n +1)]
            f.n <- (n - n0)/(n1 - n0)# in [0, 1]
        }
        ## Now "find" y-interval:
        y.0 <- sqrt(n0)* qDiptab[i.n, ]
	y.1 <- sqrt(n1)* qDiptab[i2 , ]
        sD  <- sqrt(n) * D
	P <- 1 - approx(y.0 + f.n*(y.1 - y.0), P.s, rule = 2,# <- [min, max]
			xout = sD)[["y"]]
    }
    structure(list(statistic = c(D = D), p.value = P, nobs = n,
		   alternative = "non-unimodal, i.e., at least bimodal",
		   method = method, data.name = DNAME), class = "htest")
}
