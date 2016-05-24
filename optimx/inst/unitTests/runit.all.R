# To run these tests:
#   library(optimx)
#   library(svUnit)
#   runit.all <- system.file("unitTests", "runit.all.R", package = "optimx")
#   source(runit.all); clearLog(); test.all()
#   Log()

test.all <- function() {

	## test 1

	checkTrue(require(optimx))

	## test 2

	genrose.f<- function(x, gs=NULL){ # objective function
	## One generalization of the Rosenbrock banana valley function (n parameters)
		n <- length(x)
		if(is.null(gs)) { gs=100.0 }
		fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
		return(fval)
	}

	genrose.g <- function(x, gs=NULL){
	# vectorized gradient for genrose.f
	# Ravi Varadhan 2009-04-03
		n <- length(x)
		if(is.null(gs)) { gs=100.0 }
		gg <- as.vector(rep(0, n))
		tn <- 2:n
		tn1 <- tn - 1
		z1 <- x[tn] - x[tn1]^2
		z2 <- 1 - x[tn]
		gg[tn] <- 2 * (gs * z1 - z2)
		gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
		return(gg)
	}

	genrose.h <- function(x, gs=NULL) { ## compute Hessian
	   if(is.null(gs)) { gs=100.0 }
		n <- length(x)
		hh<-matrix(rep(0, n*n),n,n)
		for (i in 2:n) {
			z1<-x[i]-x[i-1]*x[i-1]
			z2<-1.0-x[i]
			hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
			hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
			hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
			hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
		}
		return(hh)
	}

	startx<-4*seq(1:10)/3.
	ans8<-optimx(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, control=list(all.methods=TRUE, save.failures=TRUE, trace=0), gs=10)
	print(ans8)
	print(ans8[, "gevals"])
	print(ans8["spg", ])
	print(ans8, par.select = 1:3)
	print(ans8, best.only = TRUE)

	genrose.f<- function(x, gs=NULL){ # objective function
	## One generalization of the Rosenbrock banana valley function (n parameters)
		n <- length(x)
		if(is.null(gs)) { gs=100.0 }
		fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
		return(fval)
	}

	genrose.g <- function(x, gs=NULL){
	# vectorized gradient for genrose.f
	# Ravi Varadhan 2009-04-03
		n <- length(x)
		if(is.null(gs)) { gs=100.0 }
		gg <- as.vector(rep(0, n))
		tn <- 2:n
		tn1 <- tn - 1
		z1 <- x[tn] - x[tn1]^2
		z2 <- 1 - x[tn]
		gg[tn] <- 2 * (gs * z1 - z2)
		gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
		return(gg)
	}

	genrose.h <- function(x, gs=NULL) { ## compute Hessian
	   if(is.null(gs)) { gs=100.0 }
		n <- length(x)
		hh<-matrix(rep(0, n*n),n,n)
		for (i in 2:n) {
			z1<-x[i]-x[i-1]*x[i-1]
			z2<-1.0-x[i]
			hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
			hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
			hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
			hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
		}
		return(hh)
	}

	startx<-4*seq(1:10)/3.
	ans8<-optimx(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, 
		control=list(all.methods=TRUE, save.failures=TRUE, trace=0), gs=10)
	ans8.sum <- summary(ans8, order = value)[1, ]

	ans8.sum.target <- structure(list(p1 = 0.999999999920414, p2 = 0.999999999913547, 
	    p3 = 0.999999999885138, p4 = 0.999999999895513, p5 = 0.999999999967326, 
	    p6 = 1.00000000007184, p7 = 1.00000000015134, p8 = 1.00000000020547, 
	    p9 = 1.00000000022618, p10 = 1.00000000021113, value = 1, 
	    fevals = 147, gevals = 85, niter = NA_real_, convcode = 0, 
	    kkt1 = TRUE, kkt2 = TRUE, xtimes = 0.016), .Names = c("p1", 
	    "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10", "value", 
	    "fevals", "gevals", "niter", "convcode", "kkt1", "kkt2", "xtimes"
	    ), details = structure(list("Rvmmin", c(-2.90873103358124e-09, 
	    -1.04026964799964e-09, -4.07835964941952e-09, -4.75625960764814e-09, 
	    -2.02688776855642e-09, 2.5813320369885e-09, 4.34381553008344e-09, 
	    5.85728665721567e-09, 6.40616893084557e-09, -4.40224834363789e-09
	    ), structure(c(79.9999999843576, -39.9999999968166, 0, 0, 0, 
 	    0, 0, 0, 0, 0, -39.9999999968166, 101.999999983846, -39.9999999965419, 
	    0, 0, 0, 0, 0, 0, 0, 0, -39.9999999965419, 101.999999976613, 
	    -39.9999999954055, 0, 0, 0, 0, 0, 0, 0, 0, -39.9999999954055, 
	    101.99999997623, -39.9999999958205, 0, 0, 0, 0, 0, 0, 0, 0, -39.9999999958205, 
	    101.999999989285, -39.999999998693, 0, 0, 0, 0, 0, 0, 0, 0, -39.999999998693, 
	    102.000000011188, -40.0000000028736, 0, 0, 0, 0, 0, 0, 0, 0, 
	    -40.0000000028736, 102.000000028102, -40.0000000060534, 0, 0, 
	    0, 0, 0, 0, 0, 0, -40.0000000060534, 102.000000040265, -40.0000000082187, 
	    0, 0, 0, 0, 0, 0, 0, 0, -40.0000000082187, 102.000000045837, 
	    -40.0000000090471, 0, 0, 0, 0, 0, 0, 0, 0, -40.0000000090471, 
	    22), .Dim = c(10L, 10L)), c(178.067667089027, 166.6561012052, 
	    148.884023804866, 126.494474509405, 101.685420756996, 76.8972546720349, 
	    54.5792163971264, 36.9557928034118, 25.7800547148152, 1.99999408284202
	    ), "none"), .Dim = c(1L, 5L), .Dimnames = list("Rvmmin", c("method", 
	    "ngatend", "nhatend", "hev", "message"))), maximize = FALSE, npar = 10L, 
	    row.names = "Rvmmin", class = c("optimx", "data.frame"))

	# don't compare xtimes
	ans8.sum$xtimes <- ans8.sum.target$xtimes <- NULL
	checkEquals(ans8.sum, ans8.sum.target)


}

