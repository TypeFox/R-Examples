# To run these tests:
#   library(optimx)
#   library(svUnit)
#   runit.all <- system.file("unit1", "runit.1.R", package = "optimx")
#   source(runit.1); clearLog(); test.1()
#   Log()

test.1 <- function() {

	## test 1

	checkTrue(require(optimx))

	## test 2
	f1<-function(xx){ # function of one parameter
	   ((3*xx+2)*xx-5)*xx+4
	}
	ansone<-optimx(c(1), f1, control=list(all.methods=TRUE))
	ansone.sum<-summary(ansone, order= "value")[1,]




        ## test 3
	ansoneb<-optimx(c(1), f1, lower=c(-1), upper=c(10),control=list(all.methods=TRUE))
	ansoneb.sum<-summary(ansoneb, order= "value")[1,]
	ansoneb.sum.target<- structure(list(p1 = 0.555555555557684, value = 2.35390946502058, 
	    fevals = 17, gevals = NA_real_, niter = NA_real_, convcode = 0, 
	    kkt1 = TRUE, kkt2 = TRUE, xtimes = 0), .Names = c("p1", "value", 
	    "fevals", "gevals", "niter", "convcode", "kkt1", "kkt2", "xtimes"
	    ), details = structure(list("bobyqa", 5.52926311436595e-11, 
	    structure(14.0000000000441, .Dim = c(1L, 1L)), 14.0000000000441, "none"), 
	    .Dim = c(1L, 5L), .Dimnames = list("bobyqa", c("method", "ngatend", "nhatend",
	    "hev", "message" ))), maximize = FALSE, npar = 1L, row.names = "bobyqa", 
	    class = c("optimx", "data.frame"))


	# don't compare xtimes
	ansoneb.sum$xtimes <- ansoneb.sum.target$xtimes <- NULL
	checkEquals(ansoneb.sum, ansoneb.sum.target)


}

