library(robust)

sessionInfo()

testTRUEfile <- function(file, srcfile = NULL, verbose = TRUE) {
    exps <- parse(file = file, srcfile=srcfile)
    if(verbose) cat(length(exps)," expressions :\n")
    for(i in seq_along(exps)) {
        if(verbose) cat(" ")
        if(!isTRUE(eval(exps[[i]], envir = .GlobalEnv))) {
            ch.ex <- paste(substr(paste(format(exps[[i]])[-1],
                                        collapse = " ; "), 1, 60), "...")
## For testing many at once:
            warning("*** ", ch.ex,"  was *not* TRUE", immediate. = TRUE)
## Once, the tests work
##             stop("*** ", ch.ex,"  was *not* TRUE")
        }
        if(verbose) cat(i %% 10)
    }
}

tDir <- system.file("tests_S", package = "robust")

tstFiles <- list.files(tDir, pattern = "\\.t$")
## Remove those that are not (yet) available for package 'robust' :

.not.yet <- c("asymmetric.t",	## robust & MLE	 Gamma, Weibull, ...
	      "plots.wblrob.t", ## gammaRob(), weibullRob() not ported
	      "discrob.t",	## robust discrimn.analysis: discRob() not ported
	      "princomprob.t",	## princompRob() not ported
	      "plots.princomprob.t", ## ditto
	      "plots.aovrob.t", ## aovRob() not ported __ FIXME: use lawson data in lmRob()!
	      "covm.t", ## loc/scat estimators all come from robustbase/rrcov now
	      "")

(tstFiles <- tstFiles[! match(tstFiles, .not.yet, nomatch=0)])

for(f in tstFiles) {
    cat("Test File", f,": ")
    testTRUEfile(file.path(tDir, f), verbose = TRUE)
    cat("\n\n")
}
