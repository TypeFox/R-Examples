
chisqTest <- function(draws=1e5,	# number of (total) draws
                      bins=200,     	# number of equally-spaced bins
                      edge=7, 		# cutoff for binning at +/- edge
                      seed=123456789,
                      steps=50,     	# resolution (number of rows until draws)
                      generators=c("Ziggurat", "MT", "LZLLV",
                                   "GSL", "QL", "Gretl"),
                      showplot=interactive()) {

    gr <- seq(-edge,edge,length=bins+1) # bins+1 'borders' defining the grids
    ##d  <- 2*binedge/bins              # difference between grids (not used)
    pv <- diff(pnorm(gr))               # expectect probability in each cell using
                                        # difference in probability distr. at point
    pv <- pv/sum(pv)			# normalise

    op <- options("warn"=-1)            # suppress warning of chisq ties
    res <- mclapply(generators, FUN=function(g, seed) {
        #cat("Running ", g, "\n")
        mat <- ziggbin(bins, draws, g, seed, edge, steps)
        vals <- apply(mat, 1, FUN=function(row, pv) {
            z <- chisq.test(row, p=pv)$statistic
        }, pv)
        vals
    }, seed, mc.cores=getOption("mc.cores", 2L))
    options(op)

    # 'x' axis: where summed up the draws inside xiggbin()
    x <- seq(1,steps)*(draws/steps)

    names(res) <- generators
    res <- data.frame(draws=x,
                      as.data.frame(res))

    attr(res, "draws")   <- draws
    attr(res, "bins")    <- 200
    attr(res, "seed")    <- seed
    attr(res, "steps")   <- steps
    attr(res, "created") <- format(Sys.time())
    attr(res, "version") <- packageVersion("RcppZiggurat")

    if (showplot) {
        plotChiSq(res)
    }

    invisible(res)
}
