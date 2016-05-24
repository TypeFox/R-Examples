#!/usr/bin/r -t

suppressMessages(library(gcbd))
suppressMessages(library(getopt))

## default values
N <- 250				# nobs
n <- 30 				# runs
verbose <- FALSE			# verbosity
benchmark <- "matmult"			# benchmark to run

## option specification
options <- matrix(c('verbose',   'v', 0, "integer",   paste("verbose operations, default is", ifelse(verbose, "true", "false")),
                    'help',      'h', 0, "logical",   paste("help on options"),
                    'nobs',      'n', 1, "integer",   paste("number of rows and columns in matrix, default is", N),
                    'runs',      'r', 1, "integer",   paste("number of benchmark runs, default is", n),
                    'benchmark', 'b', 1, "character", paste("benchmark to run (matmult, qr, svd, lu), default is", benchmark)),
                    ncol=5,
                    byrow=TRUE)

## parse argv
opt <- getopt(options, opt=argv)

if (!is.null(opt$help)) {
    cat(getopt(options, opt=argv, command="benchmark.r", usage=TRUE))
    q(status=100)
}

N         <- ifelse(is.null(opt$nobs), N, opt$nobs)
n         <- ifelse(is.null(opt$runs), n, opt$runs)
verbose   <- ifelse(is.null(opt$verbose), verbose, opt$verbose)
benchmark <- ifelse(is.null(opt$benchmark), benchmark, opt$benchmark)

benchmark <- match.arg(benchmark, c("matmult", "qr", "svd", "lu"))

if (verbose) {
    cat("Running", benchmark, "on dimension", N, "for", n, "repeats\n")
}

requirements()

purgeAtlas()
purgeAtlas39()
purgeMKL()
purgeGoto()

cmd      <- paste('r -lgcbd,Matrix   -e"cat(', benchmark, 'Benchmark(', N, ",", n, '))"', sep="")
#cmdmagma <- paste('r -lgcbd,magma    -e"cat(', benchmark, 'Benchmarkmagma(', N, ",", n, '))"', sep="")
cmdgpu   <- paste('r -lgcbd,gputools -e"cat(', benchmark, 'Benchmarkgputools(', N, ",", n, '))"', sep="")
cmdtest  <- paste('r -lgcbd,Matrix   -e"cat(', benchmark, 'BenchmarkTEST(', N, ",", n, '))"', sep="")

ref <- as.numeric(system(cmd, intern=TRUE))

installAtlas()
atlas      <- as.numeric(system(cmd,      intern=TRUE))
purgeAtlas()

installAtlas39()
atl39      <- as.numeric(system(cmd,      intern=TRUE))
purgeAtlas39()

installMKL()
mkl      <- as.numeric(system(cmd,      intern=TRUE))
purgeMKL()

installGoto()
goto      <- as.numeric(system(cmd,      intern=TRUE))
#magmaGoto <- ifelse( hasMagma(), as.numeric(system(cmdmagma, intern=TRUE)), NA)
purgeGoto()

gpu <- ifelse( hasGputools(), as.numeric(system(cmdgpu, intern=TRUE)), NA)

res <- data.frame(type=benchmark, nobs=N, nrun=n,
                  ref=ref, atlas=atlas, atl39=atl39,
                  mkl=mkl, gotob=goto, gpu=gpu)

if (verbose) print(res)
databaseResult(res)
