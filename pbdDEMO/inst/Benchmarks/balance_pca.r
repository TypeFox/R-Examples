library(pbdDMAT, quiet = TRUE)
library(rbenchmark)

###################SETTINGS######################

init.grid()

comm.set.seed(1234, diff = TRUE)

# biggest size to check
start <- 100
stop <- 5000
step <- 100

# blocking
bldim <- 2

# normal family
mean <- 100
sd <- 1000

#################################################


# benchmark
sizes <- seq(from=start, to=stop, by=step)
if ( !(stop %in% sizes) ){
	sizes <- c(sizes, stop)
}

for (N in sizes){
	if (comm.rank()==0){
		x <- matrix(rnorm(N^2, mean=mean, sd=sd), nrow=N, ncol=N)
		srl <- system.time(prcomp(x))[3]
	}
	else {
		x <- NULL
		srl <- 0
	}
	srl <- allreduce(srl)
	barrier()
	
	dx <- as.ddmatrix(x)
	prl <- system.time(prcomp(dx))[3]
	prl <- allreduce(prl, op='max')

	if (prl < srl){
		break
	}
}

size <- N*N*8/1024
unit <- "kb"
if (log10(size) > 3){
	size <- size/1024
	unit <- "mb"
}
if (log10(size) > 3){
	size <- size/1024
	unit <- "gb"
}

comm.cat(sprintf("\n############## prcomp(x) ##############\n", comm.size()), quiet=T)
comm.cat(sprintf("dim(x):\t\t%dx%d ~ %.2f %s\n", N, N, size, unit), quiet=T)
comm.cat(sprintf("1 core R:\t%.3f seconds\n", srl), quiet=T)
comm.cat(sprintf("%d core pbdR:\t%.3f seconds\n", comm.size(), prl), quiet=T)

finalize()
