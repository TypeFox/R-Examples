library(pbdDMAT, quiet = TRUE)
library(rbenchmark)

###################SETTINGS######################

init.grid()

comm.set.seed(1234, diff = TRUE)

# biggest size to check
ncol <- 100
start <- 100
stop <- 5000
step <- 50

# blocking
bldim <- 4

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
		x <- matrix(rnorm(N*ncol, mean=mean, sd=sd), nrow=N, ncol=ncol)
		y <- matrix(rnorm(N, mean=mean, sd=sd), nrow=N, ncol=1)
		srl <- system.time(lm.fit(x=x, y=y))[3]
	}
	else {
		x <- NULL
		y <- NULL
		srl <- 0
	}
	srl <- allreduce(srl)
	barrier()
	
	dx <- as.ddmatrix(x, bldim=bldim)
	dy <- as.ddmatrix(y, bldim=bldim)
	prl <- system.time(lm.fit(x=dx, y=dy))[3]
	prl <- allreduce(prl, op='max')

	if (prl < srl){
		break
	}
}

size <- N*ncol*8/1024
unit <- "kb"
if (log10(size) > 3){
	size <- size/1024
	unit <- "mb"
}
if (log10(size) > 3){
	size <- size/1024
	unit <- "gb"
}


comm.cat(sprintf("\n############## lm.fit(x, y) ##############\n", comm.size()), quiet=T)
comm.cat(sprintf("dim(x):\t\t%dx%d ~ %.2f %s\n", N, ncol, size, unit), quiet=T)
comm.cat(sprintf("1 core R:\t%.3f seconds\n", srl), quiet=T)
comm.cat(sprintf("%d core pbdR:\t%.3f seconds\n", comm.size(), prl), quiet=T)

finalize()
