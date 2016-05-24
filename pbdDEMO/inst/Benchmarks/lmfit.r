library(pbdDMAT, quiet = TRUE)

###################SETTINGS######################
init.grid()

comm.set.seed(1234, diff = TRUE)

# size
N <- 15000
p <- 500

# blocking
bldim <- 4

# normal family
mean <- 100
sd <- 1000

# replications
reps <- 10

#################################################


# benchmark
datatimes <- system.time({
	dx <- ddmatrix("rnorm", nrow=N, ncol=p, bldim=bldim, mean=mean, sd=sd, ICTXT=0)
	dy <- ddmatrix("rnorm", nrow=N, ncol=1, bldim=bldim, mean=mean, sd=sd, ICTXT=0)
})[3]

datatimes <- allreduce(datatimes, op='max')

size <- N*p*8/1024
unit <- "kb"
if (log10(size) > 3){
	size <- size/1024
	unit <- "mb"
}
if (log10(size) > 3){
	size <- size/1024
	unit <- "gb"
}

comm.cat(sprintf("\n%.2f %s of data generated in %.3f seconds\n\n", size, unit, datatimes), quiet=T)


times <- sapply(1:reps, function(.) system.time(lm.fit(x=dx, y=dy))[3])
total <- allreduce(sum(times), op='max')
avg <- total/reps

bench <- data.frame(operation="lm.fit(dx, dy)", mean.runtime=avg, total.runtime=total)
row.names(bench) <- ""
comm.print(bench, quiet=T)

finalize()
