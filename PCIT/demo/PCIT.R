data(PCIT)

m <- m[1:200,1:200]        # just use a small subset of the data

# lets see the size and a subset of the correlation matrix
dim(m)
m[1:10,1:10]

# apply the PCIT algorithm, forcing the serial implementation
system.time(result_serial <- pcit(m, force.serial=TRUE))

# apply the PCIT algorithm, using the parallel implementation if we can, otherwise fall back to the serial implementation
system.time(result <- pcit(m))

# pcit() doesn't return the ind's in the same order when done in parallel and serial
# check that we got the same answer using both functions
all.equal(m[idx(result_serial)], m[idx(result)])

# get the matric indices for the meaningful and unmeaningful correlations
meaningful.idx <- idx(result_serial)
unmeaningful.idx <- idxInvert(nrow(m), meaningful.idx)

# create a copy of the correlation matrix and set unmeaingful correlations to zero
m.new <- m
m.new[unmeaningful.idx] <- 0

# convert adjacency matrix into an edge list
edgeList <- getEdgeList(m.new)

cc <- clusteringCoefficient(m.new)
ccp <- clusteringCoefficientPercent(m.new)
ccp

op <- par(mfrow=c(3,2))
plot(density(m[upper.tri(m)]), main="Density Plot of Raw Correlation Coefficients", xlab="Correlation Coefficient")
hist(cc, main="Connectivity Distribution", xlab="Proportion of Connections", ylab="Number of Genes")
hist(cc*length(cc), main="Connectivity Distribution", xlab="Number of Connections", ylab="Number of Genes")
# plot the distribution of all correlations superimposed by that of the meaningful corrections in black
plotCorCoeff(m, list("PCIT Significant" = meaningful.idx), col=c("black"))
# plot the distribution of all correlations superimposed by that of the meaningful corrections in black and the absolute correlations > 0.5 in red
abs.idx <- which(abs(m) > 0.5)
plotCorCoeff(m, list("abs. cor. > 0.5" = abs.idx, "PCIT Significant" = meaningful.idx), col=c("red", "black"))
# we'll change the order and use some transparent colours using rgb()
plotCorCoeff(m, list("PCIT Significant" = meaningful.idx, "abs. cor. > 0.5" = abs.idx), col=c(rgb(1,0,0,0.7), rgb(0,0,0,0.7)))
par(op)
