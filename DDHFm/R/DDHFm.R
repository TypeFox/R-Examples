"DDHFm" <-
function(x){
#
# Apply data driven Haar-Fisz to arbitrary sized set
#
# Genes are rows, replicates are columns
#
# First compute row means
#
rowmean <- apply(x, 1, mean)
#
# Work out how to sort rows (and how to sort them back again)
#
sl <- sort.list(rowmean)
usl <- sort.list(sl)
#
# Now sort rows of data matrix according to increasing row mean
#
xs <- x[sl,]
#
# Convert matrix to a vector in prep. for HF processing
#
xsv <- as.vector(t(xs))
lxsv <- length(xsv)
#
# Work out padding required.
#
J <- ceiling(log(lxsv)/log(2))
lpx <- 2^J
lpadx <- floor((lpx - lxsv)/2)
rpadx <- lpx - lxsv - lpadx
#
# Do the padding
#
padx <- c( rep(0,lpadx), xsv, rep(0, rpadx))
#
# Do the DDHFT
#
padhf <- ddhft.np.2(padx)$hft
#
# Do the unpadding
#
hf <- padhf[ (lpadx+1):(lpadx+lxsv)]
#
# Put back into matrix form 
#
xshf <- matrix(hf, nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
#
# Unsort rows
#
xshf[usl,]
}

