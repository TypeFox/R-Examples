variogram <-
function(resid, timeVar, id, binwidth, numElems = 0, inc.var = TRUE, irregular = TRUE){

# process variance
data.proc.var <- data.frame(id = id, res = resid)
idlist        <- unique(data.proc.var$id)
nsubj         <- length(idlist)
cs.nobs       <- cumsum(as.numeric(table(data.proc.var$id)))
res.save      <- c()

for(i in 1 : (nsubj - 1)){
resid1   <- data.proc.var[data.proc.var$id == idlist[i], "res"]
resid2   <- data.proc.var[(cs.nobs[i] + 1) : nrow(data.proc.var), "res"]
result   <- as.numeric(outer(resid1, resid2, function(x, y) (x - y)^2/2))
res.save <- rbind(res.save, c(length(result), sum(result)))
}
proc.var <- colSums(res.save)[2] / colSums(res.save)[1]

# variograms
v  <- unlist(lapply (tapply (resid, id, function (x) outer (x, x, function (x, y) (x - y)^2/2)), function (x) x[lower.tri (x)]))
u  <- unlist(lapply (tapply (timeVar, id, function (x) outer (x, x, function (x, y) abs(x - y))), function (x) x[lower.tri (x)]))


# for irreularly spaced data

if(irregular == TRUE){

bins <- seq(floor(min(u)), ceiling(max(u)), binwidth)

bin.means <- NULL
bin.sizes <- NULL

for (i in 1 : (length(bins) - 1)){

bin.means <- c(bin.means, mean(cbind(u, v)[cbind(u, v)[, 1] >= bins[i] & cbind(u, v)[, 1] < bins[i+1], 2]))
bin.sizes <- c(bin.sizes, length(cbind(u, v)[cbind(u, v)[, 1] >= bins[i] & cbind(u, v)[, 1] < bins[i+1], 2]))

}

bin.mids <- NULL
for (i in 1 : (length(bins) - 1)) bin.mids <- c(bin.mids, (bins[i] + bins[i+1]) / 2)

bin.mids.elim  <- bin.mids[which(bin.sizes > numElems)]
bin.means.elim <- bin.means[which(bin.sizes > numElems)]

if(inc.var == TRUE){
plot(bin.mids.elim, bin.means.elim, type = "l", xlab = "Lag", ylim = c(0, max(proc.var, max(bin.means.elim))+0.01), ylab = "Semivariogram")
lines(bin.mids.elim, rep(proc.var, length(bin.mids.elim)))
}else{
plot(bin.mids.elim, bin.means.elim, type = "l", xlab = "Lag", ylab = "Semivariogram")
}

output              <- list()
output$bin.mids     <- bin.mids
output$bin.means    <- bin.means
output$bin.sizes    <- bin.sizes
output$process.var  <- proc.var

}

if(irregular == FALSE){

lags  <- unique(u)
means <- tapply(v, u, function(x) mean(x))
sizes <- tapply(v, u, function(x) length(x))

if(inc.var == TRUE){
plot(lags, means, type = "l", xlab = "Lag", ylim = c(0, max(proc.var, max(means))+0.01), ylab = "Semivariogram")
lines(lags, rep(proc.var, length(lags)))
}else{
plot(lags, means, type = "l", xlab = "Lag", ylab = "Semivariogram")
}

output <- list()
output$lags  <- lags
output$means <- means
output$sizes <- sizes
output$process.var <- proc.var

}

output

}
