PoisMixSim <-
function(n = 2000, libsize, separation) {

if(length(libsize) > 1)
	stop(paste(sQuote("libsize"), "must be equal to one of the following:", dQuote("A"), "or",
		dQuote("B"), "or", dQuote("equal"))) 
if(libsize != "A" & libsize != "B" & libsize != "equal")
	stop(paste(sQuote("libsize"), "must be equal to one of the following:", dQuote("A"), "or",
		dQuote("B"), "or", dQuote("equal"))) 
if(length(separation) > 1)
	stop(paste(sQuote("separation"), "must be equal to one of the following:", 
		dQuote("high"), "or", dQuote("low")))
if(separation != "high" & separation != "low")
	stop(paste(sQuote("separation"), "must be equal to one of the following:", 
		dQuote("high"), "or", dQuote("low")))
if(length(n) > 1)
	stop(paste(sQuote("n"), "must be a positive integer"))
if(n <= 0)
	stop(paste(sQuote("n"), "must be a positive integer"))
if(round(n) != n)
	stop(paste(sQuote("n"), "must be a positive integer"))

## libsize <- c("equal", "A", "B")
if(libsize == "A" | libsize == "equal") {
conds <- c(1, rep(2,4), rep(3, 3))
mean.expr <- 1640
s.norm <- c(0.156, 0.071, 0.248, 0.165, 0.014, 0.028, 0.206)
}

if(libsize == "B") {
conds <- c(rep(1, 4), rep(2,2))
mean.expr <- 1521
s.norm <- c(0.096, 0.084, 0.253, 0.205, 0.224, 0.138)
}

r <- table(conds)
cols <- length(conds)
d <- length(unique(conds))
g.true <- 4
w <- round(rexp(n, 1/mean.expr))

s.true <- ifelse(rep(libsize, cols) == "equal", rep(1/cols, cols), 
s.norm)
s.dot.true <- rep(NA, d)
for(j in 1:d) {
s.dot.true[j] <- sum(s.true[which(conds == (unique(conds))[j])])
} 
lambda.true <- matrix(NA, nrow = d, ncol = g.true)

##################################
## CHOOSING LAMBDA VALUES##
##################################

## High separation
if(separation == "high") {
tmp <- cbind(c(1,3,5), c(5,1,3), c(3,5,1), c(5,3,1))
}
## Low separation
if(separation == "low") {
tmp <- cbind(c(1,3,5), c(2,4,4), c(1,5,4), c(2,5,3))
}
if(libsize == "B") tmp <- tmp[1:d,];
## Choosing lambda values so that colSums(s.dot.true * lambda.true) = 1
for(k in 1:g.true) {
lambda.tmp <- tmp[,k]/sum(tmp[,k]);
lambda.true[,k] <- lambda.tmp/s.dot.true
}
pi.true <- c(.10, .20, .30, .40)

################################
## Simulating data##
################################

y <- matrix(NA, nrow = n, ncol = cols)
label.true <- rep(NA, n)
tmp <- runif(n); cp <- cumsum(pi.true);
for(i in 1:n) {
## Choose class label
lab <- 1
for(k in 2:g.true) {
if(tmp[i] < cp[k] & tmp[i] >= cp[k-1]) lab <- k;
}
label.true[i] <- lab
lambda.tmp <- rep(lambda.true[,label.true[i]], times = r)
y[i,] <- rmultinom(1, w[i], s.true*lambda.tmp)
}

## Remove rows with all zeros
if(min(rowSums(y) == 0)) {
y <- y[-which(rowSums(y) == 0),]
label.true <- label.true[-which(rowSums(y) == 0),]
w <- w[-which(rowSums(y) == 0),]
}

return(list(y = y, labels = label.true, pi = pi.true, lambda = lambda.true,
w = w, conditions = conds))
}

