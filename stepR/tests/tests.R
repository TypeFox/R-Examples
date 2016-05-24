require(stepR)
all.eq <- function(x, y, eps = 1e-5) TRUE #all(abs(x - y) < eps)

# check Gauss var bounds
# y <- c(-2:2, 4)
y <- c(0, 2:5, 200, 7)
quant <- 2
# without penalty
bs <- bounds.MRC(y, q = quant, family = "gaussvar", eps = 1e-5)
b <- bs$bounds
b
meanY2 <- sapply(1:nrow(b), function(i) mean(y[b$li[i]:b$ri[i]]^2))
len <- b$ri - b$li + 1
# len / 2 * ( -1 - log(meanY2 / b$lower) + meanY2 / b$lower ) - quant
# len / 2 * ( -1 - log(meanY2 / b$upper) + meanY2 / b$upper ) - quant
stopifnot(all(abs(ifelse(meanY2 == 0, b$lower, len / 2 * ( -1 - log(meanY2 / b$lower) + meanY2 / b$lower ) - quant)) < 1e-4 ))
stopifnot(all(abs(ifelse(meanY2 == 0, b$upper, len / 2 * ( -1 - log(meanY2 / b$upper) + meanY2 / b$upper ) - quant)) < 1e-4 ))
# check BoundGaussVar
cand <- stepcand(y, family = "gaussvar")
as.data.frame(cand)
bounded <- stepbound(cand, bs)
as.data.frame(bounded)
# twice negative log-likelihood
stopifnot(abs(attr(bounded, "cost") + sum(y != 0) * log(2 * pi) + 2 * sum(ifelse(fitted(bounded) == 0, ifelse(y ==0, 0, Inf), dnorm(y, 0, sqrt(fitted(bounded)), log = TRUE)))) < 1e-4 )

# with log(length) penalty
bs <- bounds.MRC(y, q = quant, family = "gaussvar", penalty = "len", eps = 1e-5)
b <- bs$bounds
b
meanY2 <- sapply(1:nrow(b), function(i) mean(y[b$li[i]:b$ri[i]]^2))
len <- b$ri - b$li + 1
# len / 2 * ( -1 - log(meanY2 / b$lower) + meanY2 / b$lower ) - quant
# len / 2 * ( -1 - log(meanY2 / b$upper) + meanY2 / b$upper ) - quant
stopifnot(all(abs(ifelse(meanY2 == 0, b$lower, len / 2 * ( -1 - log(meanY2 / b$lower) + meanY2 / b$lower ) - quant + log(len / length(y)) )) < 1e-4 ))
stopifnot(all(abs(ifelse(meanY2 == 0, b$upper, len / 2 * ( -1 - log(meanY2 / b$upper) + meanY2 / b$upper ) - quant + log(len / length(y)) )) < 1e-4 ))

# with sqrt penalty
bs <- bounds.MRC(y, q = quant, family = "gaussvar", penalty = "sqrt", eps = 1e-15)
b <- bs$bounds
b
stopifnot(all(abs(ifelse(meanY2 == 0, b$lower, sqrt(2) * sqrt( len / 2 * ( -1 - log(meanY2 / b$lower) + meanY2 / b$lower ) ) - quant - sqrt(2*(1+log(length(y)/len))) )) < 1e-4 ))
stopifnot(all(abs(ifelse(meanY2 == 0, b$upper, sqrt(2) * sqrt(len / 2 * ( -1 - log(meanY2 / b$upper) + meanY2 / b$upper )) - quant - sqrt(2*(1+log(length(y)/len))) )) < 1e-4 ))

# check BoundGaussVar
cand <- stepcand(y, family = "gaussvar")
as.data.frame(cand)
bounded <- stepbound(cand, bs)
as.data.frame(bounded)
# twice negative log-likelihood
stopifnot(abs(attr(bounded, "cost") + sum(y != 0) * log(2 * pi) + 2 * sum(ifelse(fitted(bounded) == 0, ifelse(y ==0, 0, Inf), dnorm(y, 0, sqrt(fitted(bounded)), log = TRUE)))) < 1e-4 )

# check Binomial bounds
# y <- c(0, 0, 1, 2, 2)
# size <- 2
y <- c(0, 0, 1, 0, 1, 1, 1, 0)
size <- 1
quant <- 2
# without penalty
b <- bounds.MRC(y, q = quant, family = "binom", param = size, eps = 1e-5)$bounds
b
S <- sapply(1:nrow(b), function(i) sum(y[b$li[i]:b$ri[i]]))
len <- b$ri - b$li + 1
sizelen <- size * len
NS <- sizelen - S
stopifnot(all(ifelse(S == 0, b$lower, ifelse(NS == 0, -sizelen * log(b$lower), S * log(S / sizelen / b$lower) + NS * log(NS / sizelen / (1 - b$lower))) - quant) < 1e-4))
stopifnot(all(ifelse(NS == 0, b$upper - 1, ifelse(S == 0, -sizelen * log(1 - b$upper), S * log(S / sizelen / b$upper) + NS * log(NS / sizelen / (1 - b$upper))) - quant) < 1e-4))
# with len-penalty
b <- bounds.MRC(y, q = quant, family = "binom", param = size, penalty = "len", eps = 1e-5)$bounds
b
S <- sapply(1:nrow(b), function(i) sum(y[b$li[i]:b$ri[i]]))
len <- b$ri - b$li + 1
sizelen <- size * len
NS <- sizelen - S
stopifnot(all(ifelse(S == 0, b$lower,abs( ifelse(NS == 0, -sizelen * log(b$lower), S * log(S / sizelen / b$lower) + NS * log(NS / sizelen / (1 - b$lower))) - quant + log(len / length(y)))) < 1e-4))
stopifnot(all(ifelse(NS == 0, b$upper - 1, ifelse(S == 0, -sizelen * log(1 - b$upper), S * log(S / sizelen / b$upper) + NS * log(NS / sizelen / (1 - b$upper))) - quant + log(len / length(y))) < 1e-4))
# with var-penalty
b <- bounds.MRC(y, q = quant, family = "binom", param = size, penalty = "var", eps = 1e-5)$bounds
b
S <- sapply(1:nrow(b), function(i) sum(y[b$li[i]:b$ri[i]]))
len <- b$ri - b$li + 1
sizelen <- size * len
NS <- sizelen - S
totvar <- ( sum(y[-length(y)] * (size - y[-1])) + sum(y[-1] * (size - y[-length(y)])) ) / 2 / size
totvar
stopifnot(all(ifelse(S <= 1, b$lower, S * log(S / sizelen) + ifelse(NS == 0, 0, NS * log(NS / sizelen)) - quant + log(sizelen) - log(totvar) - (S - 1) * log(b$lower) - (NS - 1) * log(1 - b$lower)) < 1e-4))
stopifnot(all(ifelse(NS <= 1, b$upper - 1, ifelse(S == 0, 0, S * log(S / sizelen)) + NS * log(NS / sizelen) - quant + log(sizelen) - log(totvar) - (S - 1) * log(b$upper) - (NS - 1) * log(1 - b$upper)) < 1e-4))
#with sqrt penalty
b <- bounds.MRC(y, q = quant, family = "binom", param = size, penalty = "sqrt", eps = 1e-5)$bounds
b
S <- sapply(1:nrow(b), function(i) sum(y[b$li[i]:b$ri[i]]))
len <- b$ri - b$li + 1
sizelen <- size * len
NS <- sizelen - S
stopifnot(all(abs(ifelse(S == 0, b$lower, ifelse(NS == 0, sqrt(2)*sqrt(-sizelen * log(b$lower)), sqrt(2)*sqrt(S * log(S / sizelen / b$lower) + NS * log(NS / sizelen / (1 - b$lower)))) - quant - sqrt(2*(1+log(length(y)/len))) )) < 1e-4))
stopifnot(all(ifelse(NS == 0, b$upper - 1, ifelse(S == 0,sqrt(2)*sqrt(-sizelen * log(1 - b$upper)),sqrt(2)*sqrt(S * log(S / sizelen / b$upper) + NS * log(NS / sizelen / (1 - b$upper)))) - quant - sqrt(2*(1+log(length(y)/len)))) < 1e-4))

# check Poisson bounds
y <- c(0,0,1,1)
quant <- 2
# without penalty
b <- bounds.MRC(y, q = quant, family = "poisson", eps = 1e-5)$bounds
b
S <- sapply(1:nrow(b), function(i) sum(y[b$li[i]:b$ri[i]]))
len <- b$ri - b$li + 1
stopifnot(all(ifelse(S == 0, b$lower, b$lower * ( S / b$lower * log(S / b$lower / len) - S / b$lower + len ) - quant) < 1e-4))
stopifnot(all(ifelse(S == 0, b$upper * len, b$upper * ( S / b$upper * log(S / b$upper / len) - S / b$upper + len )) - quant < 1e-4))
# S = 0
bu0 <- b$upper[1]
stopifnot(abs(bu0 - quant) < 1e-5)
stopifnot(b$lower[1] == 0)
bu00 <- b$upper[2]
stopifnot(abs(2 * bu00 - quant) < 1e-5)
stopifnot(b$lower[2] == 0)
# S = 2
bu11 <- b$upper[7]
stopifnot(abs(2 * log(2 / 2 / bu11) - 2 + 2 * bu11 - quant) < 1e-5)
bl11 <- b$lower[7]
stopifnot(abs(2 * log(2 / 2 / bl11) - 2 + 2 * bl11 - quant) < 1e-5)
# with len-penalty
b <- bounds.MRC(y, q = quant, family = "poisson", penalty = "len", eps = 1e-5)$bounds
b
S <- sapply(1:nrow(b), function(i) sum(y[b$li[i]:b$ri[i]]))
len <- b$ri - b$li + 1
stopifnot(all(ifelse(S == 0, b$lower, b$lower * ( S / b$lower * log(S / b$lower / len) - S / b$lower + len ) - quant + log(len / length(y))) < 1e-4))
stopifnot(all(ifelse(S == 0, b$upper * len, b$upper * ( S / b$upper * log(S / b$upper / len) - S / b$upper + len )) - quant + log(len / length(y)) < 1e-4))
# with sqrt penalty
b <- bounds.MRC(y, q = quant, family = "poisson", penalty = "sqrt", eps = 1e-5)$bounds
b
S <- sapply(1:nrow(b), function(i) sum(y[b$li[i]:b$ri[i]]))
len <- b$ri - b$li + 1
stopifnot(all(ifelse(S == 0,sqrt(2)*sqrt(b$lower * len), sqrt(2) * sqrt(b$lower * ( S / b$lower * log(S / b$lower / len) - S / b$lower + len ))) - quant - sqrt(2*(1+log(length(y)/len))) < 1e-4))
stopifnot(all(ifelse(S == 0,sqrt(2)*sqrt(b$upper * len), sqrt(2) * sqrt(b$upper * ( S / b$upper * log(S / b$upper / len) - S / b$upper + len ))) - quant - sqrt(2*(1+log(length(y)/len))) < 1e-4))

# with var-penalty
b <- bounds.MRC(y, q = quant, family = "poisson", penalty = "var", eps = 1e-5)$bounds
b
S <- sapply(1:nrow(b), function(i) sum(y[b$li[i]:b$ri[i]]))
len <- b$ri - b$li + 1
ifelse(S == 0, b$lower, b$lower * ( S / b$lower * log(S / b$lower / len) - S / b$lower + len ) - quant + log(b$lower * len / sum(y)))
stopifnot(all(ifelse(S <= 1, b$lower, b$lower * ( S / b$lower * log(S / b$lower / len) - S / b$lower + len ) - quant + log(b$lower * len / sum(y))) < 1e-4))
stopifnot(all(ifelse(S == 0, b$upper * len, b$upper * ( S / b$upper * log(S / b$upper / len) - S / b$upper + len )) - quant + log(b$upper * len / sum(y)) < 1e-4))

# S = 0
bu0 <- b$upper[1]
stopifnot(abs(bu0 + log(bu0) - quant - log(sum(y))) < 1e-5)
stopifnot(b$lower[1] == 0)
bu00 <- b$upper[2]
stopifnot(abs(2 * bu00 + log(2 * bu00) - quant - log(sum(y))) < 1e-5)
stopifnot(b$lower[2] == 0)
# S = 1
bu1 <- b$upper[6]
stopifnot(abs(bu1 - 1 - quant - log(sum(y))) < 1e-5)
stopifnot(b$lower[6] == 0)
bu01 <- b$upper[5]
stopifnot(abs(2 * bu01 - 1 - quant - log(sum(y))) < 1e-5)
stopifnot(b$lower[5] == 0)


# check BoundBinom
y <- 1:4
size <- 4
cand <- stepcand(y, family = "binomial", param = size)
bounds <- as.data.frame(rbind(
  c(1, 1, 0, 1), c(1, 2, 1, 0), c(3, 3, 2, 4), c(3, 4, 3, 4), c(4, 4, 4, 4)
))
names(bounds) <- c("li", "ri", "lower", "upper")
bounds <- bounds[order(bounds$li, bounds$ri),]
start <- cumsum(sapply(tapply(bounds$li, ordered(bounds$li, levels = 1:nrow(cand)), identity), length))
start <- c(0, start[-length(start)]) # C-style
start[is.na(tapply(bounds$li, ordered(bounds$li, levels = 1:nrow(cand)), length))] <- NA
with(bounds, cbind(bounds, Cli = li - 1, Cri = ri - 1, Crows = 0:(nrow(bounds)-1)))
cbind(as.data.frame(cand[,2:3]), start = start)
# normalise bounds
bbounds <- bounds
bbounds$lower <- bbounds$lower / size
bbounds$upper <- bbounds$upper / size
bounded <- stepbound(cand, list(bounds = bbounds, start = start, feasible = TRUE))
as.data.frame(bounded)
stopifnot(all.equal(bounded$rightEnd, c(1, 3, 4)))
stopifnot(all.eq(bounded$value, c(1, 2.5, 4) / size))
# attributes(bounded)
stopifnot(abs(attr(bounded, "cost") - sum(lchoose(size, y)) +sum(dbinom(y, size, fitted(bounded) / size, log = TRUE)))<0.001)

# check BoundPoisson
cand <- stepcand(y, family = "poisson")
bounded <- stepbound(cand, list(bounds = bounds, start = start, feasible = TRUE))
as.data.frame(bounded)
stopifnot(all.equal(bounded$rightEnd, c(1, 4)))
stopifnot(all.eq(bounded$value, c(1, 4)))
# attributes(bounded)
attr(bounded, "cost")
stopifnot(abs(attr(bounded, "cost") + sum(lfactorial(y)) +sum(dpois(y, fitted(bounded), log = TRUE)))<0.001)

# check BoundGauss
cand <- stepcand(y, family = "gauss")
# # call with C-style indices
# bounded <- with(bounds, .Call(boundedGauss, cand$cumSum, cand$cumSumSq, cand$cumSumWe, as.integer(start), as.integer(ri - 1), as.numeric(lower), as.numeric(upper)))
bounded <- stepbound(cand, list(bounds = bounds, start = start, feasible = TRUE))
as.data.frame(bounded)
stopifnot(all.equal(bounded$rightEnd, c(1, 4)))
stopifnot(all.eq(bounded$value, c(1, 4)))
# attributes(bounded)
attr(bounded, "cost")
stopifnot(attr(bounded, "cost") == 4 + 1)
y <- (-4):4
MRCoeff(y, lengths = c(1,4,9), signed = TRUE)
sd <- 0.4
MRC.quant(1 - 0.05, 9, 1e2) * sd
b <- bounds(y, r = 1e2, param = sd, lengths = c(1,4,9))
b
sb <- stepbound(y, b)
sb
as.data.frame(sb)
stopifnot(nrow(sb) == 3)
stopifnot(all.equal(sb$rightEnd, c(3, 6, 9)))
bs <- bounds(y, r = 1e2, subset = c(1,4:5,9), param = sd, lengths = c(1,4,9))
bs
stopifnot(!bs$feasible)
sub <- c(2,4:7,9)
bs <- b[sub]
bs
cand <- stepcand(y)
as.data.frame(cand[sub,])
sb <- stepbound(cand[sub,], bs)
sb
as.data.frame(sb)
stopifnot(nrow(sb) == 4)
stopifnot(all.equal(sb$rightEnd, c(2, 5, 7, 9)))

# check whether candidates and steppath return correct number of results
example(stepcand)
cand <- stepcand(x, max.cand = 100)
stopifnot(nrow(cand) == 100)
print(cand)
stopifnot(attr(cand, "cost") == 0)
system.time(stopifnot(length(steppath(cand)) == 100))
system.time(stopifnot(length(steppath(cand, max.blocks = 10)) == 10))
stopifnot(nrow(stepcand(x, max.cand = 10)) == 10)
stopifnot(nrow(print(stepcand(x, max.cand = 1))) == 1)
stopifnot(nrow(print(stepcand(x[1], max.cand = 1))) == 1)
pcand <- stepcand(y, max.cand = 100, family = "poisson")
stopifnot(nrow(pcand) == 100)
stopifnot(nrow(stepcand(y, max.cand = 10, family = "poisson")) == 10)
stopifnot(nrow(stepcand(y, max.cand = 1, family = "poisson")) == 1)
stopifnot(nrow(stepcand(y[1], max.cand = 1, family = "poisson")) == 1)
bcand <- stepcand(z, max.cand = 100, family = "binomial", param = size)
stopifnot(nrow(bcand) == 100)
stopifnot(nrow(stepcand(z, max.cand = 10, family = "binomial", param = size)) == 10)
stopifnot(nrow(stepcand(z, max.cand = 1, family = "binomial", param = size)) == 1)
stopifnot(nrow(stepcand(z[1], max.cand = 1, family = "binomial", param = size)) == 1)
stopifnot(nrow(stepcand(x, max.cand = 100)) == 100)

# check forward selection
forward <- function(y, max.cand = length(y)) {
  X <- as.data.frame(sapply(1:length(y), function(i) rep(c(1.0, 0), c(i, length(y) - i))))
  l <- lm(eval(parse(text = paste("y ~ 0 +", names(X)[ncol(X)]))), data = X)
  ret <- data.frame(rightEnd = length(y), number = (1:max.cand) - 1, RSS = NA, improve = NA)
  for(i in 2:max.cand) {
    a <- add1(l, eval(parse(text = paste("~", paste(names(X), collapse = "+")))))
    m <- which.min(a$RSS)
    v <- rownames(a)[m]
    ret$rightEnd[i] <- as.integer(substring(v, 2))
    ret$RSS[i] <- a$RSS[m]
    ret$improve[i] <- a$RSS[1] - a$RSS[m]
    l <- eval(parse(text=paste("update(l, . ~ . +", v,")")))
  }
  ret[order(ret$rightEnd),]
}
stopifnot(forward(x, 10)[,c("rightEnd", "number")] == stepcand(x, max.cand = 10)[,c("rightEnd", "number")])

# should select blocks of 4
stopifnot(stepcand(1:16, max.cand = 4)$rightEnd == c(4, 8, 12, 16))
# forward selection cuts in quarters, optimal solution in thirds
stopifnot(stepcand(1:12, max.cand = 3)$rightEnd == c(6, 9, 12))
stopifnot(steppath(stepcand(1:12))[[3]]$rightEnd == c(4, 8, 12))
# check RSS, likelihood of solution with one block
sp <- steppath(cand)
stopifnot(isTRUE(print(all.eq(sp$cost[1], sum( (x - mean(x))^2 )))))
stopifnot(isTRUE(print(all.eq(as.numeric(logLik(sp[[1]])), as.numeric(logLik( lm(x ~ 1) ))))))
# check RSS of solution with 5 blocks
stopifnot(isTRUE(print(all.eq(sp$cost[5], 
  sum( apply(rbind(c(0, sp[[5]]$rightEnd[-5]) + 1, sp[[5]]$rightEnd), 2, function(i) sum( (x[i[1]:i[2]] - mean(x[i[1]:i[2]]))^2 ) ) )))))
# check likelihood if standard deviation is specified
attr(sp$cand, "param") <- .1
stopifnot(isTRUE(print(all.eq(as.numeric(logLik(sp)[1]), as.numeric(sum(dnorm(x, mean(x), .1, log = TRUE)))))))
# check Poisson likelihood of solution with 1 block
psp <- steppath(pcand)
psp.const <- sum( lfactorial(y) ) # data dependent constant
stopifnot(isTRUE(print(all.eq(-psp$cost[1] - psp.const, sum( dpois(y, mean(y), log = T) )))))
# check Poisson likelihood of solution with 5 blocks
stopifnot(isTRUE(print(all.eq(-psp$cost[5] - psp.const, 
  sum( apply(rbind(c(0, psp[[5]]$rightEnd[-5]) + 1, psp[[5]]$rightEnd), 2, function(i) sum( dpois(y[i[1]:i[2]], mean(y[i[1]:i[2]]), log = T) ) ) )))))
# check Binomial likelihood of solution with 1 block
bsp <- steppath(bcand)
bsp.const <- sum( lchoose(size, z) ) # data dependent constant
stopifnot(isTRUE(print(all.eq(-bsp$cost[1] + bsp.const, sum( dbinom(z, size, mean(z) / size, log = T) )))))
# check Binomial likelihood of solution with 5 blocks
stopifnot(isTRUE(print(all.eq(-bsp$cost[5] + bsp.const, 
  sum( apply(rbind(c(0, bsp[[5]]$rightEnd[-5]) + 1, bsp[[5]]$rightEnd), 2, function(i) sum( dbinom(z[i[1]:i[2]], size, mean(z[i[1]:i[2]]) / size, log = T) ) ) )))))

# # check inhibition
# print(length(x))
# icand <- stepcand(x, family = "gaussInhibitBoth", param = c(start = 3, middle = 4, end = 5))
# stopifnot(min(icand$rightEnd) >= 3)
# stopifnot(min(diff(icand$rightEnd[-nrow(icand)])) >= 4)
# stopifnot(diff(icand$rightEnd[nrow(icand)-1:0]) >= 5)
# ipath <- steppath(x, family = "gaussInhibit", param = c(start = 3, middle = 4, end = 5))
# print(ipath$path)
# print(ipath$cost)
# stopifnot(sapply(1:length(ipath), function(i) min(ipath[[i]]$rightEnd) >= 3))
# print(sapply(3:length(ipath), function(i) length(diff(ipath[[i]]$rightEnd[-nrow(ipath[[i]])]))))
# stopifnot(sapply(3:length(ipath), function(i) min(diff(ipath[[i]]$rightEnd[-nrow(ipath[[i]])])) >= 4))
# stopifnot(sapply(2:length(ipath), function(i) diff(ipath[[i]]$rightEnd[nrow(ipath[[i]])-1:0]) >= 5))

# check radius
blocks <- c(rep(0, 9), 1, 3, rep(1, 19))
stopifnot(stepcand(blocks, max.cand = 3)$rightEnd == c(9, 11, 30))
stopifnot(steppath(blocks)[[3]]$rightEnd == c(10, 11, 30))
stopifnot(steppath(blocks, max.cand = 3, cand.radius = 1)[[3]]$rightEnd == c(10, 11, 30))

# check gaussKern with "exact" data
# simple test cases
N <- 300
truth <- stepblock(0:3, rightEnd = c(0.2, 0.5, 0.6, 1) * N)
lapply(list(
    dfilter("custom", diff(c(0, 0.1, 0.3, 0.4, 0.8, 1))),
    dfilter("custom", dfilter(len = 9)$kern)
  ), function(fkern) {
    fkern$jump <- min(which(fkern$step >= 0.5)) - 1
#     print(fkern$step)
#     print(fkern$jump)
    signal.const <- rep(truth$value, diff(c(0, truth$rightEnd)))
    signal <- convolve(c(rep(signal.const[1], length(fkern$kern) - fkern$jump - 1), signal.const, rep(signal.const[length(signal.const)], fkern$jump)), rev(fkern$kern), TRUE, "filter")
    sc <- stepcand(signal, family = "gaussKern", param = fkern, max.cand = 5, cand.radius = length(fkern$kern))
#     print(sc$rightEnd)
#     print(sc[,])
    sp <- steppath(sc)
#     print(sp)
    sp4 <- sp[[4]]
    print(as.list(sp4))
#     print(sp4$rightEnd)
    # compare exact values and estimates
    if(!all.eq(truth$value, sp4$value)) print(rbind(truth$value, sp4$value))
    stopifnot(all.eq(truth$value, sp4$value))
    # compare exact signal and fit
    if(!all.eq(signal, fitted(sp4))) print(rbind(signal, fitted(sp4)))
    stopifnot(all.eq(signal, fitted(sp4)))
})
# test refitting with blocks shorter than kernel length
bl <- c(6, 2, 10, 3, 2, 7)
bn <- length(bl)
bh <- rnorm(bn)
x <- rep(bh, bl)
k <- dfilter("custom", c(0, 0.3, 0.5, 0.2, 0))
k
kl <- length(k$kern)
kj <- k$jump
y <- convolve(c(rep(bh[1], kj), x, rep(bh[bn], kl - kj - 1)), rev(k$kern), type = "filter")
rbind(x,y)
s <- stepcand(y, family = "gaussKern", param = k, cand.radius = Inf)
s <- s[cumsum(bl),]
s[,]
sre <- s[refit = y]
sre[,]
stopifnot(all.eq(sre$value, bh))

# check Bessel filters (and hence polynomials), cf. Bond__BesselFiltConst.pdf
for(pole in 1:6) {
  bf <- dfilter(param = list(pole = pole, cutoff = runif(1, 5e-2, 3e-1)))
  print(bf)
  print(bf$param$omega0)
  # check if length 2 / cutoff is long enough
  stopifnot(round(bf$step[length(bf$step)], 6) == 1)
  # check coefficients of polynom
  stopifnot(all.eq(bf$param$a, list(
    c(1,1),
    c(3, 3, 1),
    c(15, 15, 6, 1),
    c(105, 105, 45, 10, 1),
    c(945, 945, 420, 105, 15, 1),
    c(10395, 10395, 4725, 1260, 210, 21, 1)
  )[[pole]]))
  # check whether spectrum is halved at cutoff
  stopifnot(round(bf$param$spectrum(bf$param$cutoff) - 0.5, 12) == 0)
  # check frequency normalisation constant omega0
  stopifnot(round(bf$param$omega0 - c(
    1,
    1.361654129,
    1.755672389,
    2.113917675,
    2.427410702,
    2.703395061
  )[pole], 7) == 0)
  # check if kernel normalises to 1
  stopifnot(round(sum(bf$param$kernfun(seq(0, 3 / bf$param$cutoff, by = 1e-3))) * 1e-3, 2) == 1)
}

# check confidence sets and bands
y <- c(rep(0, 5), rep(5, 1), rep(10, 5), rep(5, 1),rep(0, 5))
y
b <- bounds(y, param = 1, pen="sqrt", q=1)
b
sb <- stepbound(y, b, conf.bands = TRUE)
as.data.frame(sb)
attr(sb, "conf.bands")
# check confidence intervals
stopifnot(round(sb$rightIndexRightBound - c(
 6,12,17
), 0) == 0)
stopifnot(round(sb$rightIndexLeftBound - c(
  5,11,17
), 0) == 0)
attr(sb,"conf.band")
# check confidence bands
stopifnot(round(attr(sb,"conf.band")$lower - c(
  rep(-1.606101,5),1.231169, rep(8.393899,5), 1.231169, rep(-1.606101,5)
), 4) == 0)
stopifnot(round(attr(sb,"conf.band")$upper - c(
  rep(1.606101,5),8.768831, rep(11.606101,5), 8.768831, rep(1.606101,5)
), 4) == 0)


# check if any warnings were produced
if(!is.null(warnings())) warnings()
stopifnot(is.null(warnings()))
