# A. Estimate the static AIDS model; show the elasticity estimates
setwd("C:/aErer"); source("r121wanJAAE.r", echo = FALSE)
hSta; names(hSta); summary(hSta)
es <- aiElas(z = hSta); names(es); es 

# B. Inputs for elasticity computation
z <- hSta; cof <- coef(z$est); vco <- vcov(z$est)
nE <- z$nExoge; nS <- z$nShare; nP <- z$nParam
av <- colMeans(z$y[, z$share])  # omitted variable may not be the last one
av.sh <- av[c(z$share[-z$nOmit], z$share[z$nOmit])]

# C1. Expenditure elasticity: estimate, variance, t, p
e.ex <- v.ex <- NULL
for(i in 1:(nS - 1)) {
  loc.beta <- nE + nP * (i - 1)  # beta location
  e.ex[i] <- 1 + cof[loc.beta] / av.sh[i]
  v.ex[i] <- vco[loc.beta, loc.beta] / (av.sh[i] ^ 2)
}
t.ex <- e.ex / sqrt(v.ex)
p.ex <- 2 * (1 - pt(q = abs(t.ex), df = df.residual(z$est)))

# C2. Expenditure elasticity: combined
ex <- data.frame(cbind(e.ex, sqrt(v.ex), t.ex, p.ex))
rownames(ex) <- names(av.sh)[-z$nOmit]
expen <- bsTab(ex)
colnames(expen) <- c("Elas.expen", "Estimate")
expen

# D1. Hicksian elasticiity: estimate, variance, t, p
e.hi <- v.hi <- matrix(data = NA, nrow = nS - 1, ncol = nS)
for (i in 1:(nS - 1)) {
  for (j in 1:nS) {
    delta <- ifelse(test = i == j, yes = 1, no = 0)
    loc.gama <- nE + nP * (i - 1) + j  # gamma location
    e.hi[i, j] <- -delta + cof[loc.gama] / av.sh[i] + av.sh[j]
    v.hi[i, j] <- vco[loc.gama, loc.gama] / (av.sh[i] ^ 2)
  }
}
t.hi <- e.hi / sqrt(v.hi)
p.hi <- 2 * (1 - pt(q = abs(t.hi), df = df.residual(z$est)))

# D2. Hicksian elasticity: combined
hick <- data.frame(matrix(data = NA, nrow = 2 * (nS - 1), ncol= nS))
for (j in 1:nS) {
  hi <- data.frame(cbind(e.hi[, j], sqrt(v.hi)[, j], t.hi[, j], p.hi[, j])) 
  rownames(hi) <- names(av.sh)[-z$nOmit]       
  hick[, j] <- bsTab(hi)[, 2]     
}
colnames(hick) <- names(av.sh)
hicks <- cbind(Elas.Hicks = bsTab(hi)[1], hick)
hicks