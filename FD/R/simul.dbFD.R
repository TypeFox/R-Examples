`simul.dbFD` <-
function(s = c(5, 10, 15, 20, 25, 30, 35, 40), t = 3, r = 10, p = 100, tr.method = c("unif", "norm", "lnorm"), abun.method = c("lnorm", "norm", "unif"), w.abun = TRUE){

if (p < max(s) ) stop("'p' must be greater than the maximum species richness level in 's'.")

l.s<-length(s)
r.rep<-rep(r,l.s)
nb.sp<-rep(s,r.rep)

tr.method <- match.arg(tr.method)
abun.method <- match.arg(abun.method)

# create pool of species (overall trait matrix)
traits <- matrix(NA, p, t)
if (tr.method == "unif") traits <- apply(traits, 2, function(p) runif(p))
if (tr.method == "norm") traits <- apply(traits, 2, function(p) rnorm(p))
if (tr.method == "lnorm") traits <- apply(traits, 2, function(p) rlnorm(p))

abun <- list(rep(0, p) )
abun <- rep(abun, r *l.s)

fill.abun <- function(x, y, z){
  set <- sample(1:length(x), size = y)
  if (z == "lnorm") x[set] <- rlnorm(length(set))
  if (z == "norm") x[set] <- rnorm(length(set))	
  if (z == "unif") x[set] <- runif(length(set))
  return(x)
}

abun <- mapply(fill.abun, abun, nb.sp, MoreArgs = list(z = abun.method) )
abun <- t(abun)

# names traits
no.tr <- 1:t
tr <- c("tr")
names.tr <- paste(tr, no.tr, sep = "")

#names species
no.sp <- 1:p
sp <- c("sp")
names.sp <- paste(sp, no.sp, sep = "")

dimnames(traits) <- list(names.sp, names.tr)

# names abundances
no.com <- 1:(r * l.s)
names.com <- paste("com", no.com, sep = "")
dimnames(abun) <- list(names.com, names.sp)

values <- dbFD(traits, abun, calc.FDiv = T, w.abun = w.abun, messages = F, calc.CWM = F)
results <- cbind(nb.sp, values$FDis, values$FRic, values$FEve, values$FDiv, values$RaoQ)
names.var <- c("nb.sp", "FDis", "FRic", "FEve", "FDiv", "RaoQ")
dimnames(results) <- list(names.com, names.var)

# calculate FDis of gamma to explore concavity criterion
abun.gamma <- apply(abun, 2, mean)
dbFD.gamma <- dbFD(traits, abun.gamma, calc.FDiv = F, w.abun = w.abun, messages = F)
FDis.gamma <- dbFD.gamma$FDis

# average FDis of all communities
FDis.mean <- mean(values$FDis)

# correlations and tests
nbsp.FDis <- cor.test(results[,1], results[,2])
cor1 <- bquote(italic(r)==.(round(nbsp.FDis$estimate, 3)))
FRic.FDis <- cor.test(results[,2], results[,3])
cor2 <- bquote(italic(r)==.(round(FRic.FDis$estimate, 3)))
FEve.FDis <- cor.test(results[, 2], results[, 4])
cor3 <- bquote(italic(r)==.(round(FEve.FDis$estimate, 3)))
FDiv.FDis <- cor.test(results[, 2], results[, 5])
cor4 <- bquote(italic(r)==.(round(FDiv.FDis$estimate, 3)))
RaoQ.FDis <- cor.test(results[, 6], results[, 2])
cor6 <- bquote(italic(r)==.(round(RaoQ.FDis$estimate, 3)))
nbsp.RaoQ <- cor.test(results[, 1], results[, 6])
cor5 <- bquote(italic(r)==.(round(nbsp.RaoQ$estimate, 3)))

# plot results

par(mar = c(5, 5, 4, 2), mfrow = c(3, 2),las = 1)

plot(results[, 3], results[, 2], xlab = "", ylab = "", pch = ".", cex.axis = 1.2)
title(xlab = "FRic", ylab = "FDis",cex.lab = 1.5)
mtext("a", line = 1, adj = 0, cex = 1.5, font = 2)
mtext(cor2, side = 3, line = 1, cex = 1.2)

plot(results[, 5], results[, 2], xlab = "", ylab = "", pch = ".", cex.axis = 1.2)
mtext("b", line = 1, adj = 0, cex = 1.5, font = 2)
title(xlab = "FDiv", ylab = "FDis", cex.lab = 1.5)
mtext(cor4, side = 3, line = 1, cex = 1.2)

plot(results[, 4], results[, 2], xlab = "", ylab = "", pch = ".", cex.axis = 1.2)
title(xlab = "FEve", ylab = "FDis", cex.lab = 1.5)
mtext("c", line = 1, adj = 0, cex = 1.5, font = 2)
mtext(cor3, side = 3, line = 1, cex = 1.2)

plot(results[, 6], results[, 2], xlab = "", ylab = "", pch = ".", cex.axis = 1.2)
title(xlab = "Rao's Q", ylab = "FDis", cex.lab = 1.5)
mtext("d", line = 1, adj = 0, cex = 1.5, font = 2)
mtext(cor6, side = 3, line = 1, cex = 1.2)

boxplot(results[, 2]~results[, 1], cex.axis = 1.2)
mtext("e", line = 1, adj = 0, cex = 1.5, font = 2)
title(xlab = "Species richness", ylab = "FDis", cex.lab = 1.5)
mtext(cor1, side = 3, line = 1, cex = 1.2)

boxplot(results[, 6]~results[, 1], cex.axis = 1.2)
mtext("f", line = 1, adj = 0, cex = 1.5, font = 2)
title(xlab = "Species richness", ylab = "Rao's Q", cex.lab = 1.5)
mtext(cor5, side = 3, line = 1, cex = 1.2)

res <- list()
res$results <- results
res$traits <- traits
res$abun <- abun
res$abun.gamma <- abun.gamma
res$FDis.gamma <- FDis.gamma
res$FDis.mean <- FDis.mean

invisible(res)
}
