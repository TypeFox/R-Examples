d <- generateFakeData()

# compute small area estimates
sae <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, method="hybrid", keep.data=TRUE)

# compute unit weights
w <- uweights(sae, forTotal=TRUE)
summary(w)  # summary statistics
plot(w)  # histogram of weights
# checks
all.equal(sum(w * sae$y), sum(EST(sae) * sae$Narea))
all.equal(colSums(w * as.matrix(sae$X)), colSums(sae$Xp * sae$Narea))
