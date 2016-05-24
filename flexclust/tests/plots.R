library("flexclust")

pdf("plots.pdf")
p05 <- bundestag(2005)

c6 <- cclust(p05, k=6, save.data=TRUE)

plot(c6)
plot(c6, hull="ell")

image(c6)

barplot(c6)
barplot(c6, bycluster=FALSE)
barplot(c6, oneplot=FALSE)

barchart(c6)
barchart(c6, shade=TRUE)

flexclust:::stripes(c6)
flexclust:::stripes(c6, type="second")
flexclust:::stripes(c6, type="all")
flexclust:::stripes(c6, type="all", beside=TRUE)

plot(flexclust:::shadow(c6))
plot(flexclust:::Silhouette(c6))

shadowStars(c6)
shadowStars(c6, varwidth=TRUE)
shadowStars(c6, varwidth=TRUE,
            panel=flexclust:::panelShadowSkeleton)
shadowStars(c6, varwidth=TRUE,
            panel=flexclust:::panelShadowViolin)
shadowStars(c6, varwidth=TRUE,
            panel=flexclust:::panelShadowBP)

dev.off()
