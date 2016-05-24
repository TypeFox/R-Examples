## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
# Subset example data to one sample
library(shazam)
db <- subset(InfluenzaDb, BARCODE == "RL013")

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Use hs1f model and normalize by junction length
dist_hs1f <- distToNearest(db, model="hs1f", first=FALSE, normalize="length", 
                           nproc=1)

# Use genotyped V assignments and 5-mer model
dist_hs5f <- distToNearest(db, vCallColumn="V_CALL_GENOTYPED", model="hs5f", 
                           first=FALSE, normalize="none", nproc=1)

## ---- eval=TRUE, warning=FALSE, fig.width=7------------------------------
# Generate m1n histogram
library(ggplot2)
p1 <- ggplot() + theme_bw() + 
    ggtitle("Distance to nearest: hs1f") + xlab("distance") +
    geom_histogram(data=dist_hs1f, aes(x=DIST_NEAREST), binwidth=0.025, 
                   fill="steelblue", color="white") +
    geom_rug(aes(x=seq(0.0, 0.5, 0.1)), color="firebrick")
plot(p1)

## ---- eval=TRUE, warning=FALSE, fig.width=7------------------------------
# Generate hs5f histogram
p2 <- ggplot() + theme_bw() + 
    ggtitle("Distance to nearest: hs5f") + xlab("distance") +
    geom_histogram(data=dist_hs5f, aes(x=DIST_NEAREST), binwidth=1, 
                   fill="steelblue", color="white") +
    geom_rug(aes(x=seq(0, 15, 1)), color="firebrick")
plot(p2)

# Zoom in to find threshold
p3 <- p2 + xlim(c(0, 15)) + scale_y_sqrt()
plot(p3)

