## ---- results="hide", message=FALSE--------------------------------------
## load "MeanShift" package
library( MeanShift )

## load `seeds` dataset
load( "seeds.RData" )

## ---- message=FALSE------------------------------------------------------
## wheat variety labels
seeds.labels <- seeds[,"variety"]

## organize data by columns
seeds.data <- t( seeds[,c( "area", "perimeter", "compactness", 
                      "length", "width", "asymmetry", 
                      "groove.length" )] )

print( dim( seeds.data ) )

## standardize the variables
seeds.data <- seeds.data / apply( seeds.data, 1, sd )

## form a set of candidate bandwidths
h.cand <- quantile( dist( t( seeds.data ) ), seq( 0.05, 0.40, by=0.05 ) )

## ---- message=FALSE------------------------------------------------------
## perform mean shift clustering with the blurring version of the algorithm
system.time( bms.clustering <- lapply( h.cand,
function( h ){ bmsClustering( seeds.data, h=h ) } ) )

## ---- echo=FALSE, results="hide"-----------------------------------------
tmp.labels3 <- bms.clustering[[3]]$labels
tmp.labels3[tmp.labels3==3] <- "pink"
tmp.labels3[tmp.labels3==4] <- 3
tmp.labels3[tmp.labels3=="pink"] <- 4
bms.clustering[[3]]$labels <- as.integer( tmp.labels3 )

bms.clustering[[3]]$components <- bms.clustering[[3]]$components[,c( 1, 2, 4, 3, 5 )]
colnames( bms.clustering[[3]]$components ) <- colnames( bms.clustering[[3]]$components )[c( 1, 2, 4, 3, 5 )]

## ----fig.width=7, fig.height=4, fig.align="center"-----------------------
## the resulting object is a list with names "components" and "labels"
class( bms.clustering[[1]] )
names( bms.clustering[[1]] )

## extract the cluster labels
ms.labels1 <- bms.clustering[[1]]$labels
print( ms.labels1 )

## extract the cluster modes/representatives
ms.modes1 <- bms.clustering[[1]]$components
print( ms.modes1 )

## plot
par( mfrow=c( 1, 2 ) )
plot( seeds.data[5,], seeds.data[6,], col=bms.clustering[[1]]$labels,
xlab=names( seeds )[5], ylab=names( seeds )[6], main="Mean shift labels",
cex=0.65, pch=16 )
plot( seeds.data[5,], seeds.data[6,], col=seeds.labels,
xlab=names( seeds )[5], ylab=names( seeds )[6], main="True labels",
cex=0.65, pch=16 )

## bandwidth h is too small -> "overclustering"

## extract the cluster labels
ms.labels6 <- bms.clustering[[6]]$labels
print( ms.labels6 )

## extract the cluster modes/representatives
ms.modes6 <- bms.clustering[[6]]$components
print( ms.modes6 )

## plot
par( mfrow=c( 1, 2 ) )
plot( seeds.data[5,], seeds.data[6,], col=bms.clustering[[8]]$labels,
xlab=names( seeds )[5], ylab=names( seeds )[6], main="Mean shift labels",
cex=0.65, pch=16 )
plot( seeds.data[5,], seeds.data[6,], col=seeds.labels,
xlab=names( seeds )[5], ylab=names( seeds )[6], main="True labels",
cex=0.65, pch=16 )

## bandwidth h is too large -> "underclustering"

## extract the cluster labels
ms.labels3 <- bms.clustering[[3]]$labels
print( ms.labels3 )

## extract the cluster modes/representatives
ms.modes3 <- bms.clustering[[3]]$components
print( ms.modes3 )

## plot
par( mfrow=c( 1, 2 ) )
plot( seeds.data[5,], seeds.data[6,], col=bms.clustering[[3]]$labels,
xlab=names( seeds )[5], ylab=names( seeds )[6], main="Mean shift labels",
cex=0.65, pch=16 )
## add estimated cluster modes to the plot
points( ms.modes3[5,], ms.modes3[6,], col=1:ncol( ms.modes3 ),
pch="+", cex=3 )
plot( seeds.data[5,], seeds.data[6,], col=seeds.labels,
xlab=names( seeds )[5], ylab=names( seeds )[6], main="True labels",
cex=0.65, pch=16 )

## just right!


