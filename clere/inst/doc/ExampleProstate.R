library(clere)
library(lasso2)
data(Prostate)
y <- Prostate[,"lpsa"]
x <- as.matrix(Prostate[,-which(colnames(Prostate)=="lpsa")])
itraining <- 1:(0.8*nrow(x))

## Training dataset
xt <- x[ itraining,] ; yt <- y[ itraining]

## Validation dataset
xv <- x[-itraining,] ; yv <- y[-itraining]

## Run analysis
## if you run this script on computer having more than 3 processors
## you want to take nstart > 5 otherwise we suggest that you leave it to 5
Seed <- 1234
mod  <- fitClere(y=yt,x=xt,g=5,analysis="aic",nstart=5,parallel=TRUE,seed=Seed,
                 sparse=TRUE,nItEM=2000,nBurn=1000,nItMC=10,dp=5,nsamp=1000)

print( summary(mod) )
pdf("Figure1.pdf")
plot(mod)
dev.off()

clusters(mod,thresold=0.7)
## lcavol lweight     age    lbph     svi     lcp gleason   pgg45 
##      2       2       1       1       1       1       1       1 

## Posterior probabilities of membership are available through the slot \code{P} in object of class \code{Clere}.
print( mod@P )

## Prediction Error
error <- mean( (yv - predict(mod,xv))^2 )
print(paste("Prediction error =",error))

