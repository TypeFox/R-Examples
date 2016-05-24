benchB3 <- function(method, prior=rep(1/4, 4), sv="4", scale=FALSE, ...)
{
data("B3", package = "klaR", envir = environment())
B3 <- get("B3")
complete <- FALSE
y <- B3$PHASEN
zyklus <- numeric(length(y))
zyklus[1] <- as.numeric(complete)
for (i in seq(along = y)[-1])
    {
    if ((y[i]==sv) && (y[(i-1)]!=sv))
        zyklus[i] <- zyklus[i-1] + 1
    else
        zyklus[i] <- zyklus[i-1]
    }
if(!complete) zyklus[zyklus == max(zyklus)] <- 0
zykNR <- zyklus

inCYobs <- (zykNR!=0) # obs in complete cycles
b3 <- B3[inCYobs,]
zykNR <- zykNR[inCYobs]

if (scale) b3 <- data.frame(PHASEN=b3$PHASEN, scale(b3[,-1]))

Nzyk <- max(zykNR)    # Number of cycles in Data
APE <- numeric(Nzyk)  # Put missclassification on test data here
MODEL <- list()       # Put the models here

for (i in 1:Nzyk)
    {
    # Training- and Testdata, new cycles
    trainobs <- (zykNR!=i)
    testobs <-  (zykNR==i)
    traindata <- b3[trainobs,]
    testdata <- b3[testobs,]
    
    # Testing
    if (!is.null(prior)) 
        MODEL[[i]] <- do.call(method, list(formula=PHASEN~., data=traindata, prior=prior, ...))
    else 
        MODEL[[i]] <- do.call(method, list(formula=PHASEN~., data=traindata, ...))
    pre <- predict(MODEL[[i]], testdata, ...)
    if (is.list(pre)) pre <- pre$class
    if (is.matrix(pre)) pre <- factor(max.col(pre), levels=levels(y))
    APE[i] <- mean(pre!=testdata$PHASEN)    
    cat("\nError Rate in", i, "th cycle: ", round(APE[i],3))
    }
cat("\n------------------------------------------\n")
names(APE) <- seq(along = APE)
cat("Mean Error Rate of method",method,":",round(mean(APE),3),"\n")
return(list(MODEL=MODEL, error=APE, l1co.error=mean(APE)))
}
