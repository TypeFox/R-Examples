"iband" <- function(object)
{

if (FALSE)
{

## Defining a vector of a fine grid of concentrations
concVec<- with(ryegrass, seq(min(conc), max(conc), length.out=150))

## Calculating predicted values including confidence intervals
predictMatrix<-predict(ryegrass.m1, newdata = data.frame(conc = concVec),
interval="confidence")

## Adding confidence limits to the plot
plot(ryegrass.m1, broken = TRUE)
lines(concVec, predictMatrix[, 2], lty = 2)
lines(concVec, predictMatrix[, 3], lty = 2)
}
}