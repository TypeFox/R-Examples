stats <- function(x) {
	c( mean = mean(x),
	   SS = sum( (x - mean(x))^2 )
	 )}
summary(coag~diet, data=coagulation, fun=stats) -> s; s
s <- unclass(s)  # now we can access as a matrix
grandMean <- mean(coagulation$coag); grandMean
groupMean <- s[coagulation$diet,2]; groupMean
SST <- sum((coagulation$coag - grandMean)^2); SST  # total variation
SSE <- sum((coagulation$coag - groupMean)^2); SSE  # w/in group variation
SSM <- sum((groupMean - grandMean)^2 ); SSM        # b/w group variation
