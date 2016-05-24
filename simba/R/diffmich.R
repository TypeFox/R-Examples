"diffmich" <- 
function(x1, y1, x2, y2, permutations=1000, a=3, b=0.5, trace=FALSE, ...) {
    xzus <- as.numeric(c(x1, x2))
    yzus <- as.numeric(c(y1, y2))
	##die beiden zu vergleichenden modelle rechnen und differenz bilden
	dma0 <- as.numeric(coef(fitmich(x1, y1, a=a, b=b))[1]) - as.numeric(coef(fitmich(x2, y2, a=a, b=b))[1])
	dmb0 <- as.numeric(coef(fitmich(x1, y1, a=a, b=b))[2]) - as.numeric(coef(fitmich(x2, y2, a=a, b=b))[2])
	permsa <- as.numeric(1:permutations)
	permsb <- as.numeric(1:permutations)
	if (trace) {cat(permutations, "perms: ")}
	for(i in 1:permutations) {
		tmp1 <- sample(length(xzus), length(x1))
		tmp2 <- sample(length(xzus), length(x2))
		xt1 <- xzus[tmp1]
		yt1 <- yzus[tmp1]
		xt2 <- xzus[tmp2]
		yt2 <- yzus[tmp2]
		permsa[i] <- as.numeric(coef(fitmich(xt1, yt1, a=a, b=b))[1]) - as.numeric(coef(fitmich(xt2, yt2, a=a, b=b))[1])
		permsb[i] <- as.numeric(coef(fitmich(xt1, yt1, a=a, b=b))[2]) - as.numeric(coef(fitmich(xt2, yt2, a=a, b=b))[2])
		if (trace) {cat(paste(i,""))}
		}
	if (dma0 >= 0) {	
		signifa <- sum(permsa >= dma0)/permutations
	}
	else {
		signifa <- sum(permsa <= dma0)/permutations
	}	
	if (signifa == 0) {
		signif <- 1/permutations
	}
	if (dmb0 >= 0) {	
		signifb <- sum(permsb >= dmb0)/permutations
	}
	else {
		signifb <- sum(permsb <= dmb0)/permutations
	}
	if (signifa == 0) {
		signif <- 1/permutations
	}
	res <- c(call=match.call())
	res$diffa <- as.numeric(dma0)
	res$signifa <- signifa
	res$diffb <- as.numeric(dmb0)
	res$signifb <- signifb
	res$permutations <- permutations
	res$permsa <- permsa
	res$permsb <- permsb
	class(res) <-"diffmich"
	res
}
