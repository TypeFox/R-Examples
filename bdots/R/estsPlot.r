ests.plot <- function(part1.list) {
	coef.id1 <- part1.list$coef.id1
  coef.id2 <- part1.list$coef.id2
	coef.id3 <- part1.list$coef.id3
	coef.id4 <- part1.list$coef.id4
	groups   <- part1.list$groups
	
	coef.agg <- as.data.frame(rbind(coef.id1, coef.id2, coef.id3, coef.id4))
	if(length(names(coef.agg)) == 4) {
		names(coef.agg) <- c("mini", "peak", "slope", "cross")
		coef.agg$Group <- c(rep(groups[1], nrow(coef.id1)), rep(groups[2], nrow(coef.id2)))
		
		par(mfrow = c(2,2))
		for(i in 1:4) {
			hist(coef.agg[,i], xlab = names(coef.agg)[i], main = names(coef.agg)[i])
		}
	} else {
		names(coef.agg) <- c("mu", "height", "sig1", "sig2", "base1", "base2")
		coef.agg$Group <- c(rep(groups[1], nrow(coef.id1)), rep(groups[2], nrow(coef.id2)))
		
		par(mfrow = c(3,2))
		for(i in 1:6) {
			hist(coef.agg[,i], xlab = names(coef.agg)[i], main = names(coef.agg)[i])
		}
	}
}