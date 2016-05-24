"diffic" <- 
function(x1, y1, x2, y2, permutations=1000, resc.x=FALSE, resc.y=FALSE, trace=FALSE, ...) {
	if (resc.x) {
		maxS <- max(mean(x1), mean(x2)) #get the higher one of the two means
		x1 <- x1+(maxS-mean(x1)) #and apply it to both data sets to make the x the same scale
		x2 <- x2+(maxS-mean(x2))
		}
	if (resc.y) {
		maxD <- max(mean(y1), mean(y2))
		y1 <- y1+(maxD-mean(y1))
		y2 <- y2+(maxD-mean(y2))
		}
	m1 <- data.frame(as.numeric(y1), as.numeric(x1))
	m2 <- data.frame(as.numeric(y2), as.numeric(x2))
	names(m1) <- c("x","y")
	names(m2) <- c("x","y")
	m1.lm <- lm(m1) #calculate the linear models
	m2.lm <- lm(m2)
	ds0 <- as.numeric(m1.lm$coefficients[1]-m2.lm$coefficients[1]) #calculate the difference in intercept
	perms <- numeric(permutations)
	if (trace) {cat(permutations, "perms: ")}
	for(i in 1:permutations) {
		tmp1 <- sample(nrow(m1), nrow(m1)/2)
		tmp2 <- sample(nrow(m2), nrow(m2)/2)
		m12 <- rbind(m1[tmp1,], m2[tmp2,])
		m21 <- rbind(m1[-tmp1,], m2[-tmp2,])
		perms[i] <- as.numeric(lm(m12)$coefficients[1]-lm(m21)$coefficients[1])
		if (trace) {cat(paste(i,""))}
		}
	if (ds0 >= 0) {	
		signif <- length(perms[perms>=ds0])/permutations
	}
	else {
		signif <- length(perms[perms<=ds0])/permutations
	}	
	if (signif == 0) {
		signif <- 1/permutations
	}
	res <- c(call=match.call())
	res$slope.diff <- as.numeric(ds0)
	res$signif <- signif
	res$permutations <- permutations
	res$perms <- perms
	class(res) <- "dsl"
	res
}