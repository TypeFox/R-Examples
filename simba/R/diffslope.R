"diffslope" <- 
function(x1, y1, x2, y2, permutations=1000, ic=FALSE, resc.x=FALSE, resc.y=TRUE, trace=FALSE, ...) {
	if (resc.x) {
		maxS <- max(mean(x1), mean(x2)) #das hoehere der beiden means herausfinden
		x1 <- x1+(maxS-mean(x1)) #und auf beide datensaetze anwenden
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
	m1.lm <- lm(m1) #die beiden linearen modelle rechnen
	m2.lm <- lm(m2)
	ds0 <- as.numeric(m1.lm$coefficients[2]-m2.lm$coefficients[2]) #deren differenz ausrechnen
	if(ic){
		m12.lmcoeff <- matrix(data=NA, nrow=permutations, ncol=2)
		m21.lmcoeff <- matrix(data=NA, nrow=permutations, ncol=2)
		dic <- as.numeric(m1.lm$coefficients[1]-m2.lm$coefficients[1])
		if (trace) {cat(permutations, "perms: ")}
		for(i in 1:permutations) {
			tmp1 <- sample(nrow(m1), nrow(m1)/2)
			tmp2 <- sample(nrow(m2), nrow(m2)/2)
			m12 <- rbind(m1[tmp1,], m2[tmp2,])
			m21 <- rbind(m1[-tmp1,], m2[-tmp2,])
			m12.lmcoeff[i,] <- as.numeric(lm(m12)$coefficients)
			m21.lmcoeff[i,] <- as.numeric(lm(m21)$coefficients)
			if (trace) {cat(paste(i,""))}
			}
		perms <- m12.lmcoeff - m21.lmcoeff
		if (ds0 >= 0) {	
			signif <- length(perms[perms[,2]>=ds0,2])/permutations
			}
		else {
			signif <- length(perms[perms[,2]<=ds0,2])/permutations
			}
		if (dic >= 0) {	
			signific <- length(perms[perms[,1]>=ds0,1])/permutations
			}
		else {
			signific <- length(perms[perms[,1]<=ds0,1])/permutations
			}
			
		if (signif == 0) {
			signif <- 1/permutations
			}
		if (signific == 0) {
			signific <- 1/permutations
			}			
		}
	else{
		perms <- vector("numeric", permutations)
		if (trace) {cat(permutations, "perms: ")}
		for(i in 1:permutations) {
			tmp1 <- sample(nrow(m1), nrow(m1)/2)
			tmp2 <- sample(nrow(m2), nrow(m2)/2)
			m12 <- rbind(m1[tmp1,], m2[tmp2,])
			m21 <- rbind(m1[-tmp1,], m2[-tmp2,])
			perms[i] <- as.numeric(lm(m12)$coefficients[2]-lm(m21)$coefficients[2])
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
		perms <- cbind(perms,perms)	
		}
	res <- c(call=match.call())
	res$slope.diff <- as.numeric(ds0)
	res$signif <- signif
	res$permutations <- permutations
	res$perms <- perms[,2]
	class(res) <- "dsl"
	if(ic) {
		res$intercept <- as.numeric(dic)
		res$signific <- as.numeric(signific)
		res$permsic <- as.numeric(perms[,1])
		class(res) <- "dsl2"
		}
	res
}