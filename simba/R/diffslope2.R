"diffslope2" <- 
function(x1, y1, x2, y2, permutations=1000, resc.x=FALSE, resc.y=TRUE, ...) {
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
	do.lm <- function(df, sub) {
	    as.numeric(lm(df[sample(nrow(df), nrow(sub)),])$coefficients[2])
	}   
	m1 <- data.frame(y1, x1)
	m2 <- data.frame(y2, x2)
	names(m1) <- c("x","y")
	names(m2) <- c("x","y")
	m1.lm <- lm(m1) #die beiden linearen modelle rechnen
	m2.lm <- lm(m2)
	ds0 <- as.numeric(m1.lm$coefficients[2]-m2.lm$coefficients[2]) #deren differenz ausrechnen
	m12 <- rbind(m1, m2)
	perm1 <- sapply(1:permutations, function(x) do.lm(m12, m1))
    perm2 <- sapply(1:permutations, function(x) do.lm(m12, m1))
    perms <- perm1-perm2
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
	res
}