invert.ci <-
function(Y,Z,prob,perms,targetp) {

	ate <- estate(Y,Z,prob=prob)
	Ys <- genouts(Y,Z,ate)
	distro <- gendist(Ys,perms,prob=prob)
	
	mindistro <- quantile(distro,mean(c(targetp,0)))
	maxdistro <- quantile(distro,mean(c(targetp,1)))
		
	bound <- ATEg <- ATEgorig <- quantile(distro,targetp)
	bw <- min(abs(mindistro-ATEg),abs(maxdistro-ATEg))
	
	Ys1 <- genouts(Y-ATEg*Z,Z,0)
	testS <- estate(Y-ATEg*Z,Z,prob=prob)
	dist1 <- gendist(Ys1,perms,prob=prob)
	pguess <- mean(dist1 >= testS)

if (pguess >= targetp) bound <- ATEg - bw
if (pguess < targetp) bound <- ATEg + bw


# see if bound is good enough; might need to go farther

	YsM <- genouts(Y-bound*Z,Z,0)
	testM <- estate(Y-bound*Z,Z,prob=prob)
	distM <- gendist(YsM,perms,prob=prob)
	pguessM <- mean(distM >= testM)

counter.max <- 100
counter <- 0

while (pguess > targetp & pguessM > targetp) {
	temp <- ATEg
	ATEg <- bound
	bound <- ATEg - bw

	YsM <- genouts(Y-bound*Z,Z,0)
	testM <- estate(Y-bound*Z,Z,prob=prob)
	distM <- gendist(YsM,perms,prob=prob)
	pguessM <- mean(distM >= testM)	
	counter <- counter + 1
	if (counter >= counter.max) stop("Cannot Reach p.")
	}
	

while (pguess < targetp & pguessM < targetp) {
	temp <- ATEg
	ATEg <- bound
	bound <- ATEg + bw

	YsM <- genouts(Y-bound*Z,Z,0)
	testM <- estate(Y-bound*Z,Z,prob=prob)
	distM <- gendist(YsM,perms,prob=prob)
	pguessM <- mean(distM >= testM)
	counter <- counter + 1
	if (counter >= counter.max) stop("Cannot Reach p.")
}


findroot <- function(ATEg,targetp) {
	Ys1 <- genouts(Y-ATEg*Z,Z,0)
	testS <- estate(Y-ATEg*Z,Z,prob=prob)
	dist1 <- gendist(Ys1,perms,prob=prob)
	return(mean(dist1 >= testS) - targetp)
	}
	
if (pguessM == targetp) {
	ATEg <- bound
	pguess <- targetp
	}
	
if (pguess != targetp) {
	lowint <- uniroot(findroot,c(bound,ATEg),targetp=targetp)
	lowintM <- lowint$root
} else lowintM <- ATEg

return(lowintM)
}
