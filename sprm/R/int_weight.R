# TODO: Add comment
# 
# Author: hoffmann
###############################################################################


int_weight <- function(dt1, dt2, ind1, ind2, y0, fun, probp1, hampelp2, hampelp3, probp4=NULL, yweights=FALSE, center1, cov1, center2, cov2){
	
	a <- ncol(dt1)
	n <- nrow(dt1) + nrow(dt2)
	
    # cov1 <- covMcd(dt1)
	wtn1 <- mahalanobis(dt1, center1, cov1)
    wtn1 <- wtn1/median(wtn1) *qchisq(0.5, a)
    # cov2 <- covMcd(dt2)
	wtn2 <- mahalanobis(dt2, center2, cov2)
    wtn2 <- wtn2/median(wtn2) *qchisq(0.5, a)

	
	probct <- qchisq(probp1, a)
	if(fun=="Fair"){
		wte1 <- 1/(1 + abs(wtn1/(probct*2)))
		wte2 <- 1/(1 + abs(wtn2/(probct*2)))
	} else if(fun=="Huber") {
		wte1 <- wtn1
		wte1[which(wtn1 <= probct)] <- 1
		wte1[which(wtn1 > probct)] <- probct/abs(wtn1[which(wtn1 > probct)])
		wte2 <- wtn2
		wte2[which(wtn2 <= probct)] <- 1
		wte2[which(wtn2 > probct)] <- probct/abs(wtn2[which(wtn2 > probct)])
	} else if(fun=="Hampel") {
		
		probct <- qchisq(probp1, a)
		hampelb <- qchisq(hampelp2, a)
		hampelr <- qchisq(hampelp3, a)
		
		wte1 <- wtn1
		wte2 <- wtn2
		
		wte1[which(wtn1 <= probct)] <- 1 
		wte1[which(wtn1 > probct & wtn1 <= hampelb)] <- probct/abs(wtn1[which(wtn1 > probct & wtn1 <= hampelb)])
		wte1[which(wtn1 > hampelb & wtn1 <= hampelr)] <- probct*(hampelr-abs(wtn1[which(wtn1 > hampelb & wtn1 <= hampelr)]))/(hampelr -hampelb)*1/abs(wtn1[which(wtn1 > hampelb & wtn1 <= hampelr)])
		wte1[which(wtn1 > hampelr)] <- 0
		wte2[which(wtn2 <= probct)] <- 1 
		wte2[which(wtn2 > probct & wtn2 <= hampelb)] <- probct/abs(wtn2[which(wtn2 > probct & wtn2 <= hampelb)])
		wte2[which(wtn2 > hampelb & wtn2 <= hampelr)] <- probct*(hampelr-abs(wtn2[which(wtn2 > hampelb & wtn2 <= hampelr)]))/(hampelr -hampelb)*1/abs(wtn2[which(wtn2 > hampelb & wtn2 <= hampelr)])
		wte2[which(wtn2 > hampelr)] <- 0
	}
	
	wte <- vector(length=n)
	wte[ind1] <- wte1
	wte[ind2] <- wte2
	
	if (yweights!=FALSE){
		
		Ts <- matrix(ncol=1, nrow=n)
		Ts[ind1,] <- scale(dt1[,1], center=center1[1], scale=sqrt(cov1[1,1])) + center1[1]
		Ts[ind2,] <- scale(dt2[,1], center=center2[1], scale=sqrt(cov2[1,1])) + center2[1]
		
		m1 <- center1[1]
		m2 <- center2[1]
		m <- (m1+m2)/2
		
		ty <- (Ts-m)*y0
		wye <- apply(as.matrix(ty), 1, function(x){biweight(x,M=0,c=min(0,qnorm(probp4, mean=median(ty), sd=1)))})
		
	} else {
		wye <- wte
	}
	
	we <- sqrt(wte * wye)
	if(any(we<1e-6)){
		w0 <- which(we<1e-6)
		we <- replace(we,list=w0,values=1e-6)
	}
	
	return(list(we=we, wte=wte, wye=wye))
}
