getQhz_Fixed <- function(Xseq,Yseq,kType,WIN,h,both=T,onlyk2=F){
    
    gridLen <- length(Xseq)
    gridCoords <- matrix(c(Xseq,Yseq),gridLen,2)
    
    if(is.na(h[1])) return(list(qhz=NA,qhz_sq=NA))

    
    qhz <- qhz_sq <- NA
    
    if(kType=="gaus"){
        if(onlyk2){
			qsam2 <- mvrnorm(10000,c(0,0),diag(2)*.5*h[1]*h[1])
			qhz_sq <- (1/(4*pi))*apply(gridCoords,1,translateCoord,mvmat=qsam2,WIN=WIN)
		} else {
			iden <- matrix(c(1,0,0,1),2,2)
			qsam <- mvrnorm(10000,c(0,0),iden*h[1]*h[1])
			if(both){
				qsam2 <- mvrnorm(10000,c(0,0),iden*.5*h[1]*h[1])
				qhz_sq <- (1/(4*pi))*apply(gridCoords,1,translateCoord,mvmat=qsam2,WIN=WIN)
			}
			qhz <- apply(gridCoords,1,translateCoord,mvmat=qsam,WIN=WIN)
		}
    } else {
		if(onlyk2){
			qsam2 <- rpoint(n=10000,function(x,y,h) {(1/(h^4))*(9/(pi^2))*(1-sqrt((x/h)^2+(y/h)^2)^2)^2*(1-sqrt((x/h)^2+(y/h)^2)^2)^2 * (abs(sqrt(((x/h)^2+(y/h)^2)))<=1)},win=owin(c(-h,h),c(-h,h)),h=h[1])
			qsam2 <- matrix(c(qsam2$x,qsam2$y),length(qsam2$x),2)
			qhz_sq <- apply(gridCoords,1,translateCoord,mvmat=qsam2,WIN=WIN) 
		} else {
			qsam <- rpoint(n=10000,function(x,y,h) {(1/(h^2))*(3/pi)*(1-sqrt((x/h)^2+(y/h)^2)^2)^2 * (abs(sqrt(((x/h)^2+(y/h)^2)))<=1)},win=owin(c(-h,h),c(-h,h)),h=h[1])
			qsam <- matrix(c(qsam$x,qsam$y),length(qsam$x),2)
			if(both){
				qsam2 <- rpoint(n=10000,function(x,y,h) {(1/(h^4))*(9/(pi^2))*(1-sqrt((x/h)^2+(y/h)^2)^2)^2*(1-sqrt((x/h)^2+(y/h)^2)^2)^2 * (abs(sqrt(((x/h)^2+(y/h)^2)))<=1)},win=owin(c(-h,h),c(-h,h)),h=h[1])
				qsam2 <- matrix(c(qsam2$x,qsam2$y),length(qsam2$x),2)
				qhz_sq <- apply(gridCoords,1,translateCoord,mvmat=qsam2,WIN=WIN)        
			}
			qhz <- apply(gridCoords,1,translateCoord,mvmat=qsam,WIN=WIN)
		}
    }
    return(list(qhz=qhz,qhz_sq=qhz_sq))
}
