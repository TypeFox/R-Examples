getQhz_Adaptive <- function(ARG,kType,WIN,both=T,onlyk2=F){
    location <- c(ARG[1],ARG[2])
    h <- ARG[3]
    qhz <- qhz_sq <- NA
    iden <- matrix(c(1,0,0,1),2,2)


    if(is.na(h)||any(is.na(location))) return(c(qhz=NA,qhz_sq=NA))
    if(h==0) return(c(qhz=1,qhz_sq=1))
    
    if(kType=="gaus"){
		if(onlyk2){
		qsam2 <- mvrnorm(10000,c(0,0),iden*.5*h*h)
		qhz_sq <- (1/(4*pi))*translateCoord(location,mvmat=qsam2,WIN=WIN)
		} else {
			qsam <- mvrnorm(10000,c(0,0),iden*h*h)
			qhz <- translateCoord(location,mvmat=qsam,WIN=WIN)
			if(both){
				qsam2 <- mvrnorm(10000,c(0,0),iden*.5*h*h)
				qhz_sq <- (1/(4*pi))*translateCoord(location,mvmat=qsam2,WIN=WIN)
			}
		}
    } else {
		if(onlyk2){
			qsam2 <- rpoint(n=10000,function(x,y,h) {(1/(h^4))*(9/(pi^2))*(1-sqrt((x/h)^2+(y/h)^2)^2)^2*(1-sqrt((x/h)^2+(y/h)^2)^2)^2 * (abs(sqrt(((x/h)^2+(y/h)^2)))<=1)},win=owin(c(-h,h),c(-h,h)),h=h)
			qsam2 <- matrix(c(qsam2$x,qsam2$y),length(qsam2$x),2)
			qhz_sq <- translateCoord(location,mvmat=qsam2,WIN=WIN)
		} else {
			qsam <- rpoint(n=10000,function(x,y,h) {(1/(h^2))*(3/pi)*(1-sqrt((x/h)^2+(y/h)^2)^2)^2 * (abs(sqrt(((x/h)^2+(y/h)^2)))<=1)},win=owin(c(-h,h),c(-h,h)),h=h)
			qsam <- matrix(c(qsam$x,qsam$y),length(qsam$x),2)   
			qhz <- translateCoord(location,mvmat=qsam,WIN=WIN)             
			if(both){
				qsam2 <- rpoint(n=10000,function(x,y,h) {(1/(h^4))*(9/(pi^2))*(1-sqrt((x/h)^2+(y/h)^2)^2)^2*(1-sqrt((x/h)^2+(y/h)^2)^2)^2 * (abs(sqrt(((x/h)^2+(y/h)^2)))<=1)},win=owin(c(-h,h),c(-h,h)),h=h)
				qsam2 <- matrix(c(qsam2$x,qsam2$y),length(qsam2$x),2)
				qhz_sq <- translateCoord(location,mvmat=qsam2,WIN=WIN)
			}
		}
    }
    return(c(qhz,qhz_sq))
}        

